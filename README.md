# MultiSpline v2.0.0

**Spline-Based Nonlinear Modeling for Multilevel and Longitudinal Data**

> Author: Subir Hait · Michigan State University · haitsubi@msu.edu

---

## Overview

MultiSpline provides a unified interface for fitting, interpreting, and
visualising nonlinear relationships in single-level and multilevel regression
models using natural cubic splines, B-splines, or GAM smooths.

Version 2 is a major upgrade. Every item raised in the Moran correspondence
has been addressed, and the package now offers a complete **estimation →
interpretation → diagnostics** framework:

| Layer | Functions |
|-------|-----------|
| **Estimation** | `nl_fit()` · `nl_knots()` |
| **Prediction** | `nl_predict()` |
| **Interpretation** | `nl_derivatives()` · `nl_turning_points()` · `nl_plot()` |
| **Model selection** | `nl_compare()` |
| **Variance explained** | `nl_r2()` · `nl_icc()` |
| **Cluster heterogeneity** | `nl_het()` |

---

## What's New in v2

### 1. Two-way and Nested Clustering  *(Moran request)*

```r
# Cross-classified (e.g. students nested in both schools and districts)
fit <- nl_fit(data = df, y = "score", x = "SES",
              cluster = c("student_id", "school_id"),
              nested  = FALSE)

# Nested (students within schools)
fit <- nl_fit(data = df, y = "score", x = "SES",
              cluster = c("student_id", "school_id"),
              nested  = TRUE)
```

Internally this generates `(1 | student_id) + (1 | school_id)` for
cross-classified models and `(1 | student_id/school_id)` for nested models,
passed to `lme4::lmer()` / `lme4::glmer()`.

---

### 2. Confidence Intervals for Spline Curves  *(Moran request)*

All model types now return proper CI bands:

| Model | CI method |
|-------|-----------|
| `lm` | Analytical (standard `predict.lm`) |
| `gam` | Analytical (`mgcv`) |
| `lmerMod` | Delta method via fixed-effects variance–covariance matrix |
| `glmerMod` | Delta method on link scale + back-transformation (default) or parametric bootstrap (`glmer_ci = "boot"`) |

```r
pred <- nl_predict(fit, se = TRUE, level = 0.95)
nl_plot(pred_df = pred, x = "age", show_ci = TRUE)
```

---

### 3. Automatic Knot / df Selection

```r
# Explore and select df explicitly
ks <- nl_knots(data = df, y = "score", x = "age",
               df_range = 2:10, criterion = "AIC")
ks$best_df   # → e.g. 5

# Or let nl_fit do it automatically
fit <- nl_fit(data = df, y = "score", x = "age",
              df = "auto", df_criterion = "AIC")
```

A diagnostic plot of AIC/BIC vs df is produced automatically.

---

### 4. Multilevel R² Decomposition  *(Moran request — r2_mlm)*

```r
nl_r2(fit)
```

**Output (LMM example):**
```
=== MultiSpline R² Decomposition ===
Model type: LMM

  Marginal  R²m = 0.1823   (fixed effects only)
  Conditional R²c = 0.4761  (fixed + all random effects)

Variance Partition (r2_mlm-style):
     Component  Variance  Proportion
Fixed effects    1.4221      0.1823
student_id       1.8340      0.2350
school_id        0.8762      0.1124
Residual         3.6590      0.4690
```

Follows the Nakagawa-Schielzeth (2013) and Nakagawa-Johnson-Schielzeth (2017)
formulas for marginal and conditional R². The level-specific variance
partition mirrors the `r2_mlm` / Rights-Sterba (2019) approach.

---

### 5. Derivative-Based Interpretation Framework

```r
pred  <- nl_predict(fit, x_seq = seq(18, 80, length.out = 200))
deriv <- nl_derivatives(pred, x = "age")

# Visualise slope (first derivative) with CI band
nl_plot(deriv_df = deriv, x = "age", type = "slope")

# Visualise curvature (second derivative)
nl_plot(deriv_df = deriv, x = "age", type = "curvature")

# Trajectory + slope side by side
nl_plot(pred_df = pred, deriv_df = deriv, x = "age", type = "combo")
```

**Turning points and inflection regions:**

```r
tp <- nl_turning_points(deriv, x = "age")
tp$turning_points      # local maxima and minima with x location and type
tp$inflection_regions  # where concavity changes
tp$slope_regions       # contiguous increasing / decreasing stretches

# Overlay turning points on trajectory plot
nl_plot(pred_df = pred, x = "age",
        show_turning_points = TRUE,
        turning_points      = tp$turning_points)
```

---

### 6. Model Comparison Workflow

```r
nl_compare(fit)                           # default: linear + poly(2) + poly(3) + spline
nl_compare(fit, polynomial_degrees = 2:5) # custom comparators
```

**Output:**
```
=== MultiSpline Model Comparison ===

    Model    AIC    BIC  LogLik  npar  Deviance  LRT_vs_linear  LRT_p
   Linear 1842.1 1863.4  -916.1     5    1020.4             NA     NA
  Poly(2) 1830.5 1855.6  -908.2     6     992.3          15.74  <0.001
  Poly(3) 1825.3 1854.3  -903.6     7     978.1          25.00  <0.001
   Spline 1818.7 1851.5  -898.3     9     961.2          35.68  <0.001

  Best model by AIC: Spline
```

---

### 7. Cluster Heterogeneity in Nonlinear Effects

```r
nl_het(fit, n_clusters_plot = 50)
```

Plots cluster-specific predicted trajectories (BLUPs) against the
population-mean curve, and performs a likelihood-ratio test comparing
random-slope vs random-intercept models to test whether the nonlinear
effect genuinely varies across clusters.

---

### 8. B-spline Basis

```r
fit_bs <- nl_fit(data = df, y = "score", x = "age",
                 method = "bs", df = 6, bs_degree = 3)
```

---

### 9. Random Spline Slopes

Allows the nonlinear trajectory to vary across clusters:

```r
fit_rs <- nl_fit(data = df, y = "score", x = "age",
                 cluster = "id", random_slope = TRUE)
```

---

## Installation

```r
# From source
install.packages("path/to/MultiSpline_0.2.0.tar.gz", repos = NULL, type = "source")

# Dependencies (install first if needed)
install.packages(c("lme4", "mgcv", "dplyr", "ggplot2", "rlang", "splines"))
# Optional (for p-values in LMM)
install.packages("lmerTest")
```

---

## Quick-Start Example

```r
library(MultiSpline)

# 1. Select df automatically
ks <- nl_knots(data = nlme::Orthodont, y = "distance",
               x = "age", df_range = 2:8)

# 2. Fit multilevel model with two-way clustering and auto df
fit <- nl_fit(
  data     = nlme::Orthodont,
  y        = "distance",
  x        = "age",
  cluster  = "Subject",
  df       = "auto",
  controls = NULL
)
fit   # prints postestimation menu

# 3. Coefficient table
nl_summary(fit)

# 4. Multilevel R²
nl_r2(fit)

# 5. Model comparison
nl_compare(fit)

# 6. Predictions with CI
pred <- nl_predict(fit, se = TRUE)
nl_plot(pred_df = pred, x = "age", show_ci = TRUE, type = "trajectory")

# 7. Derivatives and turning points
deriv <- nl_derivatives(pred, x = "age")
tp    <- nl_turning_points(deriv, x = "age")
tp$turning_points

nl_plot(deriv_df = deriv, x = "age", type = "slope")
nl_plot(deriv_df = deriv, x = "age", type = "curvature")

# 8. Cluster heterogeneity
nl_het(fit)
```

---

## Changelog

### v0.2.0 (2026)
- **NEW** Two-way and nested clustering (`nested` argument)
- **NEW** CI for `glmerMod` via delta method and parametric bootstrap
- **NEW** Automatic df / knot selection (`df = "auto"`, `nl_knots()`)
- **NEW** Multilevel R² decomposition — marginal, conditional, level-specific (`nl_r2()`)
- **NEW** First and second derivatives with CI (`nl_derivatives()`)
- **NEW** Turning points, inflection regions, slope regions (`nl_turning_points()`)
- **NEW** Built-in model comparison workflow (`nl_compare()`)
- **NEW** Cluster heterogeneity analysis with LRT (`nl_het()`)
- **NEW** B-spline basis (`method = "bs"`, `bs_degree`)
- **NEW** Random spline slopes (`random_slope = TRUE`)
- **ENHANCED** `nl_plot()` now supports `type = "slope"`, `"curvature"`, `"combo"`
- **ENHANCED** `print.nl_fit` lists all v2 postestimation functions

### v0.1.0 (2026-02-27)
- Initial release: `nl_fit`, `nl_predict`, `nl_plot`, `nl_summary`, `nl_icc`

---

## Citation

Hait, S. (2026). *MultiSpline: Spline-Based Nonlinear Modeling for Multilevel
and Longitudinal Data* (R package version 0.2.0).
https://github.com/haitsubi/MultiSpline

---

## References

- Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
  obtaining R² from generalized linear mixed-effects models. *Methods in
  Ecology and Evolution*, 4(2), 133–142.
- Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient
  of determination R² from generalized linear mixed-effects models revisited
  and expanded. *Journal of the Royal Society Interface*, 14(134).
- Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
  multilevel models. *Psychological Methods*, 24(3), 309–338.
- Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*
  (2nd ed.). Chapman & Hall/CRC.
- Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
  mixed-effects models using lme4. *Journal of Statistical Software*, 67(1).
