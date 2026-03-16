# MultiSpline

An R package for fitting, predicting, and visualizing nonlinear
relationships in single-level, multilevel, and longitudinal regression
models using spline-based methods.

## Overview
In social and health science research, nonlinear effects are very common
in clustered or longitudinal data. `MultiSpline` provides a simple and
unified interface for estimating these effects without requiring manual
construction of spline terms or interaction structures.

## Installation
```r
# Install from CRAN
install.packages("MultiSpline")
# Or install the development version from GitHub
devtools::install_github("causalfragility-lab/MultiSpline")
```

## Core Functions
| Function | Description |
|:---------|:------------|
| `nl_fit()` | Fit a nonlinear single-level or multilevel model |
| `nl_summary()` | Tidy coefficient table |
| `nl_predict()` | Generate predictions with uncertainty |
| `nl_plot()` | Visualize nonlinear effects |
| `nl_icc()` | Compute intraclass correlations |

## Example
```r
library(MultiSpline)
# Simulate data
set.seed(42)
d <- data.frame(
  schid      = rep(1:10, each = 60),
  id         = rep(1:200, each = 3),
  TimePoint  = factor(rep(1:3, times = 200)),
  SES        = rnorm(600),
  math_score = rnorm(600, mean = 50, sd = 10)
)
# Fit nonlinear multilevel model
fit <- nl_fit(
  data    = d,
  y       = "math_score",
  x       = "SES",
  time    = "TimePoint",
  cluster = c("id", "schid"),
  method  = "ns",
  df      = 4
)
# Coefficient table
nl_summary(fit)
# Intraclass correlations
nl_icc(fit)
# Predictions and plot
pred <- nl_predict(fit)
nl_plot(pred, x = "SES", time = "TimePoint")
```

## License
MIT © Subir Hait

## Author
Subir Hait, Michigan State University

## Citation
If you use MultiSpline in your research, please cite:
```r
citation("MultiSpline")
```

Or:
```
Hait, S. (2026). MultiSpline: Spline-Based Nonlinear Modeling
for Multilevel and Longitudinal Data. R package version 0.1.0.
https://cran.r-project.org/package=MultiSpline
```
