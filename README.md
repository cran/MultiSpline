# MultiSpline

**Spline-based nonlinear modeling for single-level, multilevel, and longitudinal data**

MultiSpline fits a nonlinear effect of one focal predictor with a natural
cubic spline, B-spline, or GAM smooth, in ordinary or mixed-effects models
(`lme4`), and then answers the questions that usually follow: what does the
curve look like, where is it steep or flat, does it differ between clusters,
and is it worth the extra parameters?

## Installation

```r
install.packages("MultiSpline")
```

## Workflow

| Step | Function |
|------|----------|
| Fit | `nl_fit()`, `nl_knots()` |
| Curve and slope | `nl_predict()`, `nl_derivatives()`, `nl_turning_points()`, `nl_plot()` |
| Tests and comparison | `nl_wald()`, `nl_compare()` |
| Variance | `nl_r2()`, `nl_icc()` |
| Cluster heterogeneity | `nl_het()` |

## Example

Reaction times of 18 subjects over 10 days of sleep deprivation
(`lme4::sleepstudy`), with days nested in subjects:

```r
library(MultiSpline)
fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
              cluster = "Subject", df = 3)

nl_compare(fit, polynomial_degrees = 2)
#>    Model      AIC      BIC   LogLik npar Deviance LRT_vs_linear LRT_p
#>   Linear 1802.079 1814.850 -897.039    4 1794.079            NA  <NA>
#>  Poly(2) 1802.944 1818.909 -896.472    5 1792.944         1.135 0.287
#>   Spline 1805.030 1824.187 -896.515    6 1793.030         1.049 0.592
#>
#>   Best model by AIC: Linear
```

Here the comparison says a straight line is enough. All mixed models in the
comparison are fitted by maximum likelihood.

```r
nl_derivatives(fit, x_seq = c(0, 3, 6, 9))[, c("Days", "d1", "d1_lwr", "d1_upr")]
#>   Days    d1 d1_lwr d1_upr
#> 1    0  7.70 -0.374   15.8
#> 2    3  9.45  6.397   12.5
#> 3    6 11.70  8.652   14.8
#> 4    9 12.21  4.136   20.3

nl_icc(fit)
#>   Subject  Residual
#> 0.5877478 0.4122522
```

## Levels

`cluster` covers the common cases: one grouping variable, or two that are
crossed or nested (`cluster = c("student", "school"), nested = TRUE` gives
`(1 | school/student)`). For anything else pass the random-effects part of an
`lme4` formula directly:

```r
nl_fit(data, y = "ln_wage", x = "ttl_exp",
       random = "(1 | idcode) + (1 | ind_code)")
```

## What is predicted

`nl_predict()` evaluates the fitted model on a grid of 200 values of `x`
(1st to 99th percentile) with controls held at their mean or first level.
By default random effects are set to zero (`type = "conditional"`); for
binary and other non-Gaussian outcomes, `type = "marginal"` averages over the
random-effect distribution to give population-averaged probabilities.
Intervals use the delta method with the full fixed-effects covariance
matrix. See `?nl_predict` for details.

## Version 0.3.0

Version 0.3.0 corrects several results of 0.2.0 (prediction and derivative
intervals for mixed models, REML-based model selection, R-squared and ICC
with random slopes, and the order of nested grouping factors). See
`NEWS.md`.
