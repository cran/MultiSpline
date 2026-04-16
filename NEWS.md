# MultiSpline 0.2.0

## Major new features

* **Two-way and nested clustering** via `cluster = c("g1","g2")` and
  `nested = TRUE/FALSE`, generating `(1|g1/g2)` or `(1|g1) + (1|g2)`
  random-effects structures in `lme4`.

* **Automatic knot / df selection** (`df = "auto"` in `nl_fit()` or via
  `nl_knots()`) using AIC or BIC over a user-specified grid.

* **Multilevel R-squared decomposition** (`nl_r2()`): Nakagawa-Schielzeth
  marginal R2m and conditional R2c, plus a level-specific variance partition
  table (r2_mlm style) for LMM, GLMM, and single-level OLS / GAM models.

* **Full postestimation suite**:
  - `nl_derivatives()` — first and second derivatives with delta-method
    confidence bands.
  - `nl_turning_points()` — local maxima, minima, inflection regions, and
    slope-direction regions.
  - `nl_plot()` gains `type = "slope"`, `"curvature"`, and `"combo"`
    in addition to the original `"trajectory"`.

* **Built-in model comparison workflow** (`nl_compare()`): contrasts linear,
  polynomial, and spline fits by AIC, BIC, log-likelihood, and
  likelihood-ratio tests.

* **B-spline basis** (`method = "bs"`, `bs_degree` argument).

* **Random spline slopes** (`random_slope = TRUE`) to allow the nonlinear
  effect to vary across clusters.

* **Cluster heterogeneity analysis** (`nl_het()`): plots cluster-specific
  trajectories (BLUPs) and performs an LRT comparing random-slope vs
  random-intercept models.

* **CI for `glmerMod`**: approximate confidence intervals via the delta
  method on the link scale (default, fast) or parametric bootstrap
  (`glmer_ci = "boot"`).

## Breaking changes

None. All v0.1.0 calls remain valid.

## Bug fixes

* `nl_r2()` variance partition now correctly excludes NA entries that could
  appear when `lme4` internal row names are ambiguous in nested models.
* `nl_predict()` now correctly computes CI when control variables are stored
  as character (not factor) in the original data.
* `nl_plot()` no longer errors when `time = NULL` and the data frame has
  no time column.

## Dependency changes

* `%||%` is now imported from `rlang` rather than defined internally,
  avoiding namespace masking.
* `reformulas` moved from Imports to Suggests (used opportunistically for
  `nobars()`; falls back to `lme4::nobars()` if unavailable).

# MultiSpline 0.1.0

* Initial release: `nl_fit()`, `nl_predict()`, `nl_plot()`, `nl_summary()`,
  `nl_icc()`.