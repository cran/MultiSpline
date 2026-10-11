# MultiSpline 0.3.0

## Corrections (results from 0.2.0 may change)

* `nl_predict()`: standard errors and intervals for `lmer`/`glmer` models were
  computed from a spline basis rebuilt on the prediction grid instead of the
  basis used in estimation. They are now computed with the estimation basis.
  Point predictions were not affected.
* `nl_derivatives()`: interval estimates treated neighbouring predictions as
  independent, which made slope and curvature intervals far too wide (by a
  factor of 10 to 100 in tests). Derivatives are now computed from the model
  with the full covariance matrix. The function now takes the fitted model:
  `nl_derivatives(fit)`. Passing a prediction table still returns point
  estimates, with a warning and no intervals.
* Automatic df selection (`df = "auto"`), `nl_knots()`, and `nl_compare()`
  used REML likelihoods for Gaussian mixed models, which are not comparable
  across different fixed effects. All comparisons now use maximum
  likelihood. The final model is still fitted by REML unless `REML = FALSE`.
* `nl_r2()` and `nl_icc()` ignored random-slope variances. They now include
  them, averaged over the data (Johnson 2014); `nl_icc()` warns that the ICC
  then depends on `x`.
* `nl_fit(cluster = c("student", "school"), nested = TRUE)` built
  `(1 | student/school)`, i.e. schools nested in students, contrary to the
  documentation, and so omitted the school-level variance. It now builds
  `(1 | school/student)`: the first variable is nested within the second.
* `nl_r2()`: complementary log-log models now use pi^2/6 as the
  distribution-specific variance.
* `nl_het()`: cluster-specific curves for non-Gaussian models were on the
  link scale while the population curve was on the response scale; both are
  now on the response scale. The output `slope_variance` (which was the SD
  of predicted values, not of slopes) is renamed `between_cluster_sd`.

## New features

* `nl_fit(random = )` accepts any `lme4` random-effects specification, e.g.
  `"(1 | school/student)"` or `"(1 | idcode) + (1 | indcode)"`.
* `nl_predict(type = "marginal")` gives population-averaged curves for
  GLMMs; `scale = "link"` gives the linear predictor. The documentation now
  states exactly what is predicted.
* `controls` accepts expressions such as `"factor(race)"` or `"I(age^2)"`.
* `nl_wald()`: joint Wald test that all coefficients involving `x` are zero.
* `nl_fit()` warns when `x` has too few distinct values for the requested df.
* The data used in estimation are stored in the fit, so `nl_compare()` and
  `nl_het()` no longer need to find the original data frame.
* `nl_plot(type = "combo")` returns one combined figure (via patchwork when
  installed) instead of a list, so it can be saved with `ggsave()`.
* `nl_plot()` and `nl_turning_points()` read the variable names from their
  input, so `x` no longer has to be repeated.

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