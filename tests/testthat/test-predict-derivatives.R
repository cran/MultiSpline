test_that("prediction SEs use the estimation basis", {
  d   <- make_data()
  fit <- nl_fit(d, "y", "x", cluster = "id", df = 4)
  m   <- fit$model
  B   <- splines::ns(d$x, df = 4)
  g   <- seq(fit$x_info$q01, fit$x_info$q99, length.out = 50)
  X   <- cbind(1, predict(B, g))
  V   <- as.matrix(stats::vcov(m))
  p   <- nl_predict(fit, x_seq = g)
  expect_equal(p$fit, drop(X %*% lme4::fixef(m)), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(p$se.fit, sqrt(rowSums((X %*% V) * X)), tolerance = 1e-10)
})

test_that("derivatives and their SEs match the analytic basis derivative", {
  d   <- make_data()
  fit <- nl_fit(d, "y", "x", cluster = "id", df = 4)
  m   <- fit$model
  B   <- splines::ns(d$x, df = 4)
  g   <- seq(fit$x_info$q01, fit$x_info$q99, length.out = 50)
  h   <- 1e-6
  D   <- cbind(0, (predict(B, g + h) - predict(B, g - h)) / (2 * h))
  V   <- as.matrix(stats::vcov(m))
  dd  <- nl_derivatives(fit, x_seq = g)
  expect_equal(dd$d1, drop(D %*% lme4::fixef(m)), tolerance = 1e-4)
  expect_equal(dd$d1_se, sqrt(rowSums((D %*% V) * D)), tolerance = 1e-4)
})

test_that("old table interface warns and gives no intervals", {
  d   <- make_data()
  fit <- nl_fit(d, "y", "x", df = 3)
  p   <- nl_predict(fit)
  expect_warning(dd <- nl_derivatives(p, x = "x"), "takes the fitted model")
  expect_true(all(is.na(dd$d1_lwr)))
})

test_that("lm and gam SEs agree with their own predict methods", {
  d  <- make_data()
  fl <- nl_fit(d, "y", "x", df = 3, controls = "factor(grp)")
  pl <- nl_predict(fl)
  ref <- stats::predict(fl$model, newdata = as.data.frame(pl), se.fit = TRUE)
  expect_equal(pl$se.fit, unname(ref$se.fit), tolerance = 1e-10)
  fg <- nl_fit(d, "y", "x", method = "gam")
  pg <- nl_predict(fg)
  rg <- mgcv::predict.gam(fg$model, newdata = as.data.frame(pg), se.fit = TRUE)
  expect_equal(pg$se.fit, unname(as.numeric(rg$se.fit)), tolerance = 1e-10)
})

test_that("marginal GLMM curve matches numerical integration", {
  d  <- make_data()
  fb <- nl_fit(d, "yb", "x", cluster = "id", df = 3, family = stats::binomial())
  g  <- c(1, 4, 8)
  pm <- nl_predict(fb, x_seq = g, type = "marginal")
  et <- nl_predict(fb, x_seq = g, scale = "link")$fit
  s  <- sqrt(as.data.frame(lme4::VarCorr(fb$model))$vcov[1])
  ref <- vapply(et, function(e) stats::integrate(
    function(u) stats::plogis(e + u) * stats::dnorm(u, 0, s), -Inf, Inf)$value,
    numeric(1))
  expect_equal(pm$fit, ref, tolerance = 1e-6)
  pc <- nl_predict(fb, x_seq = g)
  expect_true(all(pc$lwr >= 0 & pc$upr <= 1))
})
