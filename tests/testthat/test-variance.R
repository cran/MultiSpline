test_that("df = 'auto' and comparisons use ML", {
  d   <- make_data()
  fit <- suppressMessages(nl_fit(d, "y", "x", cluster = "id", df = "auto",
                                 df_range = 2:5))
  ml  <- vapply(2:5, function(k) stats::AIC(lme4::lmer(
    y ~ splines::ns(x, df = k) + (1 | id), d, REML = FALSE)), numeric(1))
  expect_equal(fit$df_search$criterion, ml, tolerance = 1e-6)
  expect_true(lme4::isREML(fit$model))
  invisible(utils::capture.output(cmp <- suppressMessages(
    nl_compare(fit, polynomial_degrees = integer(0)))))
  expect_equal(cmp$table$LogLik[cmp$table$Model == "Spline"],
               round(as.numeric(stats::logLik(lme4::refitML(fit$model))), 3))
})

test_that("ICC matches the latent-scale formula used by Stata", {
  d  <- make_data()
  fb <- nl_fit(d, "yb", "x", cluster = "id", df = 3, family = stats::binomial())
  s2 <- as.data.frame(lme4::VarCorr(fb$model))$vcov[1]
  expect_equal(unname(nl_icc(fb)[1]), s2 / (s2 + pi^2 / 3))
})

test_that("nested clusters: first variable is nested in the second", {
  d  <- make_data()
  f1 <- nl_fit(d, "y", "x", cluster = c("id", "sch"), nested = TRUE, df = 3)
  f2 <- nl_fit(d, "y", "x", random = "(1 | sch/id)", df = 3)
  expect_equal(as.numeric(stats::logLik(f1$model)),
               as.numeric(stats::logLik(f2$model)))
})

test_that("R2 includes random-slope variance (Johnson 2014)", {
  d  <- make_data()
  fr <- suppressWarnings(nl_fit(d, "y", "x", cluster = "id", df = 2,
                                random_slope = TRUE))
  invisible(utils::capture.output(r2 <- suppressMessages(nl_r2(fr))))
  m  <- fr$model
  S  <- lme4::VarCorr(m)$id
  Z  <- cbind(1, predict(splines::ns(d$x, df = 2), d$x))
  rv <- mean(rowSums((Z %*% S) * Z))
  sf <- stats::var(drop(lme4::getME(m, "X") %*% lme4::fixef(m)))
  expect_equal(unname(r2$r2["R2c"]), (sf + rv) / (sf + rv + stats::sigma(m)^2))
  expect_warning(nl_icc(fr), "random slopes")
})

test_that("Wald test and discrete-x warning", {
  d   <- make_data()
  fit <- nl_fit(d, "y", "x", cluster = "id", df = 4)
  w   <- nl_wald(fit)
  b   <- lme4::fixef(fit$model)[-1]
  V   <- as.matrix(stats::vcov(fit$model))[-1, -1]
  expect_equal(w$chisq, drop(t(b) %*% solve(V, b)))
  expect_equal(w$df, 4L)
  d$xd <- pmin(round(d$x / 2), 3)
  expect_warning(nl_fit(d, "y", "xd", df = 2), "distinct values")
})

test_that("nl_het returns a likelihood-ratio test", {
  d  <- make_data()
  f  <- nl_fit(d, "y", "x", cluster = "id", df = 2)
  invisible(utils::capture.output(h <- suppressWarnings(suppressMessages(
    nl_het(f, plot = FALSE)))))
  expect_s3_class(h$lrt, "data.frame")
  expect_equal(nrow(h$lrt), 2L)
  expect_true(is.finite(h$lrt$Chisq[2]))
})
