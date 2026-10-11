#' Built-in model comparison workflow
#'
#' @description
#' Compares a nonlinear spline model against simpler alternatives—a linear
#' model and optional polynomial terms—to help researchers justify the
#' spline approach over simpler specifications.
#'
#' For each model, the function reports:
#' \itemize{
#'   \item AIC and BIC
#'   \item Log-likelihood (and likelihood-ratio test versus the linear model
#'     where available)
#'   \item Number of parameters (df)
#'   \item Residual variance / deviance
#' }
#'
#' Mixed models are compared on their maximum-likelihood fits: the spline
#' model is refitted by ML if it was estimated by REML, and the comparison
#' models are fitted by ML. REML likelihoods cannot be compared across
#' models with different fixed effects. The comparison models use the same
#' random-effects structure as the spline model (random intercepts only if
#' the spline model has random spline slopes).
#'
#' @param object An \code{nl_fit} object (the spline model to compare against).
#' @param polynomial_degrees Integer vector of polynomial degrees to include
#'   as additional comparators. Default \code{c(2L, 3L)}. Set to
#'   \code{integer(0)} to skip.
#' @param digits Integer; decimal places for the output table. Default \code{3}.
#' @param return_models Logical; if \code{TRUE}, also returns the fitted
#'   comparison model objects. Default \code{FALSE}.
#'
#' @return A list (invisibly) of class \code{"nl_compare"} with:
#'   \describe{
#'     \item{\code{table}}{A data frame with the comparison statistics.}
#'     \item{\code{models}}{Named list of fitted models (when
#'       \code{return_models = TRUE}).}
#'     \item{\code{best}}{Name of the model with the lowest AIC.}
#'   }
#'   The table is pretty-printed automatically.
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_knots}}
#'
#' @examples
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' nl_compare(fit, polynomial_degrees = 2)
#'
#' @export
nl_compare <- function(
    object,
    polynomial_degrees = c(2L, 3L),
    digits             = 3L,
    return_models      = FALSE
) {

  if (!inherits(object, "nl_fit"))
    stop("`object` must be an 'nl_fit' object.", call. = FALSE)


  y        <- object$y
  x        <- object$x
  controls <- object$controls
  cluster  <- object$cluster
  nested   <- object$nested
  family   <- object$family
  method   <- object$method
  df       <- object$df
  k        <- object$k
  bs_deg   <- object$bs_degree

  ctrl_part <- if (!is.null(controls)) paste(controls, collapse = " + ") else NULL

  # Helper: fit a model and extract fit stats
  .fit_stats <- function(nm, mod) {
    tryCatch({
      aic <- stats::AIC(mod)
      bic <- stats::BIC(mod)
      ll  <- as.numeric(stats::logLik(mod))
      npar<- attr(stats::logLik(mod), "df")
      dev <- if (inherits(mod, "lm")) {
        stats::deviance(mod)
      } else if (inherits(mod, "lmerMod")) {
        stats::deviance(mod)
      } else {
        tryCatch(stats::deviance(mod), error = function(e) NA_real_)
      }
      data.frame(Model = nm, AIC = aic, BIC = bic, LogLik = ll,
                 npar = npar, Deviance = dev,
                 LRT_vs_linear = NA_real_, LRT_p = NA_real_,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      message("Could not fit model '", nm, "': ", conditionMessage(e))
      NULL
    })
  }

  # --- 1. Linear model ---
  lin_rhs  <- paste(c(x, controls), collapse = " + ")
  rand_str <- if (!is.null(object$random)) object$random else
    if (!is.null(cluster)) .build_rand_formula(cluster, nested, x, df, bs_deg,
                                               method, FALSE) else NULL
  if (isTRUE(object$random_slope))
    message("Note: the comparison models use random intercepts only; the ",
            "spline model's random slopes make its fit not directly ",
            "comparable.")
  lin_form <- if (!is.null(cluster)) {
    stats::as.formula(paste0(y, " ~ ", lin_rhs, " + ", rand_str))
  } else {
    stats::as.formula(paste0(y, " ~ ", lin_rhs))
  }

  fit_data <- object$data %||% tryCatch(
    eval(object$call$data, envir = parent.frame()),
    error = function(e) NULL)
  if (is.null(fit_data))
    stop("Original data not available; refit with MultiSpline >= 0.3.0.",
         call. = FALSE)

  lin_mod <- tryCatch({
    if (!is.null(cluster)) {
      if (identical(family$family, "gaussian"))
        lme4::lmer(lin_form, data = fit_data, REML = FALSE)
      else
        lme4::glmer(lin_form, data = fit_data, family = family)
    } else if (method == "gam") {
      stats::lm(lin_form, data = fit_data)
    } else {
      stats::lm(lin_form, data = fit_data)
    }
  }, error = function(e) { message("Linear model failed: ", e$message); NULL })

  model_list <- list()
  if (!is.null(lin_mod)) model_list[["Linear"]] <- lin_mod

  # --- 2. Polynomial models ---
  if (length(polynomial_degrees) > 0L) {
    for (deg in as.integer(polynomial_degrees)) {
      poly_main <- paste0("poly(", x, ", ", deg, ")")
      poly_rhs  <- paste(c(poly_main, controls), collapse = " + ")
      poly_form <- if (!is.null(cluster)) {
        stats::as.formula(paste0(y, " ~ ", poly_rhs, " + ", rand_str))
      } else {
        stats::as.formula(paste0(y, " ~ ", poly_rhs))
      }
      nm <- paste0("Poly(", deg, ")")
      pm <- tryCatch({
        if (!is.null(cluster)) {
          if (identical(family$family, "gaussian"))
            lme4::lmer(poly_form, data = fit_data, REML = FALSE)
          else
            lme4::glmer(poly_form, data = fit_data, family = family)
        } else {
          stats::lm(poly_form, data = fit_data)
        }
      }, error = function(e) { message(nm, " failed: ", e$message); NULL })
      if (!is.null(pm)) model_list[[nm]] <- pm
    }
  }

  # --- 3. Spline model (the nl_fit object) ---
  sp <- object$model
  if (inherits(sp, "lmerMod") && lme4::isREML(sp)) sp <- lme4::refitML(sp)
  model_list[["Spline"]] <- sp

  # --- 4. Collect fit statistics ---
  rows <- lapply(names(model_list), function(nm) {
    .fit_stats(nm, model_list[[nm]])
  })
  rows <- Filter(Negate(is.null), rows)
  tbl  <- do.call(rbind, rows)

  # --- 5. LRT vs linear ---
  if (!is.null(lin_mod) && nrow(tbl) > 1L) {
    lin_ll   <- as.numeric(stats::logLik(lin_mod))
    lin_npar <- attr(stats::logLik(lin_mod), "df")
    for (i in seq_len(nrow(tbl))) {
      nm_i <- tbl$Model[i]
      if (nm_i == "Linear") next
      mod_i <- model_list[[nm_i]]
      if (is.null(mod_i)) next
      ll_i   <- tbl$LogLik[i]
      npar_i <- tbl$npar[i]
      delta_df <- npar_i - lin_npar
      if (is.finite(ll_i) && is.finite(lin_ll) && delta_df > 0) {
        lrt_stat <- 2 * (ll_i - lin_ll)
        tbl$LRT_vs_linear[i] <- round(lrt_stat, digits)
        tbl$LRT_p[i]         <- stats::pchisq(lrt_stat, df = delta_df,
                                               lower.tail = FALSE)
      }
    }
  }

  # Nicely round
  for (col in c("AIC","BIC","LogLik","Deviance","LRT_vs_linear")) {
    if (col %in% names(tbl)) tbl[[col]] <- round(tbl[[col]], digits)
  }
  tbl$LRT_p <- ifelse(is.na(tbl$LRT_p), NA_character_,
                      format.pval(tbl$LRT_p, digits = digits, eps = 0.001))

  best <- if (nrow(tbl) > 0L) tbl$Model[which.min(tbl$AIC)] else NA_character_

  out <- structure(
    list(
      table  = tbl,
      models = if (return_models) model_list else NULL,
      best   = best
    ),
    class = "nl_compare"
  )

  print(out)
  invisible(out)
}


#' @export
print.nl_compare <- function(x, ...) {
  cat("\n=== MultiSpline Model Comparison ===\n\n")
  print(x$table, row.names = FALSE, ...)
  cat("\n  Best model by AIC:", x$best, "\n\n")
  invisible(x)
}
