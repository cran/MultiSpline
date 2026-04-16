#' Multilevel R-squared decomposition for nl_fit models
#'
#' @description
#' Computes a suite of R-squared statistics for models fitted by
#' \code{\link{nl_fit}}.
#' For single-level models the standard R-squared and adjusted R-squared are
#' returned. For multilevel models (lmerMod / glmerMod) two quantities are
#' reported: the Nakagawa-Schielzeth marginal R-squared (variance explained by
#' fixed effects only) and the conditional R-squared (fixed plus all random
#' effects), together with a level-specific variance partition table analogous
#' to the r2_mlm / Raudenbush-Bryk approach.
#'
#' @details
#' Marginal and conditional R-squared for LMMs follow the Nakagawa and
#' Schielzeth (2013) formulae extended to multiple random effects by Nakagawa,
#' Johnson and Schielzeth (2017). The fixed-effects variance
#' \eqn{\sigma^2_f} is computed as the variance of the linear predictor
#' from fixed effects only (\eqn{\hat{\mu} = X\hat{\beta}}).
#'
#' The level-specific variance partition (r2_mlm-style) decomposes the total
#' modelled variance (\eqn{\sigma^2_f + \sum \sigma^2_j + \sigma^2_\epsilon})
#' to show how much each source contributes, printed as a breakdown table.
#'
#' @param object An \code{nl_fit} object returned by \code{\link{nl_fit}}.
#' @param digits Integer; decimal places for display. Default \code{4}.
#'
#' @return A list of class \code{"nl_r2"} returned invisibly and
#'   pretty-printed automatically. It contains \code{type} (one of
#'   \code{"OLS"}, \code{"GAM"}, \code{"LMM"}, or \code{"GLMM"}),
#'   \code{r2} (a named numeric vector: \code{R2} and \code{R2_adj} for OLS;
#'   \code{R2_dev} for GAM; \code{R2m} and \code{R2c} for LMM/GLMM), and
#'   \code{variance_partition} (a data frame with columns \code{component},
#'   \code{variance}, and \code{proportion} for multilevel models, or
#'   \code{NULL} for single-level models).
#'
#' @references
#' Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
#' obtaining R-squared from generalized linear mixed-effects models.
#' Methods in Ecology and Evolution, 4(2), 133--142.
#'
#' Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient
#' of determination R-squared and intra-class correlation coefficient from
#' generalized linear mixed-effects models revisited and expanded.
#' Journal of the Royal Society Interface, 14(134), 20170213.
#'
#' Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
#' multilevel models: An integrative framework for defining R-squared measures.
#' Psychological Methods, 24(3), 309--338.
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_icc}}
#'
#' @export
nl_r2 <- function(object, digits = 4L) {

  if (!inherits(object, "nl_fit"))
    stop("`object` must be an 'nl_fit' object.", call. = FALSE)

  mod <- object$model

  # ===========================================================================
  # OLS / single-level lm
  # ===========================================================================
  if (inherits(mod, "lm") && !inherits(mod, "gam")) {
    s   <- summary(mod)
    r2  <- c(R2 = s$r.squared, R2_adj = s$adj.r.squared)
    out <- structure(
      list(type = "OLS", r2 = r2, variance_partition = NULL),
      class = "nl_r2"
    )
    print(out, digits = digits)
    return(invisible(out))
  }

  # ===========================================================================
  # GAM -- deviance-based R-squared
  # ===========================================================================
  if (inherits(mod, "gam")) {
    s   <- summary(mod)
    r2  <- c(R2_dev = s$dev.expl)
    out <- structure(
      list(type = "GAM", r2 = r2, variance_partition = NULL),
      class = "nl_r2"
    )
    print(out, digits = digits)
    return(invisible(out))
  }

  # ===========================================================================
  # LMM (lmerMod)
  # ===========================================================================
  if (inherits(mod, "lmerMod")) {

    vc   <- lme4::VarCorr(mod)
    vcdf <- as.data.frame(vc)

    # Random-intercept variances per grouping factor
    re_rows <- vcdf[
      which(!is.na(vcdf$grp) &
            vcdf$var1 == "(Intercept)" &
            is.na(vcdf$var2) &
            !is.na(vcdf$vcov)), , drop = FALSE
    ]

    if (!nrow(re_rows))
      stop(
        "No random intercept variances found. Is this truly a multilevel model?",
        call. = FALSE
      )

    re_var   <- stats::setNames(as.numeric(re_rows$vcov),
                                as.character(re_rows$grp))
    # Defensive: drop any entries where name or value is NA
    re_var   <- re_var[!is.na(names(re_var)) & !is.na(re_var)]
    sigma_e2 <- as.numeric(attr(vc, "sc"))^2   # residual variance, unname

    # Fixed-effects variance: var(X * beta)
    ff_nobars <- tryCatch(
      if (requireNamespace("reformulas", quietly = TRUE))
        reformulas::nobars(stats::formula(mod))
      else
        lme4::nobars(stats::formula(mod)),
      error = function(e) NULL
    )

    sigma_f2 <- if (!is.null(ff_nobars)) {
      mf <- stats::model.frame(mod)
      tt <- stats::delete.response(stats::terms(ff_nobars))
      X  <- tryCatch(stats::model.matrix(tt, mf), error = function(e) NULL)
      if (!is.null(X)) {
        beta   <- lme4::fixef(mod)
        common <- intersect(colnames(X), names(beta))
        Xb     <- X[, common, drop = FALSE] %*% beta[common]
        stats::var(as.numeric(Xb))
      } else NA_real_
    } else NA_real_

    total_re  <- sum(re_var, na.rm = TRUE)
    total_var <- sigma_f2 + total_re + sigma_e2

    R2m <- sigma_f2 / total_var
    R2c <- (sigma_f2 + total_re) / total_var
    r2  <- c(R2m = R2m, R2c = R2c)

    # Variance partition (r2_mlm style)
    parts <- c(re_var, Residual = sigma_e2)
    vp_df <- data.frame(
      component  = c("Fixed effects", names(parts)),
      variance   = c(sigma_f2, as.numeric(parts)),
      proportion = c(sigma_f2, as.numeric(parts)) / total_var,
      stringsAsFactors = FALSE
    )

    out <- structure(
      list(type = "LMM", r2 = r2, variance_partition = vp_df),
      class = "nl_r2"
    )
    print(out, digits = digits)
    return(invisible(out))
  }

  # ===========================================================================
  # GLMM (glmerMod)
  # ===========================================================================
  if (inherits(mod, "glmerMod")) {

    fam  <- stats::family(mod)
    vc   <- lme4::VarCorr(mod)
    vcdf <- as.data.frame(vc)

    re_rows <- vcdf[
      which(!is.na(vcdf$grp) &
            vcdf$var1 == "(Intercept)" &
            is.na(vcdf$var2) &
            !is.na(vcdf$vcov)), , drop = FALSE
    ]

    re_var   <- stats::setNames(as.numeric(re_rows$vcov),
                                as.character(re_rows$grp))
    re_var   <- re_var[!is.na(names(re_var)) & !is.na(re_var)]
    total_re <- sum(re_var, na.rm = TRUE)

    # Fixed-effects variance
    ff_nobars <- tryCatch(
      if (requireNamespace("reformulas", quietly = TRUE))
        reformulas::nobars(stats::formula(mod))
      else
        lme4::nobars(stats::formula(mod)),
      error = function(e) NULL
    )

    sigma_f2 <- if (!is.null(ff_nobars)) {
      mf <- stats::model.frame(mod)
      tt <- stats::delete.response(stats::terms(ff_nobars))
      X  <- tryCatch(stats::model.matrix(tt, mf), error = function(e) NULL)
      if (!is.null(X)) {
        beta   <- lme4::fixef(mod)
        common <- intersect(colnames(X), names(beta))
        Xb     <- X[, common, drop = FALSE] %*% beta[common]
        stats::var(as.numeric(Xb))
      } else NA_real_
    } else NA_real_

    # Distribution-specific variance (link-scale approximation)
    sigma_d2 <- switch(fam$family,
      binomial = if (fam$link == "logit") (pi^2) / 3 else 1,
      poisson  = {
        mu <- stats::fitted(mod)
        log(1 + 1 / mean(mu, na.rm = TRUE))
      },
      NA_real_
    )

    total_var <- sigma_f2 + total_re + sigma_d2
    R2m <- sigma_f2 / total_var
    R2c <- (sigma_f2 + total_re) / total_var
    r2  <- c(R2m = R2m, R2c = R2c)

    parts <- c(re_var, Distribution = sigma_d2)
    vp_df <- data.frame(
      component  = c("Fixed effects", names(parts)),
      variance   = c(sigma_f2, as.numeric(parts)),
      proportion = c(sigma_f2, as.numeric(parts)) / total_var,
      stringsAsFactors = FALSE
    )

    out <- structure(
      list(type = "GLMM", r2 = r2, variance_partition = vp_df),
      class = "nl_r2"
    )
    print(out, digits = digits)
    return(invisible(out))
  }

  stop("Unknown model class inside nl_fit.", call. = FALSE)
}


#' @export
print.nl_r2 <- function(x, digits = 4L, ...) {

  cat("\n=== MultiSpline R-squared Decomposition ===\n")
  cat("Model type:", x$type, "\n\n")

  if (x$type == "OLS") {
    cat(sprintf("  R-squared          = %.4f\n", x$r2["R2"]))
    cat(sprintf("  Adjusted R-squared = %.4f\n", x$r2["R2_adj"]))

  } else if (x$type == "GAM") {
    cat(sprintf("  Deviance explained (R-squared) = %.4f\n", x$r2["R2_dev"]))

  } else {
    # LMM / GLMM
    cat(sprintf("  Marginal  R2m = %.4f   (fixed effects only)\n",
                x$r2["R2m"]))
    cat(sprintf("  Conditional R2c = %.4f  (fixed + all random effects)\n",
                x$r2["R2c"]))

    if (!is.null(x$variance_partition) && nrow(x$variance_partition) > 0L) {
      cat("\nVariance Partition (r2_mlm-style):\n")
      vp <- x$variance_partition
      vp_disp <- data.frame(
        Component  = vp$component,
        Variance   = round(vp$variance,   digits),
        Proportion = round(vp$proportion, digits),
        stringsAsFactors = FALSE
      )
      print(vp_disp, row.names = FALSE)
    }
  }

  cat("\n")
  invisible(x)
}
