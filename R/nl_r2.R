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
#' With random slopes, each random-effects term contributes the variance of
#' its contribution to the linear predictor averaged over the estimation
#' sample, \eqn{\overline{z_i^\top \Sigma z_i}} (Johnson 2014). For GLMMs
#' the distribution-specific variance is \eqn{\pi^2/3} (logit), 1 (probit),
#' \eqn{\pi^2/6} (complementary log-log), or the lognormal approximation
#' for Poisson models (Nakagawa et al. 2017).
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
#' Johnson, P. C. D. (2014). Extension of Nakagawa & Schielzeth's R2GLMM to
#' random slopes models. Methods in Ecology and Evolution, 5(9), 944--946.
#'
#' Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
#' multilevel models: An integrative framework for defining R-squared measures.
#' Psychological Methods, 24(3), 309--338.
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_icc}}
#'
#' @examples
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' nl_r2(fit)
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
  # Mixed models (lmerMod / glmerMod)
  # ===========================================================================
  if (inherits(mod, "merMod")) {
    re_var   <- .re_components(mod)
    total_re <- sum(re_var)
    sigma_f2 <- .fixed_var(mod)
    sigma_d2 <- .dist_var(mod)
    if (!is.finite(sigma_d2))
      stop("R-squared is not implemented for this family/link.", call. = FALSE)
    if (any(attr(re_var, "slopes")))
      message("Random slopes present: their variance is averaged over the ",
              "observed data (Johnson 2014).")

    total_var <- sigma_f2 + total_re + sigma_d2
    r2 <- c(R2m = sigma_f2 / total_var,
            R2c = (sigma_f2 + total_re) / total_var)

    lmm   <- inherits(mod, "lmerMod")
    parts <- c(re_var, sigma_d2)
    names(parts)[length(parts)] <- if (lmm) "Residual" else "Distribution"
    vp_df <- data.frame(
      component  = c("Fixed effects", names(parts)),
      variance   = c(sigma_f2, as.numeric(parts)),
      proportion = c(sigma_f2, as.numeric(parts)) / total_var,
      stringsAsFactors = FALSE
    )
    out <- structure(
      list(type = if (lmm) "LMM" else "GLMM", r2 = r2,
           variance_partition = vp_df),
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
