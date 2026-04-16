#' Intraclass correlation coefficients for a multilevel nl_fit model
#'
#' @description
#' Extracts variance components from a multilevel \code{nl_fit} object and
#' computes intraclass correlation coefficients (ICCs) for each grouping
#' level plus the residual.
#'
#' The ICC for grouping factor \eqn{g} is defined as:
#' \deqn{ICC_g = \frac{\sigma^2_g}{\sum_j \sigma^2_j + \sigma^2_\epsilon}}
#'
#' @param object An \code{nl_fit} object returned by \code{\link{nl_fit}}
#'   that was fitted with one or more \code{cluster} variables (i.e., a
#'   multilevel model fitted via \code{lme4::lmer()} or
#'   \code{lme4::glmer()}).
#' @param include_residual Logical; if \code{TRUE} (default), includes the
#'   residual variance component \code{ICC_resid} in the output so that all
#'   values sum to 1.
#'
#' @return
#' A named numeric vector of ICCs, one per grouping factor (named
#' \code{ICC_<groupname>}) plus \code{Residual} for the residual variance
#' (when \code{include_residual = TRUE}). All values sum to 1.
#'
#' @seealso \code{\link{nl_fit}}
#'
#' @examples
#' \dontrun{
#' fit <- nl_fit(
#'   data    = mydata,
#'   y       = "math_score",
#'   x       = "SES",
#'   cluster = c("id", "schid"),
#'   df      = 4
#' )
#' nl_icc(fit)
#' }
#'
#' @export
nl_icc <- function(object, include_residual = TRUE) {

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!inherits(object, "nl_fit")) {
    stop(
      "`object` must be an 'nl_fit' object produced by nl_fit().",
      call. = FALSE
    )
  }

  mod <- object$model

  if (!inherits(mod, c("lmerMod", "glmerMod"))) {
    stop(
      "`nl_icc()` currently supports lme4 mixed models (lmerMod / glmerMod).",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Extract random-intercept variances
  # ---------------------------------------------------------------------------
  vc   <- lme4::VarCorr(mod)
  vcdf <- as.data.frame(vc)

  # Keep only random-intercept rows (var1 == "(Intercept)", var2 == NA)
  vcdf <- vcdf[
    !is.na(vcdf$grp) &
      vcdf$var1 == "(Intercept)" &
      is.na(vcdf$var2) &
      !is.na(vcdf$vcov),
    ,
    drop = FALSE
  ]

  if (!nrow(vcdf)) {
    stop("No random intercept variances found in the model.", call. = FALSE)
  }

  re_var <- stats::setNames(
    as.numeric(vcdf$vcov),
    as.character(vcdf$grp)
  )
  re_var <- re_var[!is.na(re_var)]

  # ---------------------------------------------------------------------------
  # Residual variance
  # ---------------------------------------------------------------------------
  resid_var <- NA_real_

  if (inherits(mod, "lmerMod")) {
    sc <- attr(vc, "sc")
    if (is.finite(sc)) resid_var <- sc^2
  } else {
    fam <- stats::family(mod)
    if (identical(fam$family, "binomial")) {
      resid_var <- if (identical(fam$link, "logit")) {
        (pi^2) / 3
      } else if (identical(fam$link, "probit")) {
        1
      } else {
        NA_real_
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Compute ICCs
  # ---------------------------------------------------------------------------
  total <- sum(re_var, na.rm = TRUE)
  if (include_residual && is.finite(resid_var)) {
    total <- total + resid_var
  }

  if (!is.finite(total) || total <= 0) {
    stop(
      "Total variance is not finite or positive; cannot compute ICC.",
      call. = FALSE
    )
  }

  icc <- re_var / total

  if (include_residual && is.finite(resid_var)) {
    icc <- c(icc, Residual = resid_var / total)
  }

  icc
}
