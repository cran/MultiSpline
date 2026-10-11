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
#' For Gaussian models \eqn{\sigma^2_\epsilon} is the residual variance.
#' For binary models it is the latent-scale variance (\eqn{\pi^2/3} for
#' logit, 1 for probit, \eqn{\pi^2/6} for complementary log-log), which
#' matches Stata's \code{estat icc}. With nested factors each value is the
#' share of that level alone; the correlation between two units that share
#' a lower-level group is the sum of the shares of that group and every
#' group above it. With random slopes the ICC depends on \code{x}; the
#' function then averages each component over the data and warns.
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
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' nl_icc(fit)
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

  re_var    <- .re_components(mod)
  resid_var <- .dist_var(mod)
  if (any(attr(re_var, "slopes")))
    warning("The model has random slopes, so the ICC varies with x. ",
            "Each component's variance is averaged over the observed data ",
            "(Johnson 2014); report it as an average ICC.", call. = FALSE)
  if (!is.finite(resid_var))
    warning("No latent residual variance for this family/link; ICCs are ",
            "shares of the random-effect variance only.", call. = FALSE)

  total <- sum(re_var)
  if (include_residual && is.finite(resid_var)) total <- total + resid_var
  if (!is.finite(total) || total <= 0)
    stop("Total variance is not finite or positive; cannot compute ICC.",
         call. = FALSE)

  icc <- stats::setNames(as.numeric(re_var) / total, names(re_var))
  if (include_residual && is.finite(resid_var))
    icc <- c(icc, Residual = resid_var / total)
  icc
}
