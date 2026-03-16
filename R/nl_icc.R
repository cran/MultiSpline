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
#'   residual variance component \code{Residual} in the output so that all
#'   values sum to 1 (when a residual variance is available).
#'
#' @return
#' A named numeric vector of ICCs, one per grouping factor (named
#' \code{ICC_<groupname>}) plus \code{Residual} for the residual variance
#' (when \code{include_residual = TRUE}). When a residual variance is available,
#' all values sum to 1.
#'
#' @examples
#' # --- Toy example (automatically tested by CRAN) ---
#' # Small multilevel data: 10 clusters, 5 obs each (50 rows total)
#' set.seed(1)
#' toy <- data.frame(
#'   outcome = rnorm(50),
#'   age     = runif(50, 18, 65),
#'   id      = rep(1:10, each = 5)
#' )
#' fit <- nl_fit(
#'   data    = toy,
#'   y       = "outcome",
#'   x       = "age",
#'   cluster = "id",
#'   df      = 2
#' )
#' nl_icc(fit)
#'
#' \donttest{
#' # Two-level clustering: students nested in schools
#' set.seed(1)
#' mydata <- data.frame(
#'   math_score = rnorm(240),
#'   SES        = rnorm(240),
#'   id         = rep(1:60, each = 4),
#'   schid      = rep(1:12, each = 20)
#' )
#' fit_ml <- nl_fit(
#'   data    = mydata,
#'   y       = "math_score",
#'   x       = "SES",
#'   cluster = c("id", "schid"),
#'   df      = 4
#' )
#' nl_icc(fit_ml)
#'
#' # Without residual variance in output
#' nl_icc(fit_ml, include_residual = FALSE)
#' }
#'
#' @seealso \code{\link{nl_fit}}
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

  # Prefix names as ICC_<group> for clarity
  names(icc) <- ifelse(
    names(icc) == "Residual",
    "Residual",
    paste0("ICC_", names(icc))
  )

  icc
}
