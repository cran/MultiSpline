#' Generate predictions from an nl_fit model
#'
#' @description
#' Evaluates the fitted curve over a grid of the focal predictor \code{x}
#' (and over \code{time} levels, if the model has a time interaction),
#' holding control variables at fixed values.
#'
#' @section What is predicted:
#' \itemize{
#'   \item The grid has 200 evenly spaced values of \code{x} between its 1st
#'     and 99th sample percentiles unless \code{x_seq} is supplied. Each row
#'     is an exact evaluation of the fitted model; nothing is interpolated.
#'     Plots join these points with straight line segments.
#'   \item Controls are held at their sample mean (numeric) or first level
#'     (factor / character) unless set in \code{controls_fixed}.
#'   \item \code{type = "conditional"} (default) sets all random effects to
#'     zero: the curve for a cluster whose random effects are at their mean.
#'     For Gaussian models this equals the population-averaged curve.
#'   \item \code{type = "marginal"} averages the response-scale mean over
#'     the estimated random-effect distribution (Gauss-Hermite quadrature).
#'     For non-identity links (e.g. logistic models) this population-averaged
#'     curve is flatter than the conditional curve. Available for random
#'     intercepts and for the random spline slopes created by
#'     \code{nl_fit(random_slope = TRUE)}.
#'   \item \code{scale = "response"} (default) returns means on the outcome
#'     scale (probabilities for binary outcomes); \code{scale = "link"}
#'     returns the linear predictor.
#'   \item Standard errors use the delta method with the full covariance
#'     matrix of the fixed effects, treating variance components as known.
#'     Intervals are formed on the link scale and transformed, so they stay
#'     inside the range of the outcome.
#' }
#'
#' @param object An \code{nl_fit} object.
#' @param x_seq Optional numeric vector of x values.
#' @param time_levels Optional vector of time levels.
#' @param controls_fixed Optional named list of fixed control values.
#' @param se Logical; include standard errors and intervals. Default
#'   \code{TRUE}.
#' @param level Confidence level. Default \code{0.95}.
#' @param type \code{"conditional"} (default) or \code{"marginal"}; see
#'   \emph{What is predicted}.
#' @param scale \code{"response"} (default) or \code{"link"}.
#' @param re_form Retained for compatibility. Only \code{NA} (random effects
#'   set to zero) is supported; use \code{\link{nl_het}} for
#'   cluster-specific curves.
#' @param glmer_ci For non-Gaussian mixed models: \code{"delta"} (default)
#'   or \code{"boot"} (parametric bootstrap via \code{lme4::bootMer()}).
#' @param n_boot Number of bootstrap replicates when \code{glmer_ci = "boot"}.
#' @param ... Reserved for future use.
#'
#' @return A data frame (tibble) with the grid columns, \code{fit},
#'   \code{se.fit}, \code{lwr}, and \code{upr}. Attributes \code{type} and
#'   \code{scale} record what was predicted.
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_plot}},
#'   \code{\link{nl_derivatives}}
#'
#' @examples
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' pred <- nl_predict(fit)
#' head(pred)
#' attr(pred, "type")
#'
#' @export
nl_predict <- function(
    object,
    x_seq          = NULL,
    time_levels    = NULL,
    controls_fixed = NULL,
    se             = TRUE,
    level          = 0.95,
    type           = c("conditional", "marginal"),
    scale          = c("response", "link"),
    re_form        = NA,
    glmer_ci       = c("delta", "boot"),
    n_boot         = 500L,
    ...
) {

  if (!inherits(object, "nl_fit"))
    stop("`object` must be of class 'nl_fit'.", call. = FALSE)
  if (!is.numeric(level) || level <= 0 || level >= 1)
    stop("`level` must be between 0 and 1.", call. = FALSE)

  type     <- match.arg(type)
  scale    <- match.arg(scale)
  glmer_ci <- match.arg(glmer_ci)
  mod      <- object$model

  if (is.null(re_form) || !(length(re_form) == 1L && is.na(re_form)))
    stop("Only `re_form = NA` is supported. For cluster-specific curves use ",
         "nl_het().", call. = FALSE)

  newdata <- .prediction_grid(object, x_seq, time_levels, controls_fixed)

  crit <- stats::qnorm(1 - (1 - level) / 2)
  fam  <- .model_family(object)

  lp      <- .lp(object, newdata)
  eta     <- lp$eta
  se_link <- .qform(lp$X, lp$V)

  if (scale == "link" || identical(fam$link, "identity")) {
    fit_vals <- eta
    se_vals  <- se_link
    lwr      <- fit_vals - crit * se_vals
    upr      <- fit_vals + crit * se_vals
  } else if (type == "conditional") {
    fit_vals <- fam$linkinv(eta)
    lwr      <- fam$linkinv(eta - crit * se_link)
    upr      <- fam$linkinv(eta + crit * se_link)
    se_vals  <- abs(fam$mu.eta(eta)) * se_link
  } else {
    mg       <- .mean_and_grad(object, newdata, "marginal", "response")
    fit_vals <- mg$mu
    se_vals  <- .qform(mg$G, mg$V)
    eta_m    <- fam$linkfun(.clamp_mu(fit_vals, fam))
    se_eta_m <- se_vals / abs(fam$mu.eta(eta_m))
    lwr      <- fam$linkinv(eta_m - crit * se_eta_m)
    upr      <- fam$linkinv(eta_m + crit * se_eta_m)
  }

  if (isTRUE(se) && inherits(mod, "glmerMod") && glmer_ci == "boot") {
    message("Computing bootstrap intervals (n_boot = ", n_boot, ") ...")
    obj_b <- object
    boot_mat <- tryCatch(
      lme4::bootMer(
        mod,
        FUN = function(m) {
          obj_b$model <- m
          .mean_and_grad(obj_b, newdata, type, scale)$mu
        },
        nsim = n_boot, use.u = FALSE, type = "parametric"
      )$t,
      error = function(e) {
        warning("Bootstrap failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(boot_mat)) {
      a       <- 1 - level
      lwr     <- apply(boot_mat, 2L, stats::quantile, a / 2, na.rm = TRUE)
      upr     <- apply(boot_mat, 2L, stats::quantile, 1 - a / 2, na.rm = TRUE)
      se_vals <- apply(boot_mat, 2L, stats::sd, na.rm = TRUE)
    }
  }

  if (!isTRUE(se)) se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))

  out        <- newdata
  out$fit    <- as.numeric(fit_vals)
  out$se.fit <- as.numeric(se_vals)
  out$lwr    <- as.numeric(lwr)
  out$upr    <- as.numeric(upr)
  out <- dplyr::as_tibble(out)
  attr(out, "type")  <- type
  attr(out, "scale") <- scale
  attr(out, "x")     <- object$x
  attr(out, "time")  <- object$time
  out
}

.clamp_mu <- function(mu, fam) {
  if (fam$family %in% c("binomial", "quasibinomial"))
    pmin(pmax(mu, 1e-12), 1 - 1e-12)
  else pmax(mu, 1e-12)
}


# ---------------------------------------------------------------------------
# Internal: build the prediction grid (x, time, controls)
# ---------------------------------------------------------------------------
.prediction_grid <- function(object, x_seq = NULL, time_levels = NULL,
                             controls_fixed = NULL) {
  x        <- object$x
  time_var <- object$time
  xi       <- object$x_info

  if (is.null(x_seq)) {
    lo <- xi$q01; hi <- xi$q99
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- xi$min; hi <- xi$max }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo)
      stop("Cannot construct x grid: x range is not finite.", call. = FALSE)
    x_seq <- seq(lo, hi, length.out = 200L)
  } else {
    if (!is.numeric(x_seq) || length(x_seq) < 1L)
      stop("`x_seq` must be a numeric vector.", call. = FALSE)
    x_seq <- as.numeric(x_seq)
  }

  levels_info <- object$levels_info
  if (!is.null(time_var) && is.null(time_levels)) {
    time_levels <- levels_info[[time_var]] %||%
      sort(unique(object$data[[time_var]]))
    if (is.null(time_levels))
      stop("Cannot determine `time_levels`. Please supply them.", call. = FALSE)
  }

  if (is.null(time_var)) {
    newdata <- stats::setNames(data.frame(x_seq), x)
  } else {
    newdata <- expand.grid(x_seq, time_levels, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
    names(newdata) <- c(x, time_var)
    if (!is.null(levels_info[[time_var]]))
      newdata[[time_var]] <- factor(newdata[[time_var]],
                                    levels = levels_info[[time_var]])
  }

  ctrl_vars <- object$control_vars %||% object$controls
  for (cn in ctrl_vars) {
    val <- controls_fixed[[cn]] %||% object$control_defaults[[cn]]
    if (is.null(val))
      stop("Cannot determine a value for control `", cn, "`.", call. = FALSE)
    newdata[[cn]] <- if (!is.null(levels_info[[cn]]))
      factor(val, levels = levels_info[[cn]]) else val
  }
  newdata
}
