#' Generate predictions from an nl_fit model
#'
#' @description
#' Creates a prediction data frame over a grid of the focal predictor \code{x}
#' (and optionally over \code{time}), holding control variables at typical
#' values. For mixed models, predictions default to population-level curves
#' (random effects set to zero).
#'
#' v2 improvements:
#' \itemize{
#'   \item \strong{CI for glmerMod}: approximate confidence intervals are
#'     computed via the parametric bootstrap or the delta method. Set
#'     \code{glmer_ci = "delta"} (default, fast) or \code{"boot"} (more
#'     accurate, slower).
#'   \item \strong{Cluster-specific predictions}: set \code{re_form = NULL}
#'     to include random effects in the predictions.
#' }
#'
#' @param object An \code{nl_fit} object.
#' @param x_seq Optional numeric vector of x values. If \code{NULL}, 200
#'   evenly-spaced points between the 1st and 99th percentiles.
#' @param time_levels Optional vector of time levels.
#' @param controls_fixed Optional named list of fixed control values.
#' @param se Logical; include SEs and CIs. Default \code{TRUE}.
#' @param level Confidence level. Default \code{0.95}.
#' @param re_form For mixed models: \code{NA} (population-level, default) or
#'   \code{NULL} (include random effects).
#' @param glmer_ci Method for glmerMod CIs: \code{"delta"} (default) or
#'   \code{"boot"} (parametric bootstrap).
#' @param n_boot Number of bootstrap replicates when \code{glmer_ci = "boot"}.
#'   Default \code{500}.
#' @param ... Reserved for future use.
#'
#' @return A data frame with columns for the focal predictor, time (if any),
#'   controls at fixed values, \code{fit}, \code{se.fit}, \code{lwr},
#'   and \code{upr}.
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_plot}},
#'   \code{\link{nl_derivatives}}
#'
#' @export
nl_predict <- function(
    object,
    x_seq       = NULL,
    time_levels = NULL,
    controls_fixed = NULL,
    se          = TRUE,
    level       = 0.95,
    re_form     = NA,
    glmer_ci    = c("delta", "boot"),
    n_boot      = 500L,
    ...
) {

  if (!inherits(object, "nl_fit"))
    stop("`object` must be of class 'nl_fit'.", call. = FALSE)
  if (!is.numeric(level) || level <= 0 || level >= 1)
    stop("`level` must be between 0 and 1.", call. = FALSE)

  glmer_ci <- match.arg(glmer_ci)
  mod      <- object$model
  x        <- object$x
  time_var <- object$time
  controls <- object$controls

  # ---------------------------------------------------------------------------
  # 1) Build x grid
  # ---------------------------------------------------------------------------
  xi <- object$x_info
  if (is.null(x_seq)) {
    lo <- xi$q01; hi <- xi$q99
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- xi$min; hi <- xi$max }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo)
      stop("Cannot construct x grid: x range is not finite.", call. = FALSE)
    x_seq <- seq(lo, hi, length.out = 200L)
  } else {
    if (!is.numeric(x_seq) || length(x_seq) < 2L)
      stop("`x_seq` must be a numeric vector of length >= 2.", call. = FALSE)
    x_seq <- as.numeric(x_seq)
  }

  # ---------------------------------------------------------------------------
  # 2) Time levels
  # ---------------------------------------------------------------------------
  levels_info <- object$levels_info
  if (!is.null(time_var) && is.null(time_levels)) {
    if (!is.null(levels_info[[time_var]])) {
      time_levels <- levels_info[[time_var]]
    } else {
      mf <- try(stats::model.frame(mod), silent = TRUE)
      if (!inherits(mf, "try-error") && time_var %in% names(mf))
        time_levels <- unique(mf[[time_var]])
      else
        stop("Cannot determine `time_levels`. Please supply them.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # 3) Build newdata
  # ---------------------------------------------------------------------------
  if (is.null(time_levels)) {
    newdata        <- data.frame(x_seq, stringsAsFactors = FALSE)
    names(newdata) <- x
  } else {
    grid           <- expand.grid(x_seq = x_seq, tl = time_levels,
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    names(grid)    <- c(x, time_var)
    newdata        <- grid
    if (!is.null(levels_info[[time_var]]))
      newdata[[time_var]] <- factor(newdata[[time_var]], levels = levels_info[[time_var]])
  }

  # Fill controls
  if (!is.null(controls)) {
    for (cn in controls) {
      val <- if (!is.null(controls_fixed[[cn]])) controls_fixed[[cn]] else
             object$control_defaults[[cn]]
      if (is.null(val)) stop("Cannot determine default for control `", cn, "`.", call. = FALSE)
      if (!is.null(levels_info[[cn]]))
        newdata[[cn]] <- factor(val, levels = levels_info[[cn]])
      else
        newdata[[cn]] <- val
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Predict by model type
  # ---------------------------------------------------------------------------
  alpha <- 1 - level
  crit  <- stats::qnorm(1 - alpha / 2)

  fit_vals <- se_vals <- lwr <- upr <- NULL

  if (inherits(mod, "lm")) {
    # ------- OLS / lm -------
    pred     <- stats::predict(mod, newdata = newdata, se.fit = se)
    if (isTRUE(se)) {
      fit_vals <- as.numeric(pred$fit)
      se_vals  <- as.numeric(pred$se.fit)
      lwr      <- fit_vals - crit * se_vals
      upr      <- fit_vals + crit * se_vals
    } else {
      fit_vals <- as.numeric(pred)
      se_vals  <- lwr <- upr <- rep(NA_real_, length(fit_vals))
    }

  } else if (inherits(mod, "gam")) {
    # ------- GAM -------
    pred <- mgcv::predict.gam(mod, newdata = newdata, se.fit = se,
                              type = "response")
    if (isTRUE(se)) {
      fit_vals <- as.numeric(pred$fit)
      se_vals  <- as.numeric(pred$se.fit)
      lwr      <- fit_vals - crit * se_vals
      upr      <- fit_vals + crit * se_vals
    } else {
      fit_vals <- as.numeric(pred)
      se_vals  <- lwr <- upr <- rep(NA_real_, length(fit_vals))
    }

  } else if (inherits(mod, "lmerMod")) {
    # ------- LMM (lmer) -------
    fit_vals <- as.numeric(stats::predict(mod, newdata = newdata, re.form = re_form))
    if (isTRUE(se)) {
      ff_nobars <- tryCatch(
        {
          if (requireNamespace("reformulas", quietly = TRUE))
            reformulas::nobars(stats::formula(mod))
          else
            lme4::nobars(stats::formula(mod))
        },
        error = function(e) NULL
      )
      if (!is.null(ff_nobars)) {
        tt <- stats::delete.response(stats::terms(ff_nobars))
        X  <- tryCatch(stats::model.matrix(tt, newdata), error = function(e) NULL)
        V  <- as.matrix(stats::vcov(mod))
        if (!is.null(X)) {
          common  <- intersect(colnames(X), colnames(V))
          if (length(common) > 0L) {
            Xc      <- X[, common, drop = FALSE]
            Vc      <- V[common, common, drop = FALSE]
            se_vals <- sqrt(pmax(rowSums((Xc %*% Vc) * Xc), 0))
            lwr     <- fit_vals - crit * se_vals
            upr     <- fit_vals + crit * se_vals
          } else {
            se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
          }
        } else {
          se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
        }
      } else {
        se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
      }
    } else {
      se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
    }

  } else if (inherits(mod, "glmerMod")) {
    # ------- GLMM (glmer) — v2: CI via delta method or bootstrap -------
    fit_vals <- as.numeric(
      stats::predict(mod, newdata = newdata, re.form = re_form, type = "response")
    )
    if (isTRUE(se)) {
      if (glmer_ci == "delta") {
        # Delta method on link scale → back-transform
        fit_link <- as.numeric(
          stats::predict(mod, newdata = newdata, re.form = re_form, type = "link")
        )
        ff <- stats::formula(mod)
        ff_nobars <- tryCatch(
          if (requireNamespace("reformulas", quietly = TRUE))
            reformulas::nobars(ff) else lme4::nobars(ff),
          error = function(e) NULL
        )
        if (!is.null(ff_nobars)) {
          tt     <- stats::delete.response(stats::terms(ff_nobars))
          X      <- tryCatch(stats::model.matrix(tt, newdata), error = function(e) NULL)
          V      <- as.matrix(stats::vcov(mod))
          if (!is.null(X)) {
            common <- intersect(colnames(X), colnames(V))
            if (length(common) > 0L) {
              Xc       <- X[, common, drop = FALSE]
              Vc       <- V[common, common, drop = FALSE]
              se_link  <- sqrt(pmax(rowSums((Xc %*% Vc) * Xc), 0))
              inv_link <- mod@resp$family$linkinv
              lwr      <- inv_link(fit_link - crit * se_link)
              upr      <- inv_link(fit_link + crit * se_link)
              se_vals  <- (upr - lwr) / (2 * crit)
            } else {
              se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
            }
          } else {
            se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
          }
        } else {
          se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
        }

      } else {
        # Parametric bootstrap CIs
        message("Computing bootstrap CIs for glmerMod (n_boot = ", n_boot, ")...")
        boot_mat <- tryCatch({
          lme4::bootMer(
            mod,
            FUN = function(m) as.numeric(
              stats::predict(m, newdata = newdata, re.form = re_form,
                             type = "response")
            ),
            nsim      = n_boot,
            use.u     = FALSE,
            type      = "parametric",
            .progress = "none"
          )$t
        }, error = function(e) {
          warning("Bootstrap failed: ", conditionMessage(e), call. = FALSE)
          NULL
        })
        if (!is.null(boot_mat)) {
          lwr     <- apply(boot_mat, 2L, stats::quantile,
                           probs = alpha / 2,     na.rm = TRUE)
          upr     <- apply(boot_mat, 2L, stats::quantile,
                           probs = 1 - alpha / 2, na.rm = TRUE)
          se_vals <- apply(boot_mat, 2L, stats::sd, na.rm = TRUE)
        } else {
          se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
        }
      }
    } else {
      se_vals <- lwr <- upr <- rep(NA_real_, length(fit_vals))
    }

  } else {
    stop("Unknown underlying model class inside `nl_fit`.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 5) Assemble output
  # ---------------------------------------------------------------------------
  out        <- newdata
  out$fit    <- fit_vals
  out$se.fit <- se_vals
  out$lwr    <- lwr
  out$upr    <- upr

  if (requireNamespace("dplyr", quietly = TRUE)) out <- dplyr::as_tibble(out)
  out
}
