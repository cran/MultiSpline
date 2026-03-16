#' Generate predictions from an nl_fit model
#'
#' @description
#' Creates a prediction data frame over a grid of the focal predictor \code{x}
#' (and optionally over \code{time}), holding control variables at typical
#' values (means for numeric; reference levels for factors). For mixed models,
#' predictions default to population-level curves (random effects excluded).
#'
#' @param object An \code{nl_fit} object returned by \code{\link{nl_fit}}.
#' @param x_seq Optional numeric vector of x values to predict over. If
#'   \code{NULL}, uses 100 evenly-spaced points between the stored 1st and
#'   99th percentiles (with fallback to the stored min/max range).
#' @param time_levels Optional vector of time levels to predict for. If
#'   \code{NULL} and a time variable is present, uses stored factor levels.
#' @param controls_fixed Optional named list giving specific values for control
#'   variables. If \code{NULL}, stored defaults are used (mean for numeric;
#'   first level for factors).
#' @param se Logical; if \code{TRUE}, includes standard errors and confidence
#'   intervals where available. Default \code{TRUE}.
#' @param level Confidence level for the interval. Default \code{0.95}.
#' @param re_form For mixed models (\code{lmerMod} / \code{glmerMod}), passed
#'   to \code{predict()}. Default \code{NA} gives population-level predictions
#'   (random effects set to zero).
#' @param ... Reserved for future use.
#'
#' @return
#' A data frame (or tibble) with columns: the focal predictor \code{x},
#' \code{time} (if applicable), any control variables (at fixed values),
#' \code{fit}, \code{se.fit}, \code{lwr}, and \code{upr}.
#'
#' @examples
#' # --- Toy example (automatically tested by CRAN) ---
#' # Single-level natural spline: fit then predict
#' set.seed(1)
#' mydata <- data.frame(
#'   outcome = rnorm(120),
#'   age     = runif(120, 18, 65),
#'   id      = rep(1:30, each = 4)
#' )
#' fit  <- nl_fit(data = mydata, y = "outcome", x = "age", df = 4)
#' pred <- nl_predict(fit)
#' head(pred)
#'
#' \donttest{
#' # Custom x grid
#' pred_custom <- nl_predict(fit, x_seq = seq(20, 60, by = 5))
#'
#' # Without standard errors
#' pred_nose <- nl_predict(fit, se = FALSE)
#'
#' # With time variable (spline x time interaction)
#' set.seed(1)
#' mydata2 <- data.frame(
#'   outcome = rnorm(120),
#'   age     = runif(120, 18, 65),
#'   id      = rep(1:30, each = 4),
#'   wave    = factor(rep(1:4, times = 30))
#' )
#' fit_t <- nl_fit(
#'   data = mydata2,
#'   y    = "outcome",
#'   x    = "age",
#'   time = "wave",
#'   df   = 4
#' )
#' pred_t <- nl_predict(fit_t)
#'
#' # Multilevel model: population-level predictions
#' fit_ml <- nl_fit(
#'   data    = mydata,
#'   y       = "outcome",
#'   x       = "age",
#'   cluster = "id",
#'   df      = 4
#' )
#' pred_ml <- nl_predict(fit_ml)
#' }
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_plot}}
#'
#' @export
nl_predict <- function(
    object,
    x_seq          = NULL,
    time_levels    = NULL,
    controls_fixed = NULL,
    se             = TRUE,
    level          = 0.95,
    re_form        = NA,
    ...
) {

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!inherits(object, "nl_fit")) {
    stop("`object` must be of class 'nl_fit'.", call. = FALSE)
  }

  if (!is.numeric(level) || length(level) != 1L ||
      level <= 0 || level >= 1) {
    stop(
      "`level` must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  mod      <- object$model
  x        <- object$x
  time_var <- object$time
  controls <- object$controls

  # ---------------------------------------------------------------------------
  # 1) Build x grid
  # ---------------------------------------------------------------------------
  xi <- object$x_info
  if (is.null(xi) || !is.list(xi)) {
    stop("Missing `x_info` in `object`.", call. = FALSE)
  }

  if (is.null(x_seq)) {
    lo <- xi$q01
    hi <- xi$q99
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      lo <- xi$min
      hi <- xi$max
    }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      stop(
        "Cannot construct x grid: stored x range is not finite.",
        call. = FALSE
      )
    }
    x_seq <- seq(from = lo, to = hi, length.out = 100)
  } else {
    if (!is.numeric(x_seq) || length(x_seq) < 2L) {
      stop("`x_seq` must be a numeric vector of length >= 2.", call. = FALSE)
    }
    x_seq <- as.numeric(x_seq)
  }

  # ---------------------------------------------------------------------------
  # 2) Time levels
  # ---------------------------------------------------------------------------
  levels_info <- object$levels_info

  if (!is.null(time_var)) {
    if (is.null(time_levels)) {
      if (!is.null(levels_info) && !is.null(levels_info[[time_var]])) {
        time_levels <- levels_info[[time_var]]
      } else {
        mf <- try(stats::model.frame(mod), silent = TRUE)
        if (!inherits(mf, "try-error") &&
            is.data.frame(mf) &&
            time_var %in% names(mf)) {
          time_levels <- unique(mf[[time_var]])
        } else {
          stop(
            "Cannot determine `time_levels`. Please supply `time_levels`.",
            call. = FALSE
          )
        }
      }
    }
  } else {
    time_levels <- NULL
  }

  # ---------------------------------------------------------------------------
  # 3) Build newdata grid
  # ---------------------------------------------------------------------------
  if (is.null(time_levels)) {
    newdata        <- data.frame(x_seq, stringsAsFactors = FALSE)
    names(newdata) <- x
  } else {
    grid <- expand.grid(
      x_seq      = x_seq,
      time_level = time_levels,
      KEEP.OUT.ATTRS   = FALSE,
      stringsAsFactors = FALSE
    )
    names(grid) <- c(x, time_var)
    newdata     <- grid

    if (!is.null(levels_info) && !is.null(levels_info[[time_var]])) {
      newdata[[time_var]] <- factor(
        newdata[[time_var]],
        levels = levels_info[[time_var]]
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Fill in control variables
  # ---------------------------------------------------------------------------
  control_defaults <- object$control_defaults

  if (!is.null(controls) && length(controls) > 0L) {
    for (cn in controls) {
      if (!is.null(controls_fixed) && !is.null(controls_fixed[[cn]])) {
        val <- controls_fixed[[cn]]
      } else if (!is.null(control_defaults) &&
                 !is.null(control_defaults[[cn]])) {
        val <- control_defaults[[cn]]
      } else {
        mf <- try(stats::model.frame(mod), silent = TRUE)
        if (!inherits(mf, "try-error") &&
            is.data.frame(mf) &&
            cn %in% names(mf)) {
          vec <- mf[[cn]]
          if (is.numeric(vec)) {
            val <- mean(vec, na.rm = TRUE)
          } else if (is.factor(vec)) {
            val <- levels(vec)[1L]
          } else {
            val <- vec[which(!is.na(vec))[1L]]
          }
        } else {
          stop(
            "Cannot determine a default value for control `", cn,
            "`. Please supply `controls_fixed`.",
            call. = FALSE
          )
        }
      }
      if (!is.null(levels_info) && !is.null(levels_info[[cn]])) {
        newdata[[cn]] <- factor(val, levels = levels_info[[cn]])
      } else {
        newdata[[cn]] <- val
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 5) Predict
  # ---------------------------------------------------------------------------
  alpha    <- 1 - level
  crit     <- stats::qnorm(1 - alpha / 2)
  fit_vals <- se_vals <- lwr <- upr <- NULL

  # --- lm ---
  if (inherits(mod, "lm")) {
    pred <- stats::predict(mod, newdata = newdata, se.fit = se)
    if (isTRUE(se)) {
      fit_vals <- as.numeric(pred$fit)
      se_vals  <- as.numeric(pred$se.fit)
      lwr      <- fit_vals - crit * se_vals
      upr      <- fit_vals + crit * se_vals
    } else {
      fit_vals <- as.numeric(pred)
      se_vals  <- rep(NA_real_, length(fit_vals))
      lwr      <- rep(NA_real_, length(fit_vals))
      upr      <- rep(NA_real_, length(fit_vals))
    }

    # --- gam ---
  } else if (inherits(mod, "gam")) {
    pred <- mgcv::predict.gam(
      mod, newdata = newdata, se.fit = se, type = "response"
    )
    if (isTRUE(se)) {
      fit_vals <- as.numeric(pred$fit)
      se_vals  <- as.numeric(pred$se.fit)
      lwr      <- fit_vals - crit * se_vals
      upr      <- fit_vals + crit * se_vals
    } else {
      fit_vals <- as.numeric(pred)
      se_vals  <- rep(NA_real_, length(fit_vals))
      lwr      <- rep(NA_real_, length(fit_vals))
      upr      <- rep(NA_real_, length(fit_vals))
    }

    # --- lmerMod ---
  } else if (inherits(mod, "lmerMod")) {
    fit_vals <- as.numeric(
      stats::predict(mod, newdata = newdata, re.form = re_form)
    )
    if (isTRUE(se)) {
      ff <- stats::formula(mod)

      ff_nobars <- if (requireNamespace("reformulas", quietly = TRUE)) {
        reformulas::nobars(ff)
      } else {
        lme4::nobars(ff)
      }

      tt      <- stats::delete.response(stats::terms(ff_nobars))
      X       <- stats::model.matrix(tt, newdata)
      V       <- as.matrix(stats::vcov(mod))
      common  <- intersect(colnames(X), colnames(V))
      if (length(common) == 0L) {
        stop(
          "Cannot compute SE: no overlapping fixed-effect columns.",
          call. = FALSE
        )
      }
      Xc      <- X[, common, drop = FALSE]
      Vc      <- V[common, common, drop = FALSE]
      se_vals <- sqrt(pmax(rowSums((Xc %*% Vc) * Xc), 0))
      lwr     <- fit_vals - crit * se_vals
      upr     <- fit_vals + crit * se_vals
    } else {
      se_vals <- rep(NA_real_, length(fit_vals))
      lwr     <- rep(NA_real_, length(fit_vals))
      upr     <- rep(NA_real_, length(fit_vals))
    }

    # --- glmerMod ---
  } else if (inherits(mod, "glmerMod")) {
    fit_vals <- as.numeric(
      stats::predict(
        mod, newdata = newdata, re.form = re_form, type = "response"
      )
    )
    se_vals <- rep(NA_real_, length(fit_vals))
    lwr     <- rep(NA_real_, length(fit_vals))
    upr     <- rep(NA_real_, length(fit_vals))
    if (isTRUE(se)) {
      warning(
        "SE/CI for glmerMod are not computed by default.",
        call. = FALSE
      )
    }

  } else {
    stop("Unknown underlying model class inside `nl_fit`.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 6) Assemble output
  # ---------------------------------------------------------------------------
  out        <- newdata
  out$fit    <- fit_vals
  out$se.fit <- se_vals
  out$lwr    <- lwr
  out$upr    <- upr

  if (requireNamespace("dplyr", quietly = TRUE)) {
    out <- dplyr::as_tibble(out)
  }

  out
}
