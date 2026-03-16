#' Fit a nonlinear (spline or GAM) single-level or multilevel model
#'
#' @description
#' Fits a nonlinear regression model for an outcome \code{y} with a focal
#' predictor \code{x}, modeled either by a natural cubic spline
#' (\code{splines::ns()}) or a GAM smooth (\code{mgcv::s()}).
#'
#' Optionally includes a time variable \code{time}. For spline fits
#' (\code{method = "ns"}), multilevel random-intercept structure can be added
#' via one or more grouping variables in \code{cluster}.
#'
#' Models fitted:
#' \itemize{
#'   \item Single-level spline: \code{stats::lm()}
#'   \item Single-level GAM: \code{mgcv::gam()}
#'   \item Multilevel spline: \code{lme4::lmer()} (Gaussian) or
#'     \code{lme4::glmer()} (non-Gaussian)
#' }
#'
#' @param data A data frame (often long format for longitudinal data).
#' @param y Outcome variable name (string).
#' @param x Focal nonlinear predictor name (string). Must be numeric.
#' @param time Optional time variable name (string). If provided, a
#'   spline-by-time interaction is included for \code{method = "ns"}; for
#'   \code{method = "gam"}, a factor \code{time} uses
#'   \code{s(x, by = time)} (group-specific smooths), while a numeric
#'   \code{time} uses \code{s(x, k = k) + time + ti(x, time, k = k)}
#'   (tensor interaction).
#' @param cluster Optional character vector of grouping variable name(s) for
#'   random intercepts, e.g., \code{NULL}, \code{"id"}, or
#'   \code{c("id", "schid")}.
#' @param controls Optional character vector of additional covariate names to
#'   include linearly.
#' @param method Either \code{"ns"} (natural spline) or \code{"gam"} (GAM
#'   smooth). Multilevel fits are currently supported only for \code{"ns"}.
#' @param df Degrees of freedom for \code{splines::ns()} when
#'   \code{method = "ns"}. Must be a single numeric value >= 1. Default
#'   is \code{4}.
#' @param k Basis dimension for \code{mgcv::s()} and \code{mgcv::ti()} when
#'   \code{method = "gam"}. Must be a single numeric value >= 3. Default
#'   is \code{5}.
#' @param family A model family object, e.g., \code{stats::gaussian()} or
#'   \code{stats::binomial()}. For multilevel fits, \code{gaussian()} uses
#'   \code{lme4::lmer()} and non-Gaussian families use \code{lme4::glmer()}.
#' @param ... Additional arguments passed to the underlying fitting function
#'   (\code{stats::lm()}, \code{mgcv::gam()}, \code{lme4::lmer()}, or
#'   \code{lme4::glmer()}).
#'
#' @return
#' An object of class \code{"nl_fit"} (a named list) containing:
#' \describe{
#'   \item{\code{model}}{The fitted model object.}
#'   \item{\code{method}}{The fitting method used (\code{"ns"} or
#'     \code{"gam"}).}
#'   \item{\code{y}}{Name of the outcome variable.}
#'   \item{\code{x}}{Name of the focal predictor.}
#'   \item{\code{time}}{Name of the time variable, or \code{NULL}.}
#'   \item{\code{cluster}}{Character vector of clustering variables, or
#'     \code{NULL}.}
#'   \item{\code{controls}}{Character vector of control variable names, or
#'     \code{NULL}.}
#'   \item{\code{df}}{Degrees of freedom used for \code{splines::ns()}.}
#'   \item{\code{k}}{Basis dimension used for \code{mgcv::s()} /
#'     \code{mgcv::ti()}.}
#'   \item{\code{family}}{The family object used.}
#'   \item{\code{formula}}{The model formula passed to the fitter.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{x_info}}{A list with quantiles and range of \code{x} for
#'     building prediction grids: \code{q01}, \code{q99}, \code{min},
#'     \code{max}, \code{n}.}
#'   \item{\code{levels_info}}{A named list of factor levels for \code{time}
#'     and any factor control variables.}
#'   \item{\code{control_defaults}}{A named list of typical values for control
#'     variables: the mean for numeric variables and the first level for
#'     factors.}
#' }
#'
#' @examples
#' # --- Toy example (automatically tested by CRAN) ---
#' # Single-level natural spline on small simulated data
#' set.seed(1)
#' mydata <- data.frame(
#'   outcome        = rnorm(120),
#'   age            = runif(120, 18, 65),
#'   id             = rep(1:30, each = 4),
#'   wave           = factor(rep(1:4, times = 30)),
#'   sex            = factor(sample(c("F", "M"), 120, replace = TRUE)),
#'   baseline_score = rnorm(120)
#' )
#'
#' fit_sl <- nl_fit(data = mydata, y = "outcome", x = "age", df = 4)
#'
#' \donttest{
#' # Single-level GAM
#' fit_gam <- nl_fit(
#'   data   = mydata,
#'   y      = "outcome",
#'   x      = "age",
#'   method = "gam",
#'   k      = 5
#' )
#'
#' # Multilevel spline with random intercepts
#' fit_ml <- nl_fit(
#'   data    = mydata,
#'   y       = "outcome",
#'   x       = "age",
#'   cluster = "id",
#'   df      = 4
#' )
#'
#' # Spline with time interaction and controls
#' fit_t <- nl_fit(
#'   data     = mydata,
#'   y        = "outcome",
#'   x        = "age",
#'   time     = "wave",
#'   controls = c("sex", "baseline_score"),
#'   df       = 4
#' )
#' }
#'
#' @seealso
#' \code{\link[splines]{ns}}, \code{\link[mgcv]{gam}},
#' \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}
#'
#' @export
nl_fit <- function(
    data,
    y,
    x,
    time     = NULL,
    cluster  = NULL,
    controls = NULL,
    method   = c("ns", "gam"),
    df       = 4,
    k        = 5,
    family   = stats::gaussian(),
    ...
) {

  # ---------------------------------------------------------------------------
  # 0) Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  .chk_name1 <- function(z, nm) {
    if (!is.character(z) || length(z) != 1L || !nzchar(z)) {
      stop("`", nm, "` must be a single non-empty string.", call. = FALSE)
    }
  }

  .chk_name1(y, "y")
  .chk_name1(x, "x")

  method <- match.arg(method)

  if (!is.null(time)) .chk_name1(time, "time")

  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) < 1L) {
      stop(
        "`cluster` must be NULL or a character vector of grouping variable names.",
        call. = FALSE
      )
    }
    if (any(!nzchar(cluster))) {
      stop("`cluster` contains empty strings.", call. = FALSE)
    }
    cluster <- unique(cluster)
  }

  if (!is.null(controls)) {
    if (!is.character(controls)) {
      stop("`controls` must be NULL or a character vector.", call. = FALSE)
    }
    controls <- unique(controls[nzchar(controls)])
    if (length(controls) == 0L) controls <- NULL
  }

  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("`df` must be a single numeric value >= 1.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 3) {
    stop("`k` must be a single numeric value >= 3.", call. = FALSE)
  }

  if (!is.list(family) || is.null(family$family) ||
      !is.character(family$family)) {
    stop(
      "`family` must be a valid family object ",
      "(e.g., stats::gaussian(), stats::binomial()).",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 1) Check variables exist in data
  # ---------------------------------------------------------------------------
  needed       <- unique(c(y, x, time, controls, cluster))
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "These variables are not in `data`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[x]])) {
    stop(
      "`x` must be numeric. Variable `", x, "` is ",
      class(data[[x]])[1L], ".",
      call. = FALSE
    )
  }
  x_vec <- data[[x]]
  if (all(!is.finite(x_vec))) {
    stop("`x` contains no finite values.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2) Build formula text for fixed effects
  # ---------------------------------------------------------------------------
  main_part <- switch(
    method,

    ns = {
      if (!is.null(time)) {
        paste0("splines::ns(", x, ", df = ", df, ") * ", time)
      } else {
        paste0("splines::ns(", x, ", df = ", df, ")")
      }
    },

    gam = {
      if (!is.null(time)) {
        # mgcv::ti() requires numeric variables only.
        # Factor/character time: use s(x, by = time) for group-specific smooths.
        # Numeric time: use s(x) + time + ti(x, time) for tensor interaction.
        time_is_factor <- is.factor(data[[time]]) ||
          is.character(data[[time]])
        if (time_is_factor) {
          paste0(
            "s(", x, ", by = ", time, ", k = ", k, ") + ",
            time
          )
        } else {
          paste0(
            "s(", x, ", k = ", k, ") + ",
            time, " + ",
            "ti(", x, ", ", time, ", k = ", k, ")"
          )
        }
      } else {
        paste0("s(", x, ", k = ", k, ")")
      }
    }
  )

  fixed_terms <- c(main_part, controls)
  fixed_part  <- paste(fixed_terms, collapse = " + ")

  is_multilevel <- !is.null(cluster) && length(cluster) > 0L

  # ---------------------------------------------------------------------------
  # 3) Fit model
  # ---------------------------------------------------------------------------
  if (!is_multilevel) {

    # For GAM, the formula must be evaluated in the mgcv namespace so that
    # s() and ti() are resolved as smooth terms, not plain list objects.
    form_env <- if (identical(method, "gam")) {
      asNamespace("mgcv")
    } else {
      parent.frame()
    }

    form  <- stats::as.formula(paste0(y, " ~ ", fixed_part), env = form_env)

    model <- if (identical(method, "ns")) {
      stats::lm(form, data = data, ...)
    } else {
      mgcv::gam(form, data = data, family = family, ...)
    }

  } else {

    if (!identical(method, "ns")) {
      stop(
        "Multilevel support is currently implemented only for method = 'ns'.",
        call. = FALSE
      )
    }

    rand_terms <- paste0("(1 | ", cluster, ")")
    rand_part  <- paste(rand_terms, collapse = " + ")
    form <- stats::as.formula(
      paste0(y, " ~ ", fixed_part, " + ", rand_part)
    )

    if (identical(family$family, "gaussian")) {
      model <- lme4::lmer(form, data = data, ...)
    } else {
      model <- lme4::glmer(form, data = data, family = family, ...)
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Store compact training metadata
  # ---------------------------------------------------------------------------
  x_q <- suppressWarnings(
    stats::quantile(x_vec, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  )
  x_r <- range(x_vec, na.rm = TRUE)
  if (any(!is.finite(x_q))) x_q <- x_r

  x_info <- list(
    q01 = unname(x_q[1L]),
    q99 = unname(x_q[2L]),
    min = unname(x_r[1L]),
    max = unname(x_r[2L]),
    n   = sum(is.finite(x_vec))
  )

  levels_info <- list()

  if (!is.null(time) && is.factor(data[[time]])) {
    levels_info[[time]] <- levels(data[[time]])
  }

  if (!is.null(controls) && length(controls) > 0L) {
    for (cn in controls) {
      if (is.factor(data[[cn]])) {
        levels_info[[cn]] <- levels(data[[cn]])
      }
    }
  }

  control_defaults <- list()
  if (!is.null(controls) && length(controls) > 0L) {
    for (cn in controls) {
      vec <- data[[cn]]
      if (is.numeric(vec)) {
        control_defaults[[cn]] <- mean(vec, na.rm = TRUE)
      } else if (is.factor(vec)) {
        control_defaults[[cn]] <- levels(vec)[1L]
      } else {
        idx <- which(!is.na(vec))[1L]
        control_defaults[[cn]] <- if (length(idx)) vec[idx] else NA
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 5) Return nl_fit object
  # ---------------------------------------------------------------------------
  structure(
    list(
      model            = model,
      method           = method,
      y                = y,
      x                = x,
      time             = time,
      cluster          = if (is_multilevel) cluster else NULL,
      controls         = controls,
      df               = df,
      k                = k,
      family           = family,
      formula          = form,
      call             = match.call(),
      x_info           = x_info,
      levels_info      = levels_info,
      control_defaults = control_defaults
    ),
    class = "nl_fit"
  )
}
