#' Fit a nonlinear (spline or GAM) single-level or multilevel model
#'
#' @description
#' Fits a nonlinear regression model for an outcome \code{y} with a focal
#' predictor \code{x}, modeled using natural cubic splines (\code{"ns"}),
#' B-splines (\code{"bs"}), or GAM smooths (\code{"gam"}).
#'
#' Version 2 additions: two-way and nested clustering via the
#' \code{nested} argument; random spline slopes via \code{random_slope};
#' B-spline basis via \code{method = "bs"}; automatic df selection via
#' \code{df = "auto"}.
#'
#' @param data A data frame (often long format for longitudinal data).
#' @param y Outcome variable name (string).
#' @param x Focal nonlinear predictor name (string). Must be numeric.
#' @param time Optional time variable name (string).
#' @param cluster Optional character vector of grouping variable name(s) for
#'   random effects, e.g. \code{NULL}, \code{"id"}, or
#'   \code{c("id", "school")} for two-way clustering.
#' @param nested Logical; only used when \code{length(cluster) == 2}. If
#'   \code{TRUE}, uses a nested specification \code{(1 | g1/g2)}; if
#'   \code{FALSE} (default), uses cross-classified
#'   \code{(1 | g1) + (1 | g2)}.
#' @param controls Optional character vector of additional covariate names
#'   to include as linear fixed effects.
#' @param method Spline basis to use: \code{"ns"} (natural cubic spline,
#'   default), \code{"bs"} (B-spline), or \code{"gam"} (GAM smooth via
#'   \code{mgcv::gam()}). Multilevel fits require \code{"ns"} or \code{"bs"}.
#' @param df Degrees of freedom for the spline basis. Supply a single integer
#'   greater than or equal to 1, or the string \code{"auto"} to trigger
#'   automatic selection by information criterion. Default \code{4}.
#' @param df_range Integer vector of candidate df values evaluated when
#'   \code{df = "auto"}. Default \code{2:8}.
#' @param df_criterion Information criterion used for automatic df selection:
#'   \code{"AIC"} (default) or \code{"BIC"}.
#' @param k Basis dimension for \code{mgcv::s()} when \code{method = "gam"}.
#'   Default \code{5}.
#' @param bs_degree Polynomial degree for \code{splines::bs()} when
#'   \code{method = "bs"}. Default \code{3} (cubic).
#' @param random_slope Logical; if \code{TRUE} and \code{cluster} is supplied,
#'   fits random spline slopes in addition to random intercepts. Currently
#'   implemented only for single-level clustering. Default \code{FALSE}.
#' @param family A family object such as \code{stats::gaussian()} (default)
#'   or \code{stats::binomial()}. For multilevel fits, \code{gaussian()} uses
#'   \code{lme4::lmer()} and all other families use \code{lme4::glmer()}.
#' @param ... Additional arguments passed to the underlying fitting function
#'   (\code{stats::lm()}, \code{mgcv::gam()}, \code{lme4::lmer()}, or
#'   \code{lme4::glmer()}).
#'
#' @return An object of class \code{"nl_fit"} (a named list). It contains
#'   the fitted \code{model} object, the \code{method}, variable names
#'   (\code{y}, \code{x}, \code{time}, \code{cluster}, \code{controls}),
#'   spline settings (\code{df}, \code{df_selected}, \code{df_search},
#'   \code{k}, \code{bs_degree}), flags (\code{nested},
#'   \code{random_slope}), the \code{family}, the model \code{formula},
#'   the \code{call}, and metadata used for prediction (\code{x_info},
#'   \code{levels_info}, \code{control_defaults}).
#'
#' @seealso \code{\link{nl_predict}}, \code{\link{nl_derivatives}},
#'   \code{\link{nl_compare}}, \code{\link{nl_r2}}, \code{\link{nl_knots}}
#'
#' @examples
#' \dontrun{
#' # Single-level natural spline with automatic df selection
#' fit <- nl_fit(data = mydata, y = "score", x = "age", df = "auto")
#'
#' # Two-way cross-classified clustering
#' fit2 <- nl_fit(
#'   data    = mydata,
#'   y       = "score",
#'   x       = "age",
#'   cluster = c("student_id", "school_id"),
#'   nested  = FALSE
#' )
#'
#' # Nested clustering (students within schools)
#' fit3 <- nl_fit(
#'   data    = mydata,
#'   y       = "score",
#'   x       = "age",
#'   cluster = c("student_id", "school_id"),
#'   nested  = TRUE
#' )
#'
#' # Random spline slopes
#' fit4 <- nl_fit(
#'   data         = mydata,
#'   y            = "score",
#'   x            = "age",
#'   cluster      = "id",
#'   random_slope = TRUE
#' )
#' }
#'
#' @export
nl_fit <- function(
    data,
    y,
    x,
    time         = NULL,
    cluster      = NULL,
    nested       = FALSE,
    controls     = NULL,
    method       = c("ns", "bs", "gam"),
    df           = 4,
    df_range     = 2:8,
    df_criterion = c("AIC", "BIC"),
    k            = 5,
    bs_degree    = 3,
    random_slope = FALSE,
    family       = stats::gaussian(),
    ...
) {

  # ---------------------------------------------------------------------------
  # 0) Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)

  .chk1 <- function(z, nm) {
    if (!is.character(z) || length(z) != 1L || !nzchar(z))
      stop("`", nm, "` must be a single non-empty string.", call. = FALSE)
  }
  .chk1(y, "y")
  .chk1(x, "x")

  method       <- match.arg(method)
  df_criterion <- match.arg(df_criterion)

  if (!is.null(time)) .chk1(time, "time")

  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) < 1L)
      stop("`cluster` must be NULL or a character vector.", call. = FALSE)
    if (any(!nzchar(cluster)))
      stop("`cluster` contains empty strings.", call. = FALSE)
    cluster <- unique(cluster)
    if (length(cluster) > 2L)
      stop("MultiSpline v2 supports up to two clustering levels.", call. = FALSE)
  }

  if (!is.null(controls)) {
    if (!is.character(controls))
      stop("`controls` must be NULL or a character vector.", call. = FALSE)
    controls <- unique(controls[nzchar(controls)])
    if (length(controls) == 0L) controls <- NULL
  }

  # df can be a number or "auto"
  auto_df <- FALSE
  if (is.character(df) && length(df) == 1L && tolower(df) == "auto") {
    auto_df <- TRUE
    df      <- 4L   # placeholder; replaced after search
  } else {
    if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1)
      stop("`df` must be a single numeric value >= 1, or \"auto\".",
           call. = FALSE)
    df <- as.integer(df)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 3)
    stop("`k` must be a single numeric value >= 3.", call. = FALSE)

  if (!is.integer(bs_degree) && !is.numeric(bs_degree))
    stop("`bs_degree` must be a positive integer.", call. = FALSE)
  bs_degree <- as.integer(bs_degree)

  if (!is.list(family) || is.null(family$family))
    stop("`family` must be a valid family object.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 1) Check variable existence in data
  # ---------------------------------------------------------------------------
  needed       <- unique(c(y, x, time, controls, cluster))
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L)
    stop("Variables not in `data`: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)

  if (!is.numeric(data[[x]]))
    stop("`x` must be numeric. Got class: ", class(data[[x]])[1L],
         call. = FALSE)

  x_vec <- data[[x]]
  if (all(!is.finite(x_vec)))
    stop("`x` contains no finite values.", call. = FALSE)

  is_multilevel <- !is.null(cluster) && length(cluster) >= 1L

  # ---------------------------------------------------------------------------
  # 2) Automatic df selection
  # ---------------------------------------------------------------------------
  df_search <- NULL

  if (auto_df && method != "gam") {
    df_cands <- as.integer(df_range)
    df_cands <- df_cands[df_cands >= 1L]
    if (length(df_cands) == 0L)
      stop("`df_range` produced no valid candidate values.", call. = FALSE)

    crit_vals <- vapply(df_cands, function(dfi) {
      tryCatch({
        tmp <- .build_and_fit(
          data = data, y = y, x = x, time = time,
          cluster = cluster, nested = nested,
          controls = controls, method = method, df = dfi,
          k = k, bs_degree = bs_degree,
          random_slope = random_slope, family = family, dots = list(...)
        )
        if (df_criterion == "AIC") stats::AIC(tmp) else stats::BIC(tmp)
      }, error = function(e) NA_real_)
    }, numeric(1L))

    df_search <- data.frame(df = df_cands, criterion = crit_vals,
                            stringsAsFactors = FALSE)
    best_idx  <- which.min(crit_vals)
    if (length(best_idx) == 0L)
      stop("Automatic df selection failed for all candidate values.",
           call. = FALSE)
    df <- df_cands[best_idx]
    message("Auto knot selection: best df = ", df,
            " (", df_criterion, " = ", round(crit_vals[best_idx], 2), ")")
  }

  # ---------------------------------------------------------------------------
  # 3) Fit the selected model
  # ---------------------------------------------------------------------------
  model <- .build_and_fit(
    data = data, y = y, x = x, time = time,
    cluster = cluster, nested = nested,
    controls = controls, method = method, df = df,
    k = k, bs_degree = bs_degree,
    random_slope = random_slope, family = family, dots = list(...)
  )

  form <- tryCatch(stats::formula(model), error = function(e) NULL)

  # ---------------------------------------------------------------------------
  # 4) Store compact training metadata
  # ---------------------------------------------------------------------------
  x_q <- suppressWarnings(
    stats::quantile(x_vec, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  )
  x_r <- range(x_vec, na.rm = TRUE)
  if (any(!is.finite(x_q))) x_q <- x_r

  x_info <- list(
    q01 = unname(x_q[1L]), q99 = unname(x_q[2L]),
    min = unname(x_r[1L]), max = unname(x_r[2L]),
    n   = sum(is.finite(x_vec))
  )

  levels_info      <- list()
  control_defaults <- list()

  if (!is.null(time) && is.factor(data[[time]]))
    levels_info[[time]] <- levels(data[[time]])

  if (!is.null(controls)) {
    for (cn in controls) {
      vec <- data[[cn]]
      if (is.factor(vec)) {
        levels_info[[cn]] <- levels(vec)
      } else if (is.character(vec)) {
        # Store unique sorted values so newdata can be properly factored
        # for model.matrix — this prevents column-name mismatches in nl_predict
        levels_info[[cn]] <- sort(unique(vec[!is.na(vec)]))
      }
      control_defaults[[cn]] <- if (is.numeric(vec)) {
        mean(vec, na.rm = TRUE)
      } else if (is.factor(vec)) {
        levels(vec)[1L]
      } else {
        vec[which(!is.na(vec))[1L]]
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
      nested           = nested,
      controls         = controls,
      df               = df,
      df_selected      = df,
      df_search        = df_search,
      k                = k,
      bs_degree        = bs_degree,
      random_slope     = random_slope,
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


# ============================================================
# Internal: construct formula and call the fitting engine
# ============================================================
.build_and_fit <- function(data, y, x, time, cluster, nested,
                           controls, method, df, k, bs_degree,
                           random_slope, family, dots) {

  is_ml <- !is.null(cluster) && length(cluster) >= 1L

  # Fixed-effects spline term
  main_part <- switch(method,

    ns = {
      sc <- paste0("splines::ns(", x, ", df = ", df, ")")
      if (!is.null(time)) paste0(sc, " * ", time) else sc
    },

    bs = {
      sc <- paste0("splines::bs(", x, ", df = ", df,
                   ", degree = ", bs_degree, ")")
      if (!is.null(time)) paste0(sc, " * ", time) else sc
    },

    gam = {
      if (!is.null(time)) {
        time_is_fac <- is.factor(data[[time]]) || is.character(data[[time]])
        if (time_is_fac)
          paste0("s(", x, ", by = ", time, ", k = ", k, ") + ", time)
        else
          paste0("s(", x, ", k = ", k, ") + ", time,
                 " + ti(", x, ", ", time, ", k = ", k, ")")
      } else {
        paste0("s(", x, ", k = ", k, ")")
      }
    }
  )

  fixed_part <- paste(c(main_part, controls), collapse = " + ")

  if (!is_ml) {
    # ---- Single-level ----
    form_env <- if (method == "gam") asNamespace("mgcv") else parent.frame()
    form     <- stats::as.formula(paste0(y, " ~ ", fixed_part),
                                  env = form_env)
    if (method == "gam") {
      do.call(mgcv::gam,
              c(list(formula = form, data = data, family = family), dots))
    } else {
      do.call(stats::lm, c(list(formula = form, data = data), dots))
    }

  } else {
    # ---- Multilevel ----
    if (method == "gam")
      stop("Multilevel support requires method = 'ns' or 'bs'.", call. = FALSE)

    rand_part <- .build_rand_formula(
      cluster      = cluster,
      nested       = nested,
      x            = x,
      df           = df,
      bs_degree    = bs_degree,
      method       = method,
      random_slope = random_slope
    )

    form <- stats::as.formula(
      paste0(y, " ~ ", fixed_part, " + ", rand_part)
    )

    if (identical(family$family, "gaussian")) {
      do.call(lme4::lmer, c(list(formula = form, data = data), dots))
    } else {
      do.call(lme4::glmer,
              c(list(formula = form, data = data, family = family), dots))
    }
  }
}


# ============================================================
# Internal: build the random-effects formula string
# ============================================================
.build_rand_formula <- function(cluster, nested, x, df, bs_degree,
                                method, random_slope) {

  n_clust <- length(cluster)

  if (!random_slope || n_clust != 1L) {
    # Random intercepts only
    if (n_clust == 1L) {
      rand_part <- paste0("(1 | ", cluster[1L], ")")
    } else {
      rand_part <- if (nested) {
        paste0("(1 | ", cluster[1L], "/", cluster[2L], ")")
      } else {
        paste0("(1 | ", cluster[1L], ") + (1 | ", cluster[2L], ")")
      }
    }
  } else {
    # Random spline slopes (single cluster variable only)
    sc <- if (method == "ns") {
      paste0("splines::ns(", x, ", df = ", df, ")")
    } else {
      paste0("splines::bs(", x, ", df = ", df,
             ", degree = ", bs_degree, ")")
    }
    rand_part <- paste0("(", sc, " | ", cluster[1L], ")")
  }

  rand_part
}
