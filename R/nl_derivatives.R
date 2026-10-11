#' First and second derivatives of the fitted curve
#'
#' @description
#' Computes the slope (first derivative) and curvature (second derivative) of
#' the fitted curve with respect to \code{x}, with delta-method standard
#' errors and pointwise confidence intervals.
#'
#' @details
#' Derivatives are taken of the model itself, not of a table of predictions.
#' For each grid value the spline basis is evaluated at \eqn{x \pm h} with
#' the knots of the fitted model, giving the derivative as a linear
#' combination \eqn{d^\top\beta} (identity link) or, for other links, via
#' the chain rule. Its variance is \eqn{g^\top V g}, where \eqn{g} is the
#' gradient with respect to the fixed effects and \eqn{V} their covariance
#' matrix. This accounts for the strong correlation between neighbouring
#' points on the curve, which a calculation from a prediction table
#' cannot do.
#'
#' Intervals are pointwise, not simultaneous: they describe uncertainty at
#' each \code{x} separately.
#'
#' @param object An \code{nl_fit} object. (Passing a data frame from
#'   \code{nl_predict()}, as in version 0.2.0, still returns derivative
#'   estimates but no intervals, with a warning.)
#' @param x_seq,time_levels,controls_fixed Grid specification, as in
#'   \code{\link{nl_predict}}.
#' @param level Confidence level. Default \code{0.95}.
#' @param type,scale As in \code{\link{nl_predict}}: derivatives of the
#'   conditional or marginal curve, on the response or link scale.
#' @param h Step for the central differences. Default \code{NULL} uses
#'   \eqn{10^{-4}} times the range of \code{x}. Because the basis is
#'   evaluated exactly, the result is insensitive to \code{h}.
#' @param ... Ignored; absorbs arguments of the 0.2.0 interface.
#'
#' @return A data frame with the grid columns, \code{fit}, \code{d1},
#'   \code{d1_se}, \code{d1_lwr}, \code{d1_upr}, \code{d2}, \code{d2_se},
#'   \code{d2_lwr}, and \code{d2_upr}.
#'
#' @seealso \code{\link{nl_turning_points}}, \code{\link{nl_plot}}
#'
#' @examples
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' deriv <- nl_derivatives(fit, x_seq = 0:9)
#' deriv[, c("Days", "d1", "d1_lwr", "d1_upr")]
#'
#' @export
nl_derivatives <- function(
    object,
    x_seq          = NULL,
    time_levels    = NULL,
    controls_fixed = NULL,
    level          = 0.95,
    type           = c("conditional", "marginal"),
    scale          = c("response", "link"),
    h              = NULL,
    ...
) {

  if (is.data.frame(object) && !inherits(object, "nl_fit"))
    return(.derivatives_from_table(object, x_name = x_seq, ...))

  if (!inherits(object, "nl_fit"))
    stop("`object` must be an 'nl_fit' object.", call. = FALSE)
  type  <- match.arg(type)
  scale <- match.arg(scale)
  if (object$method == "gam" && type == "marginal")
    type <- "conditional"

  nd   <- .prediction_grid(object, x_seq, time_levels, controls_fixed)
  x    <- object$x
  crit <- stats::qnorm(1 - (1 - level) / 2)
  if (is.null(h)) h <- 1e-4 * (object$x_info$max - object$x_info$min)
  if (!is.finite(h) || h <= 0) h <- 1e-4

  nd_p <- nd; nd_p[[x]] <- nd_p[[x]] + h
  nd_m <- nd; nd_m[[x]] <- nd_m[[x]] - h

  ev <- function(d) suppressWarnings(.mean_and_grad(object, d, type, scale))
  f0 <- ev(nd); fp <- ev(nd_p); fm <- ev(nd_m)
  V  <- f0$V

  d1   <- (fp$mu - fm$mu) / (2 * h)
  g1   <- (fp$G  - fm$G)  / (2 * h)
  d2   <- (fp$mu - 2 * f0$mu + fm$mu) / h^2
  g2   <- (fp$G  - 2 * f0$G  + fm$G)  / h^2
  se1  <- .qform(g1, V)
  se2  <- .qform(g2, V)

  out        <- nd
  out$fit    <- f0$mu
  out$d1     <- d1
  out$d1_se  <- se1
  out$d1_lwr <- d1 - crit * se1
  out$d1_upr <- d1 + crit * se1
  out$d2     <- d2
  out$d2_se  <- se2
  out$d2_lwr <- d2 - crit * se2
  out$d2_upr <- d2 + crit * se2
  out <- dplyr::as_tibble(out)
  attr(out, "type")  <- type
  attr(out, "scale") <- scale
  attr(out, "x")     <- object$x
  attr(out, "time")  <- object$time
  out
}


# 0.2.0 interface: derivatives from a prediction table. Point estimates only.
.derivatives_from_table <- function(pred_df, x_name = NULL, x = NULL,
                                    time = NULL, ...) {
  x <- x %||% x_name
  if (!is.character(x) || !x %in% names(pred_df))
    stop("Supply the nl_fit object: nl_derivatives(fit).", call. = FALSE)
  warning("nl_derivatives() now takes the fitted model: nl_derivatives(fit). ",
          "From a prediction table only point estimates can be computed; ",
          "intervals are returned as NA because the correlation between ",
          "neighbouring predictions is unknown.", call. = FALSE)
  grp <- if (!is.null(time) && time %in% names(pred_df))
    split(seq_len(nrow(pred_df)), pred_df[[time]]) else list(seq_len(nrow(pred_df)))
  pieces <- lapply(grp, function(ix) {
    d  <- pred_df[ix, , drop = FALSE]
    d  <- d[order(d[[x]]), , drop = FALSE]
    xv <- d[[x]]; yv <- d[["fit"]]
    d1 <- c(NA, diff(yv, lag = 2) / diff(xv, lag = 2), NA)
    d2 <- c(NA, diff(d1, lag = 2) / diff(xv, lag = 2), NA)
    d$d1 <- d1; d$d2 <- d2
    d$d1_se <- d$d1_lwr <- d$d1_upr <- NA_real_
    d$d2_se <- d$d2_lwr <- d$d2_upr <- NA_real_
    d
  })
  dplyr::as_tibble(do.call(rbind, pieces))
}


#' Identify turning points and inflection regions
#'
#' @description
#' From the derivative table returned by \code{\link{nl_derivatives}}, finds
#' turning points (values of \code{x} where the first derivative crosses zero),
#' inflection regions (intervals where the second derivative changes sign),
#' and slope regions (contiguous stretches where the curve is increasing,
#' decreasing, or flat).
#'
#' @param deriv_df A data frame returned by \code{\link{nl_derivatives}}.
#' @param x Character; name of the focal predictor column. Taken from the
#'   \code{nl_derivatives()} output by default.
#' @param time Optional character; name of the time column. Turning points
#'   are identified separately within each time level when supplied.
#' @param tol Numeric tolerance for near-zero detection of d1. Default
#'   \code{1e-6}.
#'
#' @return A named list with three elements.
#'   \code{turning_points} is a data frame with columns \code{x},
#'   \code{time} (if applicable), \code{fit}, and \code{type}
#'   (one of \code{"maximum"}, \code{"minimum"}, or \code{"saddle"}).
#'   \code{inflection_regions} is a data frame with columns
#'   \code{x_start}, \code{x_end}, \code{time}, and \code{direction}
#'   (\code{"concave_to_convex"} or \code{"convex_to_concave"}).
#'   \code{slope_regions} is a data frame with columns \code{x_start},
#'   \code{x_end}, \code{direction}, and \code{time}.
#'
#' @seealso \code{\link{nl_derivatives}}, \code{\link{nl_plot}}
#'
#' @examples
#' fit <- nl_fit(mtcars, y = "mpg", x = "hp", df = 3)
#' tp <- nl_turning_points(nl_derivatives(fit))
#' tp$slope_regions
#'
#' @export
nl_turning_points <- function(
    deriv_df,
    x    = attr(deriv_df, "x"),
    time = attr(deriv_df, "time"),
    tol  = 1e-6
) {

  if (!is.data.frame(deriv_df))
    stop("`deriv_df` must be a data frame from nl_derivatives().", call. = FALSE)
  if (!x %in% names(deriv_df))
    stop("`x` ('", x, "') not found in `deriv_df`.", call. = FALSE)

  # Auto-detect time column
  if (is.null(time)) {
    candidates <- c("TimePoint", "time", "Time", "wave", "Wave",
                    "occasion", "Occasion")
    hit <- candidates[candidates %in% names(deriv_df)]
    if (length(hit)) time <- hit[1L]
  }

  .find_tp <- function(df_sub, time_val = NULL) {

    xv  <- df_sub[[x]]
    d1  <- df_sub[["d1"]]
    d2  <- df_sub[["d2"]]
    fit <- df_sub[["fit"]]
    n   <- length(xv)

    # --- Turning points: d1 crosses zero ---
    tp_list <- list()

    if (n >= 2L) {
      for (i in seq_len(n - 1L)) {
        if (is.finite(d1[i]) && is.finite(d1[i + 1L])) {

          if (d1[i] * d1[i + 1L] < 0) {
            x_cross   <- xv[i] - d1[i] * (xv[i + 1L] - xv[i]) /
              (d1[i + 1L] - d1[i])
            fit_cross <- fit[i] + (fit[i + 1L] - fit[i]) *
              (x_cross - xv[i]) / (xv[i + 1L] - xv[i])
            d2_cross  <- if (is.finite(d2[i])) d2[i] else NA_real_
            tp_type   <- if (!is.na(d2_cross)) {
              if (d2_cross < 0) "maximum" else "minimum"
            } else {
              if (d1[i] > 0) "maximum" else "minimum"
            }
            tp_list[[length(tp_list) + 1L]] <- data.frame(
              x_cross = x_cross,
              fit     = fit_cross,
              type    = tp_type,
              stringsAsFactors = FALSE
            )

          } else if (abs(d1[i]) < tol) {
            tp_type <- if (is.finite(d2[i])) {
              if (d2[i] < 0) "maximum" else if (d2[i] > 0) "minimum" else "saddle"
            } else "saddle"
            tp_list[[length(tp_list) + 1L]] <- data.frame(
              x_cross = xv[i], fit = fit[i], type = tp_type,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    tp_df <- if (length(tp_list)) {
      res <- do.call(rbind, tp_list)
      names(res)[names(res) == "x_cross"] <- x
      if (!is.null(time_val)) res[[time]] <- time_val
      res
    } else {
      base <- data.frame(
        matrix(ncol = 3, nrow = 0,
               dimnames = list(NULL, c(x, "fit", "type")))
      )
      if (!is.null(time_val)) base[[time]] <- character(0)
      base
    }

    # --- Inflection regions: d2 changes sign ---
    infl_list <- list()

    if (n >= 2L && "d2" %in% names(df_sub)) {
      for (i in seq_len(n - 1L)) {
        if (is.finite(d2[i]) && is.finite(d2[i + 1L]) &&
            d2[i] * d2[i + 1L] < 0) {
          dir <- if (d2[i] < 0) "concave_to_convex" else "convex_to_concave"
          infl_list[[length(infl_list) + 1L]] <- data.frame(
            x_start   = xv[i],
            x_end     = xv[i + 1L],
            direction = dir,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    infl_df <- if (length(infl_list)) {
      res <- do.call(rbind, infl_list)
      if (!is.null(time_val)) res[[time]] <- time_val
      res
    } else {
      base <- data.frame(
        x_start   = numeric(0),
        x_end     = numeric(0),
        direction = character(0),
        stringsAsFactors = FALSE
      )
      if (!is.null(time_val)) base[[time]] <- character(0)
      base
    }

    # --- Slope regions: contiguous increasing / decreasing / flat ---
    slope_dir <- ifelse(
      is.na(d1), "flat",
      ifelse(d1 > tol, "increasing",
             ifelse(d1 < -tol, "decreasing", "flat"))
    )
    rle_res  <- rle(slope_dir)
    ends     <- cumsum(rle_res$lengths)
    starts   <- c(1L, ends[-length(ends)] + 1L)

    slope_df <- data.frame(
      x_start   = xv[starts],
      x_end     = xv[ends],
      direction = rle_res$values,
      stringsAsFactors = FALSE
    )
    if (!is.null(time_val)) slope_df[[time]] <- time_val

    list(tp = tp_df, infl = infl_df, slope = slope_df)
  }

  # Run within each time group if applicable
  if (!is.null(time) && time %in% names(deriv_df)) {
    groups    <- split(deriv_df, deriv_df[[time]])
    res_list  <- lapply(names(groups), function(lv) .find_tp(groups[[lv]], lv))
    tp_all    <- do.call(rbind, lapply(res_list, `[[`, "tp"))
    infl_all  <- do.call(rbind, lapply(res_list, `[[`, "infl"))
    slope_all <- do.call(rbind, lapply(res_list, `[[`, "slope"))
  } else {
    res       <- .find_tp(deriv_df)
    tp_all    <- res$tp
    infl_all  <- res$infl
    slope_all <- res$slope
  }

  list(
    turning_points     = tp_all,
    inflection_regions = infl_all,
    slope_regions      = slope_all
  )
}
