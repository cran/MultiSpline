#' Compute first and second derivatives of the fitted spline curve
#'
#' @description
#' Uses numerical differentiation on the prediction grid to compute the
#' first derivative (slope / marginal effect), the second derivative
#' (curvature), and propagated confidence bands for both, using the
#' delta method from the \code{se.fit} column of \code{nl_predict()}.
#' Results are passed to \code{\link{nl_turning_points}} to identify
#' local turning points and inflection regions.
#'
#' @param pred_df A data frame returned by \code{\link{nl_predict}}.
#' @param x Character; name of the focal predictor column.
#' @param time Optional character; name of the time column. If present,
#'   derivatives are computed separately within each time level.
#' @param h Numeric; step size for finite-difference approximation.
#'   Default \code{NULL} uses 0.1 percent of the x range.
#' @param level Confidence level for derivative CIs. Default \code{0.95}.
#'
#' @return A data frame with columns \code{x} (the focal predictor),
#'   \code{time} (if applicable), \code{fit} (predicted outcome),
#'   \code{d1} (first derivative), \code{d1_lwr} and \code{d1_upr}
#'   (lower and upper CI for the first derivative), \code{d2} (second
#'   derivative), and \code{d2_lwr} and \code{d2_upr} (CI for the second
#'   derivative).
#'
#' @seealso \code{\link{nl_turning_points}}, \code{\link{nl_plot}}
#'
#' @export
nl_derivatives <- function(
    pred_df,
    x,
    time  = NULL,
    h     = NULL,
    level = 0.95
) {

  if (!is.data.frame(pred_df))
    stop("`pred_df` must be a data frame from nl_predict().", call. = FALSE)
  if (!x %in% names(pred_df))
    stop("`x` ('", x, "') not found in `pred_df`.", call. = FALSE)
  if (!"fit" %in% names(pred_df))
    stop("`pred_df` must contain a 'fit' column.", call. = FALSE)

  has_se <- "se.fit" %in% names(pred_df) && any(!is.na(pred_df$se.fit))
  crit   <- stats::qnorm(1 - (1 - level) / 2)

  # Auto-detect time column
  if (is.null(time)) {
    candidates <- c("TimePoint", "time", "Time", "wave", "Wave",
                    "occasion", "Occasion")
    hit <- candidates[candidates %in% names(pred_df)]
    if (length(hit)) time <- hit[1L]
  }

  .deriv_group <- function(df_sub) {
    xv  <- df_sub[[x]]
    yv  <- df_sub[["fit"]]
    sev <- if (has_se) df_sub[["se.fit"]] else rep(0, nrow(df_sub))
    n   <- length(xv)

    # Ensure sorted by x
    ord <- order(xv)
    xv  <- xv[ord]; yv <- yv[ord]; sev <- sev[ord]

    if (is.null(h)) h <- diff(range(xv, na.rm = TRUE)) * 0.001
    if (!is.finite(h) || h <= 0) h <- 1e-4

    # --- First derivative: central differences ---
    d1 <- rep(NA_real_, n)

    if (n >= 3L) {
      for (i in seq(2L, n - 1L)) {
        dx <- xv[i + 1L] - xv[i - 1L]
        if (is.finite(dx) && dx > 0)
          d1[i] <- (yv[i + 1L] - yv[i - 1L]) / dx
      }
      dx_f <- xv[2L] - xv[1L]
      dx_b <- xv[n] - xv[n - 1L]
      if (is.finite(dx_f) && dx_f > 0) d1[1L] <- (yv[2L] - yv[1L]) / dx_f
      if (is.finite(dx_b) && dx_b > 0) d1[n]  <- (yv[n] - yv[n - 1L]) / dx_b
    } else if (n == 2L) {
      dx <- xv[2L] - xv[1L]
      if (is.finite(dx) && dx > 0) d1[] <- (yv[2L] - yv[1L]) / dx
    }

    # --- Second derivative: finite difference of d1 ---
    d2 <- rep(NA_real_, n)

    if (n >= 3L) {
      for (i in seq(2L, n - 1L)) {
        dx <- xv[i + 1L] - xv[i - 1L]
        if (is.finite(dx) && dx > 0)
          d2[i] <- (d1[i + 1L] - d1[i - 1L]) / dx
      }
      dx_f <- xv[2L] - xv[1L]
      dx_b <- xv[n] - xv[n - 1L]
      if (is.finite(dx_f) && dx_f > 0)
        d2[1L] <- (d1[2L] - d1[1L]) / dx_f
      if (is.finite(dx_b) && dx_b > 0)
        d2[n]  <- (d1[n] - d1[n - 1L]) / dx_b
    }

    # --- Delta-method SE for d1 ---
    # Var(dy/dx) ~ (se_{i+1}^2 + se_{i-1}^2) / dx^2
    d1_se <- rep(NA_real_, n)

    if (has_se && n >= 3L) {
      for (i in seq(2L, n - 1L)) {
        dx <- xv[i + 1L] - xv[i - 1L]
        if (is.finite(dx) && dx > 0)
          d1_se[i] <- sqrt(sev[i + 1L]^2 + sev[i - 1L]^2) / abs(dx)
      }
      dx_f <- xv[2L] - xv[1L]
      dx_b <- xv[n] - xv[n - 1L]
      if (is.finite(dx_f) && dx_f > 0)
        d1_se[1L] <- sqrt(sev[2L]^2 + sev[1L]^2) / abs(dx_f)
      if (is.finite(dx_b) && dx_b > 0)
        d1_se[n]  <- sqrt(sev[n]^2 + sev[n - 1L]^2) / abs(dx_b)
    }

    # --- Delta-method SE for d2 (propagated from d1_se) ---
    d2_se <- rep(NA_real_, n)

    if (has_se && n >= 3L) {
      for (i in seq(2L, n - 1L)) {
        dx <- xv[i + 1L] - xv[i - 1L]
        if (is.finite(dx) && dx > 0 &&
            is.finite(d1_se[i + 1L]) && is.finite(d1_se[i - 1L]))
          d2_se[i] <- sqrt(d1_se[i + 1L]^2 + d1_se[i - 1L]^2) / abs(dx)
      }
    }

    data.frame(
      x_val  = xv,
      fit    = yv,
      d1     = d1,
      d1_lwr = d1 - crit * d1_se,
      d1_upr = d1 + crit * d1_se,
      d2     = d2,
      d2_lwr = d2 - crit * d2_se,
      d2_upr = d2 + crit * d2_se,
      stringsAsFactors = FALSE
    )
  }

  # Run within each time group if applicable
  if (!is.null(time) && time %in% names(pred_df)) {
    groups <- split(pred_df, pred_df[[time]])
    parts  <- lapply(names(groups), function(lv) {
      res         <- .deriv_group(groups[[lv]])
      res[[time]] <- lv
      res
    })
    out <- do.call(rbind, parts)
  } else {
    out <- .deriv_group(pred_df)
  }

  # Rename x_val back to the original predictor name
  names(out)[names(out) == "x_val"] <- x

  if (requireNamespace("dplyr", quietly = TRUE)) out <- dplyr::as_tibble(out)
  out
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
#' @param x Character; name of the focal predictor column.
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
#' @export
nl_turning_points <- function(
    deriv_df,
    x,
    time = NULL,
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
