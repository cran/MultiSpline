#' Plot predictions and derivatives from nl_fit models
#'
#' @description
#' Visualises results from \code{\link{nl_predict}} or
#' \code{\link{nl_derivatives}}.  The \code{type} argument selects the plot:
#' \itemize{
#'   \item \code{"trajectory"}: predicted outcome vs \code{x} with optional
#'     confidence ribbon (v1 behaviour, now default).
#'   \item \code{"slope"}: first derivative (marginal effect) vs \code{x}.
#'   \item \code{"curvature"}: second derivative vs \code{x}.
#'   \item \code{"combo"}: trajectory + slope panels side by side.
#' }
#'
#' @param pred_df A data frame from \code{\link{nl_predict}} (for
#'   \code{type = "trajectory"} or \code{"combo"}).
#' @param deriv_df A data frame from \code{\link{nl_derivatives}} (for
#'   \code{type = "slope"}, \code{"curvature"}, or \code{"combo"}).
#' @param x Character; name of the focal predictor column in \code{pred_df} /
#'   \code{deriv_df}.
#' @param time Optional character; name of the time column. Auto-detected if
#'   \code{NULL}.
#' @param type Plot type: \code{"trajectory"} (default), \code{"slope"},
#'   \code{"curvature"}, or \code{"combo"}.
#' @param show_ci Logical; add confidence ribbons. Default \code{TRUE}.
#' @param ci_level Numeric; confidence level when deriving CI from
#'   \code{se.fit}. Default \code{0.95}.
#' @param show_turning_points Logical; overlay turning-point markers on the
#'   trajectory plot. Default \code{FALSE}.
#' @param turning_points Data frame from \code{\link{nl_turning_points}} (only
#'   used when \code{show_turning_points = TRUE}).
#' @param zero_line Logical; for slope / curvature plots, draw a dashed
#'   horizontal line at zero. Default \code{TRUE}.
#' @param y_lab Y-axis label.
#' @param x_lab X-axis label. Defaults to \code{x}.
#' @param title Optional plot title.
#' @param legend_title Optional legend title.
#'
#' @return A \code{ggplot} object (or a combined plot list for \code{"combo"}).
#'
#' @seealso \code{\link{nl_predict}}, \code{\link{nl_derivatives}},
#'   \code{\link{nl_turning_points}}
#'
#' @importFrom rlang .data
#' @export
nl_plot <- function(
    pred_df          = NULL,
    deriv_df         = NULL,
    x,
    time             = NULL,
    type             = c("trajectory", "slope", "curvature", "combo"),
    show_ci          = TRUE,
    ci_level         = 0.95,
    show_turning_points = FALSE,
    turning_points   = NULL,
    zero_line        = TRUE,
    y_lab            = NULL,
    x_lab            = NULL,
    title            = NULL,
    legend_title     = NULL
) {

  type <- match.arg(type)

  # Auto-detect time
  .detect_time <- function(df) {
    if (is.null(df)) return(NULL)
    candidates <- c("TimePoint","time","Time","wave","Wave","occasion","Occasion")
    hit <- candidates[candidates %in% names(df)]
    if (length(hit)) hit[1L] else NULL
  }
  if (is.null(time)) {
    time <- .detect_time(pred_df) %||% .detect_time(deriv_df)
  }

  if (is.null(x_lab)) x_lab <- x
  if (is.null(legend_title)) legend_title <- time

  # Helper: build CI from se.fit if lwr/upr absent
  .ensure_ci <- function(df) {
    if (!is.data.frame(df)) return(df)
    has_lwr <- all(c("lwr","upr") %in% names(df))
    has_se  <- "se.fit" %in% names(df) && any(!is.na(df$se.fit))
    if (!has_lwr && has_se) {
      crit    <- stats::qnorm(1 - (1 - ci_level) / 2)
      df$lwr  <- df$fit - crit * df$se.fit
      df$upr  <- df$fit + crit * df$se.fit
    }
    df
  }

  # ---- Shared plot base ----
  .base_aes <- function(df, y_col) {
    if (!is.null(time) && time %in% names(df)) {
      df <- df[order(df[[time]], df[[x]]), , drop = FALSE]
    } else {
      df <- df[order(df[[x]]), , drop = FALSE]
    }
    if (!is.null(time) && time %in% names(df)) {
      ggplot2::ggplot(df, ggplot2::aes(
        x = .data[[x]], y = .data[[y_col]],
        color = .data[[time]], group = .data[[time]]
      ))
    } else {
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y_col]]))
    }
  }

  .add_ribbon <- function(p, df, lo_col, hi_col) {
    if (!all(c(lo_col, hi_col) %in% names(df))) return(p)
    if (!is.null(time) && time %in% names(df)) {
      p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]],
                     fill = .data[[time]], group = .data[[time]]),
        alpha = 0.14, color = NA)
    } else {
      p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
        alpha = 0.18, color = NA)
    }
  }

  .theme <- function() {
    list(
      ggplot2::geom_line(linewidth = 0.9),
      ggplot2::theme_minimal(base_size = 12),
      ggplot2::theme(
        legend.position  = if (!is.null(time)) "right" else "none",
        plot.title       = ggplot2::element_text(face = "bold")
      )
    )
  }

  # =========================================================
  # TRAJECTORY plot
  # =========================================================
  .plot_traj <- function() {
    if (is.null(pred_df)) stop("`pred_df` required for type='trajectory'.", call. = FALSE)
    df <- .ensure_ci(as.data.frame(pred_df))
    p  <- .base_aes(df, "fit")
    if (isTRUE(show_ci) && all(c("lwr","upr") %in% names(df)))
      p <- .add_ribbon(p, df, "lwr", "upr")
    p <- p + .theme() +
      ggplot2::labs(x = x_lab,
                    y = if (!is.null(y_lab)) y_lab else "Predicted outcome",
                    title = title, color = legend_title, fill = legend_title)

    if (isTRUE(show_turning_points) && !is.null(turning_points) &&
        nrow(turning_points) > 0L && x %in% names(turning_points)) {
      p <- p +
        ggplot2::geom_vline(
          data    = turning_points,
          mapping = ggplot2::aes(xintercept = .data[[x]],
                                  linetype    = .data[["type"]]),
          color = "grey30", linewidth = 0.6
        ) +
        ggplot2::scale_linetype_manual(
          values = c(maximum = "dashed", minimum = "dotted", saddle = "dotdash")
        )
    }
    p
  }

  # =========================================================
  # SLOPE (first derivative) plot
  # =========================================================
  .plot_slope <- function() {
    if (is.null(deriv_df)) stop("`deriv_df` required for type='slope'.", call. = FALSE)
    df <- as.data.frame(deriv_df)
    p  <- .base_aes(df, "d1")
    if (isTRUE(show_ci) && all(c("d1_lwr","d1_upr") %in% names(df)))
      p <- .add_ribbon(p, df, "d1_lwr", "d1_upr")
    if (isTRUE(zero_line))
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                    color = "grey50", linewidth = 0.5)
    p <- p + .theme() +
      ggplot2::labs(
        x     = x_lab,
        y     = if (!is.null(y_lab)) y_lab else paste0("dy / d(", x, ")"),
        title = title %||% paste0("Marginal Effect (First Derivative)  of ", x),
        color = legend_title, fill = legend_title
      )
    p
  }

  # =========================================================
  # CURVATURE (second derivative) plot
  # =========================================================
  .plot_curv <- function() {
    if (is.null(deriv_df)) stop("`deriv_df` required for type='curvature'.", call. = FALSE)
    df <- as.data.frame(deriv_df)
    if (!"d2" %in% names(df)) stop("`deriv_df` must contain column 'd2'.", call. = FALSE)
    p  <- .base_aes(df, "d2")
    if (isTRUE(show_ci) && all(c("d2_lwr","d2_upr") %in% names(df)))
      p <- .add_ribbon(p, df, "d2_lwr", "d2_upr")
    if (isTRUE(zero_line))
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                    color = "grey50", linewidth = 0.5)
    p <- p + .theme() +
      ggplot2::labs(
        x     = x_lab,
        y     = if (!is.null(y_lab)) y_lab else paste0("d2(y) / d(", x, ")^2"),
        title = title %||% paste0("Curvature (Second Derivative) of ", x),
        color = legend_title, fill = legend_title
      )
    p
  }

  switch(type,
    trajectory = .plot_traj(),
    slope      = .plot_slope(),
    curvature  = .plot_curv(),
    combo      = {
      p1 <- .plot_traj()
      p2 <- .plot_slope()
      list(trajectory = p1, slope = p2)
    }
  )
}
