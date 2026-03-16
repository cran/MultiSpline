#' Plot predictions from nl_predict()
#'
#' @description
#' Visualizes predicted nonlinear effects returned by
#' \code{\link{nl_predict}}. If a time variable is present (or detected),
#' draws separate colored curves for each time level with optional confidence
#' ribbons.
#'
#' @param pred_df A data frame returned by \code{\link{nl_predict}}. Must
#'   contain at minimum columns \code{fit} and the focal predictor \code{x}.
#'   For confidence bands, provide either \code{lwr}/\code{upr} or
#'   \code{se.fit}.
#' @param x Character string naming the focal predictor column in
#'   \code{pred_df}.
#' @param time Optional character string naming the time variable column in
#'   \code{pred_df}. If \code{NULL}, the function will try to auto-detect a
#'   sensible time column (e.g., \code{"TimePoint"}, \code{"time"},
#'   \code{"wave"}) if present.
#' @param show_ci Logical; if \code{TRUE}, adds a confidence ribbon. Default
#'   \code{TRUE}.
#' @param ci_level Numeric in (0, 1); confidence level used when computing CI
#'   from \code{se.fit}. Default \code{0.95}.
#' @param y_lab Character label for the y-axis. Default
#'   \code{"Predicted outcome"}.
#' @param x_lab Character label for the x-axis. Defaults to the value of
#'   \code{x}.
#' @param title Optional plot title string.
#' @param legend_title Optional legend title string. Defaults to the value of
#'   \code{time}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # --- Toy example (automatically tested by CRAN) ---
#' # Single-level spline: fit, predict, then plot
#' set.seed(1)
#' mydata <- data.frame(
#'   outcome = rnorm(120),
#'   age     = runif(120, 18, 65),
#'   id      = rep(1:30, each = 4)
#' )
#' fit  <- nl_fit(data = mydata, y = "outcome", x = "age", df = 4)
#' pred <- nl_predict(fit)
#' nl_plot(pred, x = "age")
#'
#' \donttest{
#' # Custom axis labels and title
#' nl_plot(
#'   pred,
#'   x     = "age",
#'   y_lab = "Predicted Math Score",
#'   x_lab = "Age (years)",
#'   title = "Nonlinear Effect of Age"
#' )
#'
#' # Without confidence ribbon
#' nl_plot(pred, x = "age", show_ci = FALSE)
#'
#' # With time variable: separate curves per wave
#' set.seed(1)
#' mydata2 <- data.frame(
#'   outcome = rnorm(120),
#'   age     = runif(120, 18, 65),
#'   id      = rep(1:30, each = 4),
#'   wave    = factor(rep(1:4, times = 30))
#' )
#' fit_t  <- nl_fit(
#'   data = mydata2,
#'   y    = "outcome",
#'   x    = "age",
#'   time = "wave",
#'   df   = 4
#' )
#' pred_t <- nl_predict(fit_t)
#' nl_plot(pred_t, x = "age", time = "wave", legend_title = "Wave")
#' }
#'
#' @seealso \code{\link{nl_predict}}, \code{\link{nl_fit}}
#'
#' @importFrom rlang .data
#' @export
nl_plot <- function(
    pred_df,
    x,
    time         = NULL,
    show_ci      = TRUE,
    ci_level     = 0.95,
    y_lab        = "Predicted outcome",
    x_lab        = NULL,
    title        = NULL,
    legend_title = NULL
) {

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!is.data.frame(pred_df)) {
    stop(
      "`pred_df` must be a data.frame returned by nl_predict().",
      call. = FALSE
    )
  }
  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    stop("`x` must be a single non-empty character string.", call. = FALSE)
  }
  if (!x %in% names(pred_df)) {
    stop("`x` ('", x, "') is not a column in `pred_df`.", call. = FALSE)
  }
  if (!"fit" %in% names(pred_df)) {
    stop("`pred_df` must contain a column named 'fit'.", call. = FALSE)
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      ci_level <= 0 || ci_level >= 1) {
    stop(
      "`ci_level` must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Auto-detect time if not supplied
  # ---------------------------------------------------------------------------
  if (is.null(time)) {
    candidates <- c(
      "TimePoint", "time", "Time", "wave", "Wave", "occasion", "Occasion"
    )
    hit <- candidates[candidates %in% names(pred_df)]
    if (length(hit) > 0L) {
      time <- hit[1L]
    }
  } else {
    if (!is.character(time) || length(time) != 1L || !nzchar(time)) {
      stop(
        "`time` must be NULL or a single non-empty character string.",
        call. = FALSE
      )
    }
    if (!time %in% names(pred_df)) {
      warning(
        "`time` ('", time, "') is not a column in `pred_df`. ",
        "Ignoring `time`.",
        call. = FALSE
      )
      time <- NULL
    }
  }

  if (is.null(x_lab))        x_lab        <- x
  if (is.null(legend_title)) legend_title <- time

  # ---------------------------------------------------------------------------
  # Build CI from se.fit if lwr/upr not already present
  # ---------------------------------------------------------------------------
  has_lwr_upr <- all(c("lwr", "upr") %in% names(pred_df))
  has_se      <- "se.fit" %in% names(pred_df) &&
    any(!is.na(pred_df$se.fit))

  if (isTRUE(show_ci) && !has_lwr_upr && has_se) {
    crit        <- stats::qnorm(1 - (1 - ci_level) / 2)
    pred_df$lwr <- pred_df$fit - crit * pred_df$se.fit
    pred_df$upr <- pred_df$fit + crit * pred_df$se.fit
    has_lwr_upr <- TRUE
  }

  # ---------------------------------------------------------------------------
  # Order rows to avoid jagged lines
  # ---------------------------------------------------------------------------
  if (!is.null(time)) {
    pred_df <- pred_df[order(pred_df[[time]], pred_df[[x]]), , drop = FALSE]
  } else {
    pred_df <- pred_df[order(pred_df[[x]]), , drop = FALSE]
  }

  # ---------------------------------------------------------------------------
  # Base aesthetics
  # ---------------------------------------------------------------------------
  if (is.null(time)) {
    p <- ggplot2::ggplot(
      pred_df,
      ggplot2::aes(
        x = .data[[x]],
        y = .data[["fit"]]
      )
    )
  } else {
    p <- ggplot2::ggplot(
      pred_df,
      ggplot2::aes(
        x     = .data[[x]],
        y     = .data[["fit"]],
        color = .data[[time]],
        group = .data[[time]]
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Confidence ribbon
  # ---------------------------------------------------------------------------
  if (isTRUE(show_ci) && has_lwr_upr) {

    if (is.null(time)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["lwr"]],
          ymax = .data[["upr"]]
        ),
        alpha = 0.15,
        color = NA
      )
    } else {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin  = .data[["lwr"]],
          ymax  = .data[["upr"]],
          fill  = .data[[time]],
          group = .data[[time]]
        ),
        alpha = 0.12,
        color = NA
      )
    }

  } else if (isTRUE(show_ci) && !has_lwr_upr) {
    warning(
      "`show_ci = TRUE` but no CI columns (lwr/upr) and no usable se.fit ",
      "found. Plotting without CI.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Lines and theme
  # ---------------------------------------------------------------------------
  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = title,
      color = legend_title,
      fill  = legend_title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = if (!is.null(time)) "right" else "none",
      plot.title      = ggplot2::element_text(face = "bold")
    )

  p
}
