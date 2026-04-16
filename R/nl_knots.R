#' Automatic knot / degrees-of-freedom selection for spline models
#'
#' @description
#' Performs a grid search over candidate degrees of freedom (df) for a
#' natural cubic spline (\code{"ns"}) or B-spline (\code{"bs"}) model and
#' selects the best value by an information criterion (AIC or BIC). A
#' diagnostic plot of the criterion against df is returned.
#'
#' This function is called internally by \code{\link{nl_fit}} when
#' \code{df = "auto"}, but it can also be called directly for exploration
#' before fitting the final model.
#'
#' @param data A data frame.
#' @param y Outcome variable name (string).
#' @param x Focal predictor name (string).
#' @param time Optional time variable name.
#' @param cluster Optional character vector of cluster variable names.
#' @param nested Logical; nested clustering. Default \code{FALSE}.
#' @param controls Optional character vector of control variable names.
#' @param method Either \code{"ns"} (default) or \code{"bs"}.
#' @param df_range Integer vector of candidate df values. Default \code{2:10}.
#' @param criterion Either \code{"AIC"} (default) or \code{"BIC"}.
#' @param family A family object. Default \code{stats::gaussian()}.
#' @param plot Logical; if \code{TRUE}, prints a diagnostic plot of criterion
#'   vs df. Default \code{TRUE}.
#' @param ... Additional arguments passed to the underlying fitter.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{best_df}}{The df value with the lowest criterion value.}
#'     \item{\code{search_table}}{Data frame with columns \code{df} and
#'       \code{criterion}.}
#'     \item{\code{criterion}}{The criterion used (\code{"AIC"} or
#'       \code{"BIC"}).}
#'     \item{\code{plot}}{A \code{ggplot} object (if \code{plot = TRUE}).}
#'   }
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_compare}}
#'
#' @export
nl_knots <- function(
    data,
    y,
    x,
    time      = NULL,
    cluster   = NULL,
    nested    = FALSE,
    controls  = NULL,
    method    = c("ns", "bs"),
    df_range  = 2:10,
    criterion = c("AIC", "BIC"),
    family    = stats::gaussian(),
    plot      = TRUE,
    ...
) {

  method    <- match.arg(method)
  criterion <- match.arg(criterion)
  df_cands  <- as.integer(df_range)
  df_cands  <- df_cands[df_cands >= 1L]
  if (length(df_cands) == 0L) stop("`df_range` produced no valid values.", call. = FALSE)

  dots <- list(...)

  crit_vals <- vapply(df_cands, function(dfi) {
    tryCatch({
      mod <- .build_and_fit(
        data = data, y = y, x = x, time = time,
        cluster = cluster, nested = nested,
        controls = controls, method = method, df = dfi,
        k = 5L, bs_degree = 3L,
        random_slope = FALSE, family = family, dots = dots
      )
      if (criterion == "AIC") stats::AIC(mod) else stats::BIC(mod)
    }, error = function(e) NA_real_)
  }, numeric(1L))

  tbl     <- data.frame(df = df_cands, criterion = crit_vals,
                        stringsAsFactors = FALSE)
  best_df <- df_cands[which.min(crit_vals)]

  p <- NULL
  if (plot) {
    p <- ggplot2::ggplot(
      tbl[is.finite(tbl$criterion), ],
      ggplot2::aes(x = .data[["df"]], y = .data[["criterion"]])
    ) +
      ggplot2::geom_line(color = "#2171B5", linewidth = 0.9) +
      ggplot2::geom_point(size = 2.5, color = "#2171B5") +
      ggplot2::geom_vline(xintercept = best_df, linetype = "dashed",
                          color = "#D7191C", linewidth = 0.8) +
      ggplot2::annotate("text", x = best_df, y = max(tbl$criterion, na.rm = TRUE),
                        label = paste0("Best df = ", best_df),
                        hjust = -0.1, vjust = 1.2, color = "#D7191C", size = 3.5) +
      ggplot2::labs(
        x     = "Degrees of Freedom (df)",
        y     = criterion,
        title = paste0("Knot Selection: ", criterion, " by df  (",
                       method, " spline)")
      ) +
      ggplot2::theme_minimal(base_size = 12)
    print(p)
  }

  message("Best df = ", best_df, "  (", criterion, " = ",
          round(crit_vals[which.min(crit_vals)], 2), ")")

  invisible(list(
    best_df      = best_df,
    search_table = tbl,
    criterion    = criterion,
    plot         = p
  ))
}
