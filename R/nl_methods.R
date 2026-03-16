#' Print method for nl_fit objects
#'
#' @description
#' Compact console display for objects returned by \code{\link{nl_fit}}.
#' Shows key metadata (method, outcome, predictor, time, clustering, family,
#' and controls) and the fitted model formula.
#'
#' @param x An object of class \code{nl_fit}.
#' @param ... Further arguments (currently ignored).
#'
#' @return \code{x} invisibly.
#'
#' @method print nl_fit
#' @export
print.nl_fit <- function(x, ...) {

  if (!inherits(x, "nl_fit")) {
    stop("`x` must be an object of class 'nl_fit'.", call. = FALSE)
  }

  cat("MultiSpline fit  (nl_fit object)\n")
  cat("----------------------------------\n")
  cat("Method:  ", x$method, "\n", sep = "")
  cat("Outcome: ", x$y,      "\n", sep = "")
  cat("Focal x: ", x$x,      "\n", sep = "")

  if (!is.null(x$time)) {
    cat("Time:    ", x$time, "\n", sep = "")
  } else {
    cat("Time:     (none)\n")
  }

  if (is.null(x$cluster) || length(x$cluster) == 0L) {
    cat("Cluster:  (none; single-level)\n")
  } else {
    cat("Cluster: ", paste(x$cluster, collapse = " > "), "\n", sep = "")
  }

  if (!is.null(x$family) && !is.null(x$family$family)) {
    cat("Family:  ", x$family$family, "\n", sep = "")
  }

  if (!is.null(x$controls) && length(x$controls) > 0L) {
    cat("Controls:", paste(x$controls, collapse = ", "), "\n")
  }

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }

  cat("\nModel formula:\n")
  if (!is.null(x$formula)) {
    print(x$formula)
  } else if (!is.null(x$model)) {
    print(stats::formula(x$model))
  } else {
    cat("(formula unavailable)\n")
  }

  cat("\nNext steps:\n")
  cat("  nl_summary(fit)        # tidy coefficient table\n")
  cat("  summary(fit)           # same via S3 summary method\n")
  cat("  nl_predict(fit, ...)   # prediction grid\n")
  cat("  nl_plot(pred, ...)     # visualization\n")

  if (!is.null(x$cluster) && length(x$cluster) > 0L) {
    cat("  nl_icc(fit)            # intraclass correlation\n")
  }

  invisible(x)
}


#' Summary method for nl_fit objects
#'
#' @description
#' Produces a tidy coefficient table via \code{\link{nl_summary}} and wraps
#' it in a \code{summary_nl_fit} object for pretty printing.
#'
#' @param object An \code{nl_fit} object returned by \code{\link{nl_fit}}.
#' @param digits Number of decimal places for rounding. Default \code{3}.
#' @param pvals Logical; if \code{TRUE}, attempts to include p-values.
#'   Default \code{TRUE}.
#' @param df_method For \code{lmerMod}: \code{"satterthwaite"} (requires
#'   \pkg{lmerTest}) or \code{"none"}. Default \code{"satterthwaite"}.
#' @param ... Further arguments passed to \code{\link{nl_summary}}.
#'
#' @return
#' An object of class \code{summary_nl_fit} containing \code{call},
#' \code{formula}, \code{method}, and \code{table}.
#'
#' @method summary nl_fit
#' @export
summary.nl_fit <- function(
    object,
    digits    = 3,
    pvals     = TRUE,
    df_method = c("satterthwaite", "none"),
    ...
) {

  if (!inherits(object, "nl_fit")) {
    stop("`object` must be an object of class 'nl_fit'.", call. = FALSE)
  }

  df_method <- match.arg(df_method)

  tab <- nl_summary(
    object,
    digits    = digits,
    pvals     = pvals,
    df_method = df_method,
    ...
  )

  out <- list(
    call    = object$call %||% NULL,
    formula = tryCatch(
      stats::formula(object$model),
      error = function(e) NULL
    ),
    method  = object$method,
    table   = tab
  )

  class(out) <- "summary_nl_fit"
  out
}


#' Print method for summary_nl_fit objects
#'
#' @description
#' Prints a \code{summary_nl_fit} object returned by
#' \code{\link{summary.nl_fit}}.
#'
#' @param x A \code{summary_nl_fit} object.
#' @param ... Further arguments passed to \code{print}.
#'
#' @return \code{x} invisibly.
#'
#' @method print summary_nl_fit
#' @export
print.summary_nl_fit <- function(x, ...) {

  if (!inherits(x, "summary_nl_fit")) {
    stop("`x` must be an object of class 'summary_nl_fit'.", call. = FALSE)
  }

  cat("Summary of nl_fit object\n")
  cat("  Method:", x$method, "\n")

  if (!is.null(x$formula)) {
    cat("\nFormula:\n")
    print(x$formula)
  }

  cat("\nCoefficient table:\n")
  if (!is.null(x$table)) {
    print(as.data.frame(x$table), row.names = FALSE, ...)
  } else {
    cat("(table unavailable)\n")
  }

  invisible(x)
}




