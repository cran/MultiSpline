#' Joint Wald test of the focal predictor
#'
#' @description
#' Tests whether \code{x} is associated with the outcome at all, i.e. whether
#' every fixed-effect coefficient of a term involving \code{x} (all spline
#' basis columns, and their interactions with \code{time}) is zero. The
#' individual spline coefficients reported by \code{\link{nl_summary}} have
#' no separate interpretation; this joint test is the meaningful one. To
#' test nonlinearity (spline versus a straight line) use
#' \code{\link{nl_compare}}.
#'
#' @param object An \code{nl_fit} object fitted with \code{method = "ns"} or
#'   \code{"bs"}.
#'
#' @return A one-row data frame with the Wald chi-square statistic, its
#'   degrees of freedom, and the p-value. Variance components are treated as
#'   known, so for small numbers of clusters the test is approximate.
#'
#' @seealso \code{\link{nl_compare}}, \code{\link{nl_summary}}
#'
#' @examples
#' fit <- nl_fit(lme4::sleepstudy, y = "Reaction", x = "Days",
#'               cluster = "Subject", df = 3)
#' nl_wald(fit)
#'
#' @export
nl_wald <- function(object) {
  if (!inherits(object, "nl_fit"))
    stop("`object` must be an 'nl_fit' object.", call. = FALSE)
  if (object$method == "gam")
    stop("Use summary(fit$model) for smooth-term tests in GAMs.", call. = FALSE)
  mod <- object$model
  X   <- .train_X(mod)
  tl  <- attr(.fixed_terms(mod), "term.labels")
  has_x <- vapply(tl, function(z) object$x %in% all.vars(str2lang(z)),
                  logical(1))
  cols <- which(attr(X, "assign") %in% which(has_x))
  b    <- .coef_fixed(mod)
  keep <- intersect(colnames(X)[cols], names(b))
  V    <- .vcov_fixed(mod)[keep, keep, drop = FALSE]
  stat <- drop(crossprod(b[keep], solve(V, b[keep])))
  k    <- length(keep)
  data.frame(term = object$x, chisq = stat, df = k,
             p.value = stats::pchisq(stat, k, lower.tail = FALSE))
}
