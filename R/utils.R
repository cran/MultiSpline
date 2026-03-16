# ---------------------------------------------------------------------------
# MultiSpline: Package-level imports and internal utilities
# ---------------------------------------------------------------------------

#' @importFrom splines ns
#' @importFrom rlang .data
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# NULL-coalescing operator (internal, not exported)
# Returns `a` if not NULL; otherwise returns `b`.
`%||%` <- function(a, b) if (!is.null(a)) a else b
