#' Tidy coefficient table for an nl_fit model
#'
#' @description
#' Produces a tidy coefficient table for a model fitted by
#' \code{\link{nl_fit}}. For linear mixed models (\code{lmerMod}), optional
#' p-values and denominator degrees of freedom are obtained via
#' \pkg{lmerTest} (Satterthwaite method). For GLMMs (\code{glmerMod}),
#' z-tests are reported from \code{summary()}.
#'
#' @param object An \code{nl_fit} object returned by \code{\link{nl_fit}}.
#' @param digits Integer; number of decimal places for rounding. If
#'   \code{NULL}, no rounding is applied. Default \code{3}.
#' @param pvals Logical; if \code{TRUE}, attempts to include p-values.
#'   Default \code{TRUE}.
#' @param df_method For \code{lmerMod} with \code{pvals = TRUE}:
#'   \code{"satterthwaite"} (requires \pkg{lmerTest}) or \code{"none"}.
#'   Default \code{"satterthwaite"}.
#'
#' @return
#' A data frame (or tibble) with columns: \code{Term}, \code{Estimate},
#' \code{Std.Error}, \code{df} (if available), \code{statistic}, and
#' \code{p.value} (if available).
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{summary.nl_fit}}
#'
#' @export
nl_summary <- function(
    object,
    digits    = 3,
    pvals     = TRUE,
    df_method = c("satterthwaite", "none")
) {

  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!inherits(object, "nl_fit")) {
    stop(
      "`object` must be an 'nl_fit' object produced by nl_fit().",
      call. = FALSE
    )
  }

  df_method <- match.arg(df_method)
  mod       <- object$model

  maybe_round <- function(x) {
    if (is.null(digits) || !is.numeric(x)) return(x)
    round(x, digits)
  }

  out <- NULL

  # ---------------------------------------------------------------------------
  # IMPORTANT: "gam" must be checked BEFORE "lm" because mgcv::gam objects
  # inherit from c("gam", "glm", "lm") — checking "lm" first would match GAM
  # models incorrectly and use s$coefficients instead of s$p.table/s$s.table.
  # ---------------------------------------------------------------------------

  # ---------------------------------------------------------------------------
  # 1) gam  <-- MUST come before lm
  # ---------------------------------------------------------------------------
  if (inherits(mod, "gam")) {

    s <- summary(mod)

    # -- Parametric terms (p.table) --
    # Guard against NULL or zero-row table (e.g. smooth-only models)
    p_out <- NULL
    if (!is.null(s$p.table) && nrow(s$p.table) > 0L) {
      ptab  <- as.data.frame(s$p.table)
      t_col <- grep("t value|z value", names(ptab), value = TRUE)[1L]
      p_col <- if ("Pr(>|t|)" %in% names(ptab)) {
        "Pr(>|t|)"
      } else if ("Pr(>|z|)" %in% names(ptab)) {
        "Pr(>|z|)"
      } else {
        NA_character_
      }

      p_out <- data.frame(
        Term      = rownames(ptab),
        Estimate  = ptab[["Estimate"]],
        Std.Error = ptab[["Std. Error"]],
        df        = NA_real_,
        statistic = if (!is.na(t_col)) ptab[[t_col]] else NA_real_,
        p.value   = if (!is.na(p_col)) ptab[[p_col]] else NA_real_,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }

    # -- Smooth terms (s.table) --
    # Guard against NULL or zero-row table
    s_out <- NULL
    if (!is.null(s$s.table) && nrow(s$s.table) > 0L) {
      stab     <- as.data.frame(s$s.table)
      stat_col <- if ("F" %in% names(stab)) {
        "F"
      } else if ("Chi.sq" %in% names(stab)) {
        "Chi.sq"
      } else {
        names(stab)[1L]
      }
      p_col <- if ("p-value" %in% names(stab)) {
        "p-value"
      } else if ("p.value" %in% names(stab)) {
        "p.value"
      } else {
        NA_character_
      }

      s_out <- data.frame(
        Term      = rownames(stab),
        Estimate  = NA_real_,
        Std.Error = NA_real_,
        df        = if ("Ref.df" %in% names(stab)) stab[["Ref.df"]] else
          NA_real_,
        statistic = stab[[stat_col]],
        p.value   = if (!is.na(p_col)) stab[[p_col]] else NA_real_,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }

    # Safely combine: either or both sections may be NULL
    parts <- Filter(Negate(is.null), list(p_out, s_out))
    out   <- if (length(parts) == 0L) {
      data.frame(
        Term      = character(0),
        Estimate  = numeric(0),
        Std.Error = numeric(0),
        df        = numeric(0),
        statistic = numeric(0),
        p.value   = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, parts)
    }

    # ---------------------------------------------------------------------------
    # 2) lm  <-- after gam
    # ---------------------------------------------------------------------------
  } else if (inherits(mod, "lm")) {

    s  <- summary(mod)
    cf <- as.data.frame(s$coefficients)

    out <- data.frame(
      Term      = rownames(cf),
      Estimate  = cf[["Estimate"]],
      Std.Error = cf[["Std. Error"]],
      df        = stats::df.residual(mod),
      statistic = cf[["t value"]],
      p.value   = cf[["Pr(>|t|)"]],
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    # ---------------------------------------------------------------------------
    # 3) lmerMod
    # ---------------------------------------------------------------------------
  } else if (inherits(mod, "lmerMod")) {

    if (isTRUE(pvals) &&
        df_method == "satterthwaite" &&
        requireNamespace("lmerTest", quietly = TRUE)) {

      m_lt <- lmerTest::as_lmerModLmerTest(mod)
      cf   <- as.data.frame(summary(m_lt)$coefficients)

      out <- data.frame(
        Term      = rownames(cf),
        Estimate  = cf[["Estimate"]],
        Std.Error = cf[["Std. Error"]],
        df        = cf[["df"]],
        statistic = cf[["t value"]],
        p.value   = cf[["Pr(>|t|)"]],
        row.names = NULL,
        stringsAsFactors = FALSE
      )

    } else {

      fe <- lme4::fixef(mod)
      se <- sqrt(diag(as.matrix(stats::vcov(mod))))

      out <- data.frame(
        Term      = names(fe),
        Estimate  = unname(fe),
        Std.Error = unname(se),
        df        = NA_real_,
        statistic = unname(fe / se),
        p.value   = NA_real_,
        row.names = NULL,
        stringsAsFactors = FALSE
      )

      if (isTRUE(pvals) &&
          df_method == "satterthwaite" &&
          !requireNamespace("lmerTest", quietly = TRUE)) {
        warning(
          "Install 'lmerTest' to obtain Satterthwaite df and p-values ",
          "for lmerMod.",
          call. = FALSE
        )
      }
    }

    # ---------------------------------------------------------------------------
    # 4) glmerMod
    # ---------------------------------------------------------------------------
  } else if (inherits(mod, "glmerMod")) {

    s  <- summary(mod)
    cf <- as.data.frame(s$coefficients)

    stat_col <- if ("z value" %in% names(cf)) {
      "z value"
    } else if ("t value" %in% names(cf)) {
      "t value"
    } else {
      names(cf)[3L]
    }
    p_col <- if ("Pr(>|z|)" %in% names(cf)) {
      "Pr(>|z|)"
    } else if ("Pr(>|t|)" %in% names(cf)) {
      "Pr(>|t|)"
    } else {
      NA_character_
    }

    out <- data.frame(
      Term      = rownames(cf),
      Estimate  = cf[["Estimate"]],
      Std.Error = cf[["Std. Error"]],
      df        = NA_real_,
      statistic = cf[[stat_col]],
      p.value   = if (!is.na(p_col)) cf[[p_col]] else NA_real_,
      row.names = NULL,
      stringsAsFactors = FALSE
    )

  } else {
    stop(
      "Unknown underlying model class inside `nl_fit`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Drop columns if pvals not requested
  # ---------------------------------------------------------------------------
  if (!isTRUE(pvals)) {
    out$p.value <- NULL
    out$df      <- NULL
  }

  # ---------------------------------------------------------------------------
  # Rounding
  # ---------------------------------------------------------------------------
  out$Estimate  <- maybe_round(out$Estimate)
  out$Std.Error <- maybe_round(out$Std.Error)

  if ("df" %in% names(out)) {
    out$df <- maybe_round(out$df)
  }
  if ("statistic" %in% names(out)) {
    out$statistic <- maybe_round(out$statistic)
  }
  if ("p.value" %in% names(out) && !is.null(digits)) {
    out$p.value <- round(out$p.value, digits + 2L)
  }

  if (requireNamespace("dplyr", quietly = TRUE)) {
    out <- dplyr::as_tibble(out)
  }

  out
}
