#' Cluster heterogeneity in nonlinear effects
#'
#' @description
#' Quantifies and visualises how much the nonlinear relationship between
#' \code{x} and \code{y} varies across clusters, using a model with random
#' spline slopes.
#'
#' The function:
#' \enumerate{
#'   \item Refits the model with random spline slopes (if not already fitted
#'     with \code{random_slope = TRUE}).
#'   \item Extracts cluster-specific predicted curves (BLUPs).
#'   \item Returns a heterogeneity summary: variance of slopes across clusters
#'     at each x value, and a plot of the cluster-specific trajectories.
#'   \item Performs a likelihood-ratio test comparing the random-slope model
#'     against a random-intercept-only model to assess whether heterogeneity
#'     is statistically significant.
#' }
#'
#' @param object An \code{nl_fit} object fitted with a single \code{cluster}
#'   variable (and either \code{random_slope = TRUE} or \code{FALSE}; in the
#'   latter case the model is refit internally).
#' @param n_clusters_plot Maximum number of cluster curves to display in the
#'   trajectory plot. Default \code{30L}. Use \code{Inf} for all clusters.
#' @param x_seq Optional numeric vector of x values for prediction.
#' @param level Confidence level for the population-mean ribbon. Default
#'   \code{0.95}.
#' @param plot Logical; print the trajectory plot. Default \code{TRUE}.
#' @param seed Integer seed for reproducibility when sub-sampling clusters.
#'   Default \code{42L}.
#'
#' @return A list (invisibly) with:
#'   \describe{
#'     \item{\code{trajectory_df}}{Long data frame of cluster-specific
#'       predicted values.}
#'     \item{\code{slope_variance}}{Named numeric: SD of first-derivative
#'       estimates across clusters at each x value.}
#'     \item{\code{lrt}}{LRT comparing random-slope vs random-intercept model
#'       (\code{NULL} if \code{random_slope = TRUE} was supplied).}
#'     \item{\code{plot}}{A \code{ggplot} object.}
#'   }
#'
#' @seealso \code{\link{nl_fit}}, \code{\link{nl_derivatives}}
#'
#' @export
nl_het <- function(
    object,
    n_clusters_plot = 30L,
    x_seq           = NULL,
    level           = 0.95,
    plot            = TRUE,
    seed            = 42L
) {

  if (!inherits(object, "nl_fit"))
    stop("`object` must be an 'nl_fit' object.", call. = FALSE)
  if (is.null(object$cluster) || length(object$cluster) != 1L)
    stop("`nl_het()` currently requires exactly one cluster variable.", call. = FALSE)

  cluster_var <- object$cluster[1L]
  mod_ri      <- object$model   # random-intercept model (or the supplied one)

  # Retrieve original data from the model frame
  mf    <- tryCatch(stats::model.frame(mod_ri), error = function(e) NULL)
  if (is.null(mf))
    stop("Cannot extract model frame.", call. = FALSE)

  # Attempt to recover full data (with cluster column)
  orig_data <- tryCatch(
    eval(object$call$data, envir = parent.frame()),
    error = function(e) NULL
  )
  fit_data <- if (!is.null(orig_data) && cluster_var %in% names(orig_data)) {
    orig_data
  } else {
    # Try attaching cluster from model object
    tryCatch({
      re_df  <- as.data.frame(lme4::ranef(mod_ri))
      cl_ids <- rownames(lme4::ranef(mod_ri)[[cluster_var]])
      stop("Cannot recover original data with cluster column. ",
           "Please pass the data via `nl_fit(data = ...)` using a named object.",
           call. = FALSE)
    }, error = function(e) NULL)
  }
  if (is.null(fit_data))
    stop("Cannot recover original data.", call. = FALSE)

  # ---  Refit with random slopes if needed  ---
  lrt <- NULL
  if (!isTRUE(object$random_slope)) {
    message("Refitting with random spline slopes to assess heterogeneity ...")
    mod_rs <- tryCatch(
      .build_and_fit(
        data = fit_data, y = object$y, x = object$x,
        time = object$time, cluster = object$cluster,
        nested = FALSE, controls = object$controls,
        method = object$method, df = object$df,
        k = object$k, bs_degree = object$bs_degree,
        random_slope = TRUE, family = object$family, dots = list()
      ),
      error = function(e) {
        message("Random-slope model failed: ", e$message,
                "\nReturning heterogeneity based on random-intercept BLUPs.")
        NULL
      }
    )
    if (!is.null(mod_rs)) {
      lrt <- tryCatch(
        stats::anova(mod_ri, mod_rs),
        error = function(e) NULL
      )
      mod_rs_use <- mod_rs
    } else {
      mod_rs_use <- mod_ri
    }
  } else {
    mod_rs_use <- mod_ri
  }

  # ---  x grid  ---
  xi <- object$x_info
  if (is.null(x_seq)) {
    lo   <- xi$q01; hi <- xi$q99
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- xi$min; hi <- xi$max }
    x_seq <- seq(lo, hi, length.out = 100L)
  }

  # ---  Cluster IDs  ---
  re         <- lme4::ranef(mod_rs_use)[[cluster_var]]
  all_ids    <- rownames(re)
  n_total    <- length(all_ids)
  set.seed(seed)
  plot_ids   <- if (is.finite(n_clusters_plot) && n_clusters_plot < n_total) {
    sample(all_ids, min(as.integer(n_clusters_plot), n_total))
  } else all_ids

  # ---  Cluster-specific predictions  ---
  x_var    <- object$x
  controls <- object$controls

  # Build newdata for each cluster
  traj_list <- lapply(plot_ids, function(id) {
    nd          <- data.frame(x_seq)
    names(nd)   <- x_var
    nd[[cluster_var]] <- id
    if (!is.null(controls)) {
      for (cn in controls)
        nd[[cn]] <- object$control_defaults[[cn]]
    }
    # Preserve factor levels for cluster variable
    cl_vals <- fit_data[[cluster_var]]
    if (is.factor(cl_vals))
      nd[[cluster_var]] <- factor(id, levels = levels(cl_vals))

    nd$fit <- tryCatch(
      as.numeric(stats::predict(mod_rs_use, newdata = nd, re.form = NULL,
                                allow.new.levels = TRUE)),
      error = function(e) rep(NA_real_, nrow(nd))
    )
    nd$cluster_id <- id
    nd[, c(x_var, "fit", "cluster_id")]
  })
  traj_df <- do.call(rbind, traj_list)

  # Population-level ribbon
  pop_pred <- nl_predict(object, x_seq = x_seq, se = TRUE, level = level,
                         re_form = NA)

  # Slope variance across clusters
  slope_sd <- tapply(traj_df$fit, traj_df[[x_var]], stats::sd, na.rm = TRUE)

  # ---  Plot  ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data    = traj_df,
      mapping = ggplot2::aes(
        x     = .data[[x_var]],
        y     = .data[["fit"]],
        group = .data[["cluster_id"]]
      ),
      alpha    = 0.25,
      linewidth = 0.35,
      color    = "steelblue"
    ) +
    ggplot2::geom_ribbon(
      data    = pop_pred,
      mapping = ggplot2::aes(
        x    = .data[[x_var]],
        ymin = .data[["lwr"]],
        ymax = .data[["upr"]]
      ),
      alpha = 0.30,
      fill  = "#D7191C"
    ) +
    ggplot2::geom_line(
      data    = pop_pred,
      mapping = ggplot2::aes(
        x = .data[[x_var]],
        y = .data[["fit"]]
      ),
      linewidth = 1.1,
      color    = "#D7191C"
    ) +
    ggplot2::labs(
      x     = x_var,
      y     = paste("Predicted", object$y),
      title = paste0("Cluster Heterogeneity in the Effect of ", x_var),
      subtitle = paste0(
        "Blue = cluster-specific curves  |  Red = population mean \u00b1 ",
        round(level * 100), "% CI\n",
        "Displaying ", length(plot_ids), " of ", n_total, " clusters"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, color = "grey40"))

  if (plot) print(p)

  if (!is.null(lrt)) {
    cat("\n--- LRT: Random slopes vs Random intercepts ---\n")
    print(lrt)
  }

  invisible(list(
    trajectory_df  = traj_df,
    slope_variance = slope_sd,
    lrt            = lrt,
    plot           = p
  ))
}
