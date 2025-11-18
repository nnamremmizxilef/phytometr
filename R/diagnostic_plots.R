# ===========================================================
# DIAGNOSTIC PLOTS FOR PHYTOMETR
# Comprehensive, holobiont-aware, simulation-aware diagnostics
# ===========================================================

#' Plot environmental space
#'
#' @param env_space An object created by define_env_space().
#'
#' @return A ggplot object showing realised environmental gradients.
#' @export
plot_env_space <- function(env_space) {

  if (!inherits(env_space, "env_space")) {
    stop("Input must be an 'env_space' object.")
  }

  df   <- env_space$realised
  vars <- env_space$variables

  # Melt to long format
  df_long <- reshape2::melt(
    df,
    id.vars      = "group_id",
    measure.vars = vars
  )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x     = group_id,
      y     = value,
      color = variable,
      group = variable
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "Realised Environmental Space",
      x     = "Group / Site",
      y     = "Environmental value",
      color = "Variable"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}

#' Plot individual-level environmental distributions
#'
#' @param ind_env Data frame with individual-level environment.
#'   May contain a \code{group_id} column plus one or more
#'   environmental variables.
#'
#' @return A ggplot object with density plots per variable.
#' @export
plot_env_individual <- function(ind_env) {

  ind_env <- as.data.frame(ind_env)

  vars <- setdiff(colnames(ind_env), "group_id")
  if (length(vars) == 0L) {
    stop("`ind_env` must contain at least one environmental variable.")
  }

  df_long <- tidyr::pivot_longer(
    ind_env,
    cols      = tidyselect::all_of(vars),
    names_to  = "variable",
    values_to = "value"
  )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x    = value,
      fill = if ("group_id" %in% names(ind_env)) factor(group_id) else NULL
    )
  ) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Individual-level Environmental Distributions",
      x     = "Value",
      y     = "Density",
      fill  = if ("group_id" %in% names(ind_env)) "Group ID" else NULL
    )
}

#' Plot host model diagnostics
#'
#' @param host_model Object returned by \code{define_host_model()}.
#'
#' @return A ggplot object with a simple textual summary.
#' @export
plot_host_model <- function(host_model) {

  if (is.null(host_model$type)) {
    stop("`host_model` must have a 'type' element (e.g. 'clone' or 'natural').")
  }

  if (host_model$type == "clone") {
    text <- "Clone system:\nHost variance = 0\n(no host-driven noise)"
  } else {
    text <- sprintf(
      "Natural population\nHost variance fraction = %.2f\nStructure = %s",
      host_model$var_host,
      host_model$pop_structure
    )
  }

  df <- data.frame(x = 1, y = 1, label = text)

  ggplot2::ggplot(df, ggplot2::aes(x, y, label = label)) +
    ggplot2::geom_text(size = 5) +
    ggplot2::theme_void()
}

#' Holobiont component diagnostics
#'
#' PCA + correlation heatmaps for simulated holobiont components.
#'
#' @param components_df Data frame returned by
#'   \code{simulate_holobiont_components()}. Must contain an \code{id}
#'   column plus one or more component columns.
#' @param env Optional environmental data (data frame or matrix) with the
#'   same number of rows as \code{components_df}. Used to compute
#'   env–component correlations.
#' @param labels Optional vector of labels for points in the PCA plot.
#'   Defaults to \code{components_df$id}.
#' @param color_by Optional column name in \code{components_df} or vector
#'   of length \code{nrow(components_df)} used to colour points in PCA.
#'
#' @return A list of ggplot objects. Always returns \code{PCA} and
#'   \code{CorMatrix}; if \code{env} is supplied, also returns
#'   \code{EnvEffects}.
#'
#' @export
plot_holobiont_components <- function(components_df,
                                      env      = NULL,
                                      labels   = NULL,
                                      color_by = NULL) {

  components_df <- as.data.frame(components_df)

  if (!("id" %in% names(components_df))) {
    stop("`components_df` must contain an 'id' column.")
  }

  # Use only true component columns for PCA / correlations
  comp_only <- components_df[ , setdiff(names(components_df), "id"), drop = FALSE]

  ## --- PCA -------------------------------------------------------------

  pca <- stats::prcomp(comp_only, center = TRUE, scale. = TRUE)

  scores <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  # labels
  if (is.null(labels)) {
    scores$label <- components_df$id
  } else {
    if (length(labels) != nrow(components_df)) {
      stop("`labels` must have length nrow(components_df).")
    }
    scores$label <- labels
  }

  # colours
  if (is.null(color_by)) {
    scores$group <- factor(components_df$id)
    colour_lab   <- "ID"
  } else if (is.character(color_by) &&
             length(color_by) == 1L &&
             color_by %in% names(components_df)) {

    scores$group <- components_df[[color_by]]
    colour_lab   <- color_by

  } else if (length(color_by) == nrow(components_df)) {

    scores$group <- color_by
    colour_lab   <- "group"

  } else {
    stop("`color_by` must be NULL, a column name in components_df, or a vector of length nrow(components_df).")
  }

  p_pca <-
    ggplot2::ggplot(
      scores,
      ggplot2::aes(x = PC1, y = PC2, colour = group, label = label)
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    ggplot2::labs(
      title  = "Holobiont Component PCA",
      colour = colour_lab
    ) +
    ggplot2::theme_minimal(base_size = 14)

  ## --- Component correlation heatmap (without id) ---------------------

  cor_mat <- stats::cor(comp_only, use = "pairwise.complete.obs")

  cor_long <- reshape2::melt(
    cor_mat,
    varnames   = c("Var1", "Var2"),
    value.name = "value"
  )

  p_cor <-
    ggplot2::ggplot(
      cor_long,
      ggplot2::aes(x = Var1, y = Var2, fill = value)
    ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      name   = "Correlation",
      low    = "#f7fbff",
      high   = "#08306b",
      limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = "Holobiont Component Correlation Matrix",
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  out <- list(PCA = p_pca, CorMatrix = p_cor)

  ## --- Env effects heatmap (optional) ---------------------------------

  if (!is.null(env)) {
    env_df <- as.data.frame(env)

    if (nrow(env_df) != nrow(components_df)) {
      stop("`env` must have the same number of rows as `components_df`.")
    }

    env_cor <- stats::cor(env_df, comp_only, use = "pairwise.complete.obs")

    env_cor_long <- reshape2::melt(
      env_cor,
      varnames   = c("Environment", "Component"),
      value.name = "Correlation"
    )

    p_env <-
      ggplot2::ggplot(
        env_cor_long,
        ggplot2::aes(x = Environment, y = Component, fill = Correlation)
      ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        name   = "Correlation",
        low    = "#f7fbff",
        high   = "#08306b",
        limits = c(-1, 1)
      ) +
      ggplot2::labs(
        title = "Env–Holobiont Correlation Matrix",
        x     = NULL,
        y     = NULL
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    out$EnvEffects <- p_env
  }

  out
}

#' Plot power landscape over design choices
#'
#' @param power_results Data frame with at least columns
#'   \code{n_groups}, \code{n_per_group}, and \code{power_env}.
#'
#' @return A ggplot object.
#' @export
plot_power_landscape <- function(power_results) {

  ggplot2::ggplot(
    power_results,
    ggplot2::aes(
      x    = n_groups,
      y    = n_per_group,
      fill = power_env
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Power") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Power Landscape",
      x     = "Number of groups",
      y     = "Individuals per group"
    )
}

#' Plot variance decomposition for holobiont components
#'
#' @param power_results A data frame (or list) containing at least
#'   \code{var_host}, the fraction of variance attributed to host.
#'   Only the first row is used.
#' @param components The component specification list used for simulation,
#'   where each component has an element \code{r2_env}.
#'
#' @return A stacked barplot (ggplot object) showing variance fractions
#'   for environment, host, and residual per component.
#'
#' @export
plot_variance_decomposition <- function(power_results, components) {

  var_host <- power_results$var_host[1]

  var_tab <- lapply(names(components), function(comp) {

    r2_env_vec   <- components[[comp]]$r2_env
    r2_env_total <- if (length(r2_env_vec) == 1) r2_env_vec else sum(r2_env_vec)

    var_resid <- 1 - r2_env_total - var_host

    data.frame(
      component = comp,
      var_env   = r2_env_total,
      var_host  = var_host,
      var_resid = var_resid
    )
  })

  var_tab <- do.call(rbind, var_tab)

  df_long <- tidyr::pivot_longer(
    var_tab,
    cols      = c("var_env", "var_host", "var_resid"),
    names_to  = "source",
    values_to = "value"
  )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x    = component,
      y    = value,
      fill = source
    )
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Source") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Variance Decomposition per Component",
      x     = "Component",
      y     = "Proportion of variance"
    )
}
