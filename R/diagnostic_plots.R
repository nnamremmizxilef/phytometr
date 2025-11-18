# ===========================================================
# DIAGNOSTIC PLOTS FOR PHYTOMETR
# Comprehensive, holobiont-aware, simulation-aware diagnostics
# ===========================================================

#' Plot environmental space
#'
#' @param env_space An object created by define_env_space().
#'
#' @return A ggplot object showing planned vs realised environmental gradients.
#' @export
plot_env_space <- function(env_space) {

  if (!inherits(env_space, "env_space")) {
    stop("Input must be an 'env_space' object.")
  }

  df <- env_space$realised
  vars <- env_space$variables

  # Melt to long format
  df_long <- reshape2::melt(df, id.vars = "group_id",
                            measure.vars = vars)

  ggplot2::ggplot(df_long, ggplot2::aes(
    x = group_id, y = value, color = variable, group = variable
  )) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "Realised Environmental Space",
      x = "Group / Site",
      y = "Environmental Value",
      color = "Variable"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}

#' Plot individual-level environmental distributions
#' @export
plot_env_individual <- function(ind_env) {
  library(ggplot2)
  var_names <- setdiff(colnames(ind_env), "group_id")

  # Facet histograms of each variable
  df_long <- tidyr::pivot_longer(ind_env, all_of(var_names))

  p <- ggplot(df_long, aes(value, fill = factor(group_id))) +
    geom_density(alpha = 0.4) +
    facet_wrap(~name, scales = "free") +
    theme_bw() +
    labs(title = "Individual-level Environmental Distributions",
         fill = "Group ID")

  return(p)
}

#' Plot host model diagnostics
#' @export
plot_host_model <- function(host_model) {
  if (host_model$type == "clone") {
    text <- "Clone system: Host variance = 0\n(No host-driven noise)"
  } else {
    text <- sprintf(
      "Natural population\nHost variance fraction = %.2f\nStructure = %s",
      host_model$var_host, host_model$pop_structure
    )
  }

  df <- data.frame(x = 1, y = 1, label = text)

  ggplot(df, aes(x, y, label = label)) +
    geom_text(size = 6) +
    theme_void()
}

#' Plot holobiont component diagnostics (unique feature)
#' Shows structure, correlations, PCA, env-effect signatures
#' @export
plot_holobiont_components <- function(components_df, env = NULL, labels = NULL, color_by = NULL) {
  library(ggplot2)

  comp_names <- colnames(components_df)

  # PCA overview
  pc <- prcomp(components_df, scale. = TRUE)
  pc_df <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2])

  # Add labels if provided
  if (!is.null(labels)) {
    pc_df$label <- labels
  } else {
    pc_df$label <- rownames(components_df)
  }

  # Add color grouping if provided
  if (!is.null(color_by)) {
    pc_df$group <- color_by
    p1 <- ggplot(pc_df, aes(PC1, PC2, color = group, label = label)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_text(size = 2.5, hjust = -0.1, vjust = 0.5, show.legend = FALSE) +
      theme_bw() +
      labs(title = "Holobiont Component PCA", color = "Group")
  } else {
    p1 <- ggplot(pc_df, aes(PC1, PC2, label = label)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_text(size = 2.5, hjust = -0.1, vjust = 0.5) +
      theme_bw() +
      labs(title = "Holobiont Component PCA")
  }

  # Component correlation heatmap
  corr <- stats::cor(components_df)
  corr_df <- as.data.frame(as.table(corr))
  colnames(corr_df) <- c("Var1","Var2","value")
  p2 <- ggplot(corr_df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme_bw() +
    labs(title = "Holobiont Component Correlation Matrix")

  # Env effects (if env provided)
  if (!is.null(env)) {
    env_df <- as.data.frame(env)
    names(env_df) <- colnames(env)
    env_cor <- cor(env_df, components_df)
    env_cor_df <- as.data.frame(as.table(env_cor))
    colnames(env_cor_df) <- c("Environment","Component","Correlation")
    p3 <- ggplot(env_cor_df, aes(Environment, Component, fill = Correlation)) +
      geom_tile() +
      scale_fill_gradient2() +
      theme_bw() +
      labs(title = "Env–Holobiont Correlation Matrix")
    return(list(PCA = p1, CorMatrix = p2, EnvEffects = p3))
  }
  return(list(PCA = p1, CorMatrix = p2))
}

#' Plot power landscape (design exploration)
#' @export
plot_power_landscape <- function(power_results) {
  library(ggplot2)

  p <- ggplot(power_results,
              aes(n_groups, n_per_group, fill = power_env)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    labs(title = "Power Landscape",
         x = "Number of groups",
         y = "Individuals per group",
         fill = "Power")

  return(p)
}

#' Plot variance decomposition for holobiont components
#' A unique + important diagnostic!
#' @export
plot_variance_decomposition <- function(power_results, components) {

  # Build variance table
  var_tab <- lapply(names(components), function(comp) {
    r2_env <- components[[comp]]$r2_env
    r2_env_total <- if (length(r2_env) == 1) r2_env else sum(r2_env)

    data.frame(
      component = comp,
      var_env = r2_env_total,
      var_host = power_results$var_host[1],
      var_resid = 1 - r2_env_total - power_results$var_host[1]
    )
  })

  var_tab <- do.call(rbind, var_tab)

  df_long <- tidyr::pivot_longer(
    var_tab,
    cols = c("var_env","var_host","var_resid"),
    names_to = "source",
    values_to = "value"
  )

  ggplot(df_long, aes(component, value, fill = source)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(title = "Variance Decomposition per Component",
         y = "Proportion of Variance")
}
