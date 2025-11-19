#' Holobiont component diagnostics (PCA + correlations)
#'
#' @param components_df Data frame from simulate_holobiont_components().
#'   Must contain an `id` column; may contain `group_id`.
#' @param env Optional environmental data (data frame or matrix) with the
#'   same number of rows as `components_df`. Used to compute
#'   environment-component correlations.
#' @param labels Optional vector of labels for points in the PCA plot.
#'   Defaults to `components_df$id`.
#' @param color_by Optional column name in `components_df` or vector
#'   of length `nrow(components_df)` used to colour points.
#'
#' @return A list of ggplot objects: `PCA`, `CorMatrix`, and
#'   optionally `EnvEffects` if `env` is supplied.
#' @export
plot_holobiont_components <- function(components_df,
                                      env      = NULL,
                                      labels   = NULL,
                                      color_by = NULL) {

  components_df <- as.data.frame(components_df)

  if (!("id" %in% names(components_df))) {
    stop("`components_df` must contain an 'id' column.")
  }

  # Only true component columns (exclude id / group_id) for PCA & cor
  comp_only <- components_df[ ,
    setdiff(names(components_df), c("id", "group_id")),
    drop = FALSE
  ]

  ## ---------------------------------------------------------------
  ## 1) PCA with group outlines
  ## ---------------------------------------------------------------

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

  # colours (always discrete factor)
  if (is.null(color_by) && "group_id" %in% names(components_df)) {

    scores$group <- factor(components_df$group_id)
    colour_lab   <- "group_id"

  } else if (is.null(color_by)) {

    scores$group <- factor(components_df$id)
    colour_lab   <- "id"

  } else if (is.character(color_by) &&
             length(color_by) == 1L &&
             color_by %in% names(components_df)) {

    scores$group <- factor(components_df[[color_by]])
    colour_lab   <- color_by

  } else if (length(color_by) == nrow(components_df)) {

    scores$group <- factor(color_by)
    colour_lab   <- "group"

  } else {
    stop("`color_by` must be NULL, a column name in components_df, or a vector of length nrow(components_df).")
  }

  p_pca <-
    ggplot2::ggplot(
      scores,
      ggplot2::aes(x = PC1, y = PC2, colour = group, label = label)
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = group),
      level       = 0.68,
      linewidth   = 0.5,
      alpha       = 0.7,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    ggplot2::labs(
      title  = "Holobiont Component PCA",
      colour = colour_lab
    ) +
    ggplot2::theme_minimal(base_size = 14)

  ## ---------------------------------------------------------------
  ## 2) Component-component correlation heatmap (with r & p)
  ## ---------------------------------------------------------------

  cor_mat <- stats::cor(comp_only, use = "pairwise.complete.obs")

  # p-values via cor.test
  k <- ncol(comp_only)
  p_mat <- matrix(NA_real_, nrow = k, ncol = k,
                  dimnames = list(colnames(comp_only), colnames(comp_only)))

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      x <- comp_only[, i]
      y <- comp_only[, j]
      ct <- suppressWarnings(stats::cor.test(x, y))
      p_mat[i, j] <- ct$p.value
    }
  }

  cor_long <- reshape2::melt(
    cor_mat,
    varnames   = c("Var1", "Var2"),
    value.name = "r"
  )
  p_long <- reshape2::melt(
    p_mat,
    varnames   = c("Var1", "Var2"),
    value.name = "p"
  )

  cor_long <- merge(cor_long, p_long, by = c("Var1", "Var2"))
  cor_long$label <- sprintf("r=%.2f\np=%.3f", cor_long$r, cor_long$p)

  p_cor <-
    ggplot2::ggplot(
      cor_long,
      ggplot2::aes(x = Var1, y = Var2, fill = r)
    ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      name   = "Correlation",
      low    = "#67001f",
      mid    = "white",
      high   = "#053061",
      limits = c(-1, 1)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 3
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

  ## ---------------------------------------------------------------
  ## 3) Environment-component correlation heatmap (with r & p)
  ## ---------------------------------------------------------------

  if (!is.null(env)) {
    env_df <- as.data.frame(env)

    if (nrow(env_df) != nrow(components_df)) {
      stop("`env` must have the same number of rows as `components_df`.")
    }

    env_cor <- stats::cor(env_df, comp_only, use = "pairwise.complete.obs")

    m <- ncol(env_df)
    p_env_mat <- matrix(
      NA_real_, nrow = m, ncol = k,
      dimnames = list(colnames(env_df), colnames(comp_only))
    )

    for (i in seq_len(m)) {
      for (j in seq_len(k)) {
        x <- env_df[, i]
        y <- comp_only[, j]
        ct <- suppressWarnings(stats::cor.test(x, y))
        p_env_mat[i, j] <- ct$p.value
      }
    }

    env_cor_long <- reshape2::melt(
      env_cor,
      varnames   = c("Environment", "Component"),
      value.name = "r"
    )
    p_env_long <- reshape2::melt(
      p_env_mat,
      varnames   = c("Environment", "Component"),
      value.name = "p"
    )

    env_cor_long <- merge(env_cor_long, p_env_long,
                          by = c("Environment", "Component"))
    env_cor_long$label <- sprintf("r=%.2f\np=%.3f",
                                  env_cor_long$r, env_cor_long$p)

    p_env <-
      ggplot2::ggplot(
        env_cor_long,
        ggplot2::aes(x = Environment, y = Component, fill = r)
      ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        name   = "Correlation",
        low    = "#67001f",
        mid    = "white",
        high   = "#053061",
        limits = c(-1, 1)
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        size = 3
      ) +
      ggplot2::labs(
        title = "Environment-Holobiont Correlation Matrix",
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
