#' Holobiont component diagnostics (PCA + correlations)
#'
#' @param components_df Data frame from simulate_holobiont_components().
#'   Must contain an \code{id} column; may contain \code{group_id}.
#' @param env Optional environmental data (data frame or matrix) with the
#'   same number of rows as \code{components_df}. Used to compute
#'   env–component correlations.
#' @param labels Optional vector of labels for points in the PCA plot.
#'   Defaults to \code{components_df$id}.
#' @param color_by Optional column name in \code{components_df} or vector
#'   of length \code{nrow(components_df)} used to colour points.
#'
#' @return A list of ggplot objects: \code{PCA}, \code{CorMatrix}, and
#'   optionally \code{EnvEffects} if \code{env} is supplied.
#' @export
plot_holobiont_components <- function(components_df,
                                      env      = NULL,
                                      labels   = NULL,
                                      color_by = NULL) {

  components_df <- as.data.frame(components_df)

  if (!("id" %in% names(components_df))) {
    stop("`components_df` must contain an 'id' column.")
  }

  # Use only component columns (exclude id and group_id) for PCA/cor
  comp_only <- components_df[ , setdiff(names(components_df), c("id", "group_id")), drop = FALSE]

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

  # colours: if color_by is NULL and group_id exists, use group_id
  if (is.null(color_by) && "group_id" %in% names(components_df)) {

    scores$group <- factor(components_df$group_id)
    colour_lab   <- "Group"

  } else if (is.null(color_by)) {

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

  ## --- Component correlation heatmap (without id/group_id) ------------

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
