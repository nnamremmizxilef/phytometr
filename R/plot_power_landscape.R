#' Plot power landscape over design grid
#'
#' @description
#' Heatmap of detection power over the design grid
#' (number of groups × individuals per group). By default it
#' shows *all* design combinations that were requested in
#' \code{simulate_holobiont_power()}, even if some of them
#' were not successfully simulated.
#'
#' @param power_results Data frame returned by simulate_holobiont_power().
#'   Must contain columns \code{n_groups}, \code{n_per_group}, and
#'   \code{power_env}.
#'
#' @return A ggplot object.
#' @export
plot_power_landscape <- function(power_results) {
  if (!all(c("n_groups", "n_per_group", "power_env") %in%
           colnames(power_results))) {
    stop("power_results must contain 'n_groups', 'n_per_group', and 'power_env'.")
  }

  # Try to recover the *full* design grid from attributes
  n_groups_all <- attr(power_results, "n_groups_use")
  n_per_all    <- attr(power_results, "n_per_group")

  # Fallback: use unique values in the data if attributes are missing
  if (is.null(n_groups_all)) n_groups_all <- sort(unique(power_results$n_groups))
  if (is.null(n_per_all))    n_per_all    <- sort(unique(power_results$n_per_group))

  # Build full grid of requested designs
  full_grid <- expand.grid(
    n_groups    = n_groups_all,
    n_per_group = n_per_all
  )

  # Merge power results onto the full grid
  plot_dat <- merge(
    full_grid,
    power_results[, c("n_groups", "n_per_group", "power_env")],
    by = c("n_groups", "n_per_group"),
    all.x = TRUE
  )

  # Ensure numeric ordering on axes
  plot_dat$n_groups    <- as.factor(plot_dat$n_groups)
  plot_dat$n_per_group <- as.factor(plot_dat$n_per_group)

  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = n_groups, y = n_per_group, fill = power_env)
  ) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::scale_fill_viridis_c(
      option  = "D",
      direction = 1,
      limits = c(0, 1),
      na.value = "grey95",
      name = "Power"
    ) +
    ggplot2::labs(
      title = "Power landscape over design grid",
      x = "Number of groups / conditions",
      y = "Individuals per group"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}
