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
