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
