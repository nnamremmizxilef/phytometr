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
