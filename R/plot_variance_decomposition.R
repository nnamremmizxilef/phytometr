#' Plot variance decomposition for holobiont components
#'
#' This is a conceptual diagnostic: given the R² you specified for
#' environment–component relationships, it shows how much variance is
#' attributed to environment, host, and residual noise for each component.
#'
#' In the current Phase 1 implementation, host variance is not yet simulated,
#' so the default is \code{var_host = 0}. You can change this to explore
#' "what-if" scenarios (e.g. \code{var_host = 0.2}).
#'
#' @param power_results A data frame returned by
#'   \code{simulate_holobiont_power()}. Only used for its component names
#'   (the power values themselves are not directly used here).
#' @param components The component specification list used for simulation
#'   (the same object you pass to \code{simulate_holobiont_components()}).
#'   Each component must have an element \code{r2_env}.
#' @param var_host Numeric scalar giving the fraction of variance you want
#'   to attribute to host effects (default 0 for clone systems).
#'
#' @return A stacked barplot (ggplot object) showing variance fractions
#'   for environment, host, and residual per component.
#'
#' @export
plot_variance_decomposition <- function(power_results,
                                        components,
                                        var_host = 0) {

  if (!is.list(components) || is.null(names(components))) {
    stop("`components` must be the named list used for simulation (with r2_env entries).")
  }

  var_tab <- lapply(names(components), function(comp) {

    r2_env_vec   <- components[[comp]]$r2_env
    r2_env_total <- if (length(r2_env_vec) == 1) r2_env_vec else sum(r2_env_vec)

    if (r2_env_total + var_host > 1) {
      warning(sprintf(
        "For component '%s', env + host variance exceeds 1 (%.2f); truncating residual to 0.",
        comp, r2_env_total + var_host
      ))
      var_resid <- 0
    } else {
      var_resid <- 1 - r2_env_total - var_host
    }

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
    ggplot2::aes(x = component, y = value, fill = source)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_fill_brewer(
      palette = "Set2",
      name    = "Source",
      labels  = c(
        var_env   = "Environment",
        var_host  = "Host",
        var_resid = "Residual"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Variance Decomposition per Component",
      x     = "Component",
      y     = "Proportion of variance"
    )
}
