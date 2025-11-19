#' Plot variance decomposition for holobiont components
#'
#' Shows, for each component, the relative contributions of:
#' \itemize{
#'   \item environmental variance (sum of R² across all env variables),
#'   \item host variance (from the design/power results),
#'   \item residual variance.
#' }
#'
#' @param power_results A data.frame returned by
#'   \code{simulate_holobiont_power()} containing, at minimum,
#'   columns \code{n_groups}, \code{n_per_group}, \code{power_env},
#'   and an attribute or column giving \code{var_host} (host variance
#'   fraction). For Phase 1 clone systems this is typically zero.
#' @param components The same \code{components} list used in
#'   \code{simulate_holobiont_power()}, i.e. a named list of component
#'   specifications with \code{r2_env}.
#'
#' @return A ggplot object showing stacked bars of variance components.
#' @export
plot_variance_decomposition <- function(power_results, components) {
  if (!is.list(components) || length(components) == 0L) {
    stop("'components' must be a non-empty list.")
  }
  if (!is.data.frame(power_results)) {
    stop("'power_results' must be a data.frame.")
  }

  # For now, take the first design's var_host if present, else 0
  var_host <- if ("var_host" %in% colnames(power_results)) {
    power_results$var_host[1]
  } else if (!is.null(attr(power_results, "var_host"))) {
    attr(power_results, "var_host")
  } else {
    0
  }

  var_tab <- lapply(names(components), function(comp) {
    spec <- components[[comp]]
    r2_env <- spec$r2_env
    if (any(r2_env < 0)) {
      stop("Component '", comp,
           "': 'r2_env' must be >= 0 for all variables.")
    }
    var_env <- sum(r2_env)
    var_resid <- 1 - var_env - var_host
    if (var_resid < 0) var_resid <- 0

    data.frame(
      component = comp,
      var_env   = var_env,
      var_host  = var_host,
      var_resid = var_resid
    )
  })

  var_tab <- do.call(rbind, var_tab)

  df_long <- tidyr::pivot_longer(
    var_tab,
    cols = c("var_env", "var_host", "var_resid"),
    names_to = "source",
    values_to = "value"
  )

  ggplot2::ggplot(df_long, ggplot2::aes(component, value, fill = source)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Variance Decomposition per Component",
      y      = "Proportion of Variance",
      x      = "Holobiont component",
      fill   = "Source"
    )
}
