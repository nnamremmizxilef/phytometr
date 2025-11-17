#' Plot power analysis results
#'
#' @param x Object of class 'gea_power_analysis'
#' @param add_target_line Logical. Add horizontal line at target power? (default: TRUE)
#' @param add_error_bars Logical. Add error bars for power estimates? (default: TRUE)
#' @param ... Additional arguments passed to ggplot
#' @return ggplot object
#' @export
plot.gea_power_analysis <- function(x, add_target_line = TRUE, 
                                   add_error_bars = TRUE, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  params <- attr(x, "parameters")
  target_power <- attr(x, "target_power")
  
  # Base plot
  p <- ggplot2::ggplot(x, ggplot2::aes(x = n_individuals, y = power)) +
    ggplot2::geom_line(linewidth = 1, color = "#2c7fb8") +
    ggplot2::geom_point(size = 3, color = "#2c7fb8")
  
  # Add error bars if requested
  if (add_error_bars) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = power - power_se, ymax = power + power_se),
      width = (max(x$n_individuals) - min(x$n_individuals)) * 0.02,
      color = "#2c7fb8"
    )
  }
  
  # Add target power line
  if (add_target_line) {
    p <- p + ggplot2::geom_hline(
      yintercept = target_power,
      linetype = "dashed",
      color = "#d7301f",
      linewidth = 0.8
    )
  }
  
  # Labels and theme
  p <- p +
    ggplot2::labs(
      x = "Sample size (n individuals)",
      y = "Statistical power",
      title = sprintf("GEA Power Analysis: %s", params$response_variable),
      subtitle = sprintf(
        "Effect size: %.2f | Method: %s | %d causal loci",
        params$effect_size, params$gea_method, params$n_causal
      )
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

#' Plot power and FDR together
#'
#' @param power_results Object of class 'gea_power_analysis'
#' @return ggplot object
#' @export
plot_power_fdr <- function(power_results) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  # Reshape data for plotting
  plot_data <- data.frame(
    n_individuals = rep(power_results$n_individuals, 2),
    value = c(power_results$power, power_results$fdr),
    se = c(power_results$power_se, power_results$fdr_se),
    metric = rep(c("Power", "FDR"), each = nrow(power_results))
  )
  
  params <- attr(power_results, "parameters")
  
  p <- ggplot2::ggplot(plot_data, 
                       ggplot2::aes(x = n_individuals, y = value, 
                                   color = metric, fill = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = value - se, ymax = value + se),
                        alpha = 0.2, color = NA) +
    ggplot2::scale_color_manual(values = c("Power" = "#2c7fb8", "FDR" = "#d7301f")) +
    ggplot2::scale_fill_manual(values = c("Power" = "#2c7fb8", "FDR" = "#d7301f")) +
    ggplot2::labs(
      x = "Sample size (n individuals)",
      y = "Value",
      title = "Power and False Discovery Rate",
      subtitle = sprintf("Effect size: %.2f | Method: %s", 
                        params$effect_size, params$gea_method),
      color = "Metric",
      fill = "Metric"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Plot number of loci detected
#'
#' @param power_results Object of class 'gea_power_analysis'
#' @return ggplot object
#' @export
plot_n_detected <- function(power_results) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  params <- attr(power_results, "parameters")
  
  p <- ggplot2::ggplot(power_results, 
                       ggplot2::aes(x = n_individuals, y = n_detected)) +
    ggplot2::geom_line(linewidth = 1, color = "#41ab5d") +
    ggplot2::geom_point(size = 3, color = "#41ab5d") +
    ggplot2::geom_hline(yintercept = params$n_causal, 
                       linetype = "dashed", color = "#d7301f") +
    ggplot2::annotate("text", 
                     x = min(power_results$n_individuals), 
                     y = params$n_causal,
                     label = sprintf("True causal loci (n=%d)", params$n_causal),
                     hjust = 0, vjust = -0.5, size = 3.5) +
    ggplot2::labs(
      x = "Sample size (n individuals)",
      y = "Mean number of loci detected",
      title = "Number of Loci Detected vs. Sample Size"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
  
  return(p)
}