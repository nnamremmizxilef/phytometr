#' Validate simulation parameters
#'
#' @param ... Named simulation parameters
#' @return List of validated parameters
#' @keywords internal
validate_simulation_parameters <- function(...) {
  
  params <- list(...)
  
  # Check sample sizes - stop if negative, warn if small
  if (any(params$n_individuals < 1)) {
    stop("n_individuals must be positive")
  }
  
  if (any(params$n_individuals < 10)) {
    warning("Sample sizes < 10 may produce unreliable results")
  }
  
  # Check genetic parameters
  if (params$n_causal > params$n_loci) {
    stop("n_causal cannot exceed n_loci")
  }
  
  if (params$n_causal < 1) {
    stop("n_causal must be at least 1")
  }
  
  # Check effect size BEFORE running simulations
  if ("effect_size" %in% names(params)) {
    if (params$effect_size < 0) {
      stop("effect_size must be non-negative")
    }
  }
  
  # Check heritability parameters
  total_var <- params$holobiont_heritability + params$env_effect_on_holobiont
  if (total_var > 1) {
    warning(sprintf("holobiont_heritability (%.2f) + env_effect_on_holobiont (%.2f) > 1. Rescaling.",
                   params$holobiont_heritability, params$env_effect_on_holobiont))
    scale_factor <- 0.9 / total_var
    params$holobiont_heritability <- params$holobiont_heritability * scale_factor
    params$env_effect_on_holobiont <- params$env_effect_on_holobiont * scale_factor
  }
  
  # Check alpha
  if ("alpha" %in% names(params)) {
    if (params$alpha <= 0 || params$alpha >= 1) {
      stop("alpha must be between 0 and 1")
    }
  }
  
  # Check replicates
  if (params$n_replicates < 10) {
    warning("n_replicates < 10 may produce unreliable power estimates")
  }
  
  return(params)
}

#' Calculate required sample size for target power
#'
#' @param power_results Object of class 'gea_power_analysis'
#' @param target_power Target power threshold (default: 0.8)
#' @return Integer sample size needed, or NA if target not achieved
#' @export
calculate_required_n <- function(power_results, target_power = 0.8) {
  
  if (!inherits(power_results, "gea_power_analysis")) {
    stop("Input must be a gea_power_analysis object")
  }
  
  # Find first sample size achieving target power
  achieved <- power_results$power >= target_power
  
  if (!any(achieved)) {
    warning(sprintf("Target power of %.2f not achieved in tested sample sizes", 
                   target_power))
    return(NA)
  }
  
  required_n <- min(power_results$n_individuals[achieved])
  
  message(sprintf("Required sample size for %.0f%% power: n = %d", 
                 target_power * 100, required_n))
  
  return(required_n)
}

#' Summary method for gea_power_analysis
#' @param object A gea_power_analysis object
#' @param ... Additional arguments (ignored)
#' @export
summary.gea_power_analysis <- function(object, ...) {
  
  params <- attr(object, "parameters")
  target <- attr(object, "target_power")
  
  cat("GEA Power Analysis Summary\n")
  cat("==========================\n\n")
  
  cat("Simulation parameters:\n")
  cat(sprintf("  Number of loci: %d\n", params$n_loci))
  cat(sprintf("  Causal loci: %d\n", params$n_causal))
  cat(sprintf("  Effect size: %.3f\n", params$effect_size))
  cat(sprintf("  Population structure: %s\n", params$pop_structure))
  cat(sprintf("  Response variable: %s\n", params$response_variable))
  cat(sprintf("  GEA method: %s\n", params$gea_method))
  cat(sprintf("  Replicates per sample size: %d\n\n", 
             object$n_replicates[1]))
  
  cat("Power results:\n")
  print(object[, c("n_individuals", "power", "fdr", "n_detected")], 
        row.names = FALSE)
  
  cat("\n")
  required_n <- calculate_required_n(object, target)
  
  if (!is.na(required_n)) {
    cat(sprintf("\nRequired sample size for %.0f%% power: %d individuals\n", 
               target * 100, required_n))
  }
  
  invisible(object)
}