#' Run single replicate of power simulation
#'
#' @param n Sample size for this replicate
#' @param params Simulation parameters
#' @return List with power metrics
#' @keywords internal
run_single_replicate <- function(n, params) {
  
  tryCatch({
    # Create holobiont data structure
    holobiont <- create_holobiont_data(n, host_name = "simulated_host")
    
    # Simulate host genotypes
    holobiont <- simulate_host_genetics(holobiont, params)
    
    # Simulate environment
    holobiont <- simulate_environment_data(holobiont, params)
    
    # Simulate response variable
    holobiont <- simulate_response(holobiont, params)
    
    # Run GEA analysis
    gea_results <- run_gea_analysis(holobiont, params)
    
    # Calculate power metrics
    metrics <- calc_power_metrics(gea_results, holobiont$causal_loci, params)
    
    return(metrics)
    
  }, error = function(e) {
    warning("Replicate failed: ", e$message)
    return(list(power = NA, fdr = NA, n_detected = NA))
  })
}

#' Calculate power metrics from GEA results
#'
#' @param gea_results Results from run_gea_analysis
#' @param true_causal Integer vector of true causal loci positions
#' @param params Simulation parameters
#' @return List with power, FDR, and number detected
#' @keywords internal
calc_power_metrics <- function(gea_results, true_causal, params) {
  
  pvals <- gea_results$pvals
  
  # Apply multiple testing correction
  if (params$correction_method == "bonferroni") {
    pvals_adj <- p.adjust(pvals, method = "bonferroni")
  } else if (params$correction_method == "fdr") {
    pvals_adj <- p.adjust(pvals, method = "fdr")
  } else {
    pvals_adj <- pvals  # No correction
  }
  
  # Identify detected loci
  detected <- which(pvals_adj < params$alpha)
  
  # Calculate metrics
  if (length(detected) == 0) {
    tpr <- 0
    fdr <- 0
  } else {
    # True positive rate (power)
    tpr <- sum(detected %in% true_causal) / length(true_causal)
    
    # False discovery rate
    n_false_positives <- sum(!(detected %in% true_causal))
    fdr <- n_false_positives / length(detected)
  }
  
  return(list(
    power = tpr,
    fdr = fdr,
    n_detected = length(detected),
    n_true_positives = sum(detected %in% true_causal),
    n_false_positives = sum(!(detected %in% true_causal))
  ))
}

#' Aggregate results across replicates
#'
#' @param replicates List of replicate results
#' @param params Simulation parameters
#' @return List with aggregated metrics
#' @keywords internal
aggregate_replicates <- function(replicates, params) {
  
  # Extract metrics from each replicate
  power_vals <- sapply(replicates, function(x) x$power)
  fdr_vals <- sapply(replicates, function(x) x$fdr)
  n_detected_vals <- sapply(replicates, function(x) x$n_detected)
  
  # Remove NAs from failed replicates
  power_vals <- power_vals[!is.na(power_vals)]
  fdr_vals <- fdr_vals[!is.na(fdr_vals)]
  n_detected_vals <- n_detected_vals[!is.na(n_detected_vals)]
  
  # Calculate summary statistics
  list(
    mean_power = mean(power_vals),
    se_power = sd(power_vals) / sqrt(length(power_vals)),
    median_power = median(power_vals),
    mean_fdr = mean(fdr_vals),
    se_fdr = sd(fdr_vals) / sqrt(length(fdr_vals)),
    mean_detected = mean(n_detected_vals),
    n_replicates_succeeded = length(power_vals)
  )
}

#' Format power results for output
#'
#' @param power_list List of power results for each sample size
#' @param sample_sizes Vector of sample sizes tested
#' @param params Simulation parameters
#' @return Data frame of class 'gea_power_analysis'
#' @keywords internal
format_power_results <- function(power_list, sample_sizes, params) {
  
  results_df <- data.frame(
    n_individuals = sample_sizes,
    power = sapply(power_list, function(x) x$mean_power),
    power_se = sapply(power_list, function(x) x$se_power),
    power_median = sapply(power_list, function(x) x$median_power),
    fdr = sapply(power_list, function(x) x$mean_fdr),
    fdr_se = sapply(power_list, function(x) x$se_fdr),
    n_detected = sapply(power_list, function(x) x$mean_detected),
    n_replicates = sapply(power_list, function(x) x$n_replicates_succeeded)
  )
  
  # Add parameters as attributes
  attr(results_df, "parameters") <- params
  attr(results_df, "target_power") <- params$target_power
  
  # Set class
  class(results_df) <- c("gea_power_analysis", "data.frame")
  
  return(results_df)
}