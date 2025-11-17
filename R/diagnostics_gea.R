#' GEA-specific diagnostic functions
#'
#' @name diagnostics_gea
#' @keywords internal
NULL

#' Assess GEA assumptions
#'
#' Check whether data meet assumptions for genotype-environment association analysis.
#' Tests for genotype-environment independence, population structure confounding,
#' and signal-to-noise ratio.
#'
#' @param holobiont A holobiont_data object
#' @param env_variable Character. Name of environmental variable to test (default: "main")
#' @param method Character. GEA method to assess for: "rda", "lfmm", "baypass" (default: "rda")
#' @param verbose Logical. Print detailed results? (default: TRUE)
#' @return List of assumption checks with class 'gea_assumptions'
#'
#' @examples
#' \dontrun{
#' assumptions <- assess_gea_assumptions(holobiont, env_variable = "temperature")
#' print(assumptions)
#' plot(assumptions)
#' }
#'
#' @export
assess_gea_assumptions <- function(holobiont, 
                                  env_variable = "main",
                                  method = c("rda", "lfmm", "baypass"),
                                  verbose = TRUE) {
  
  method <- match.arg(method)
  
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  
  results <- list()
  
  # 1. Check genotype-environment independence
  results$ge_independence <- test_ge_independence(holobiont, env_variable)
  
  # 2. Assess population structure confounding
  results$structure_confounding <- assess_structure_confounding(holobiont, env_variable)
  
  # 3. Estimate signal-to-noise ratio
  results$signal_noise <- estimate_signal_noise(holobiont, env_variable)
  
  # 4. Check for outliers
  results$outliers <- detect_outliers(holobiont, env_variable)
  
  # 5. Assess environmental gradient quality
  results$gradient_quality <- assess_gradient_quality(holobiont, env_variable)
  
  # 6. Method-specific checks
  if (method == "rda") {
    results$method_specific <- check_rda_assumptions(holobiont, env_variable)
  } else if (method == "lfmm") {
    results$method_specific <- check_lfmm_assumptions(holobiont, env_variable)
  }
  
  # Overall assessment
  results$overall_suitable <- assess_overall_suitability(results)
  
  class(results) <- "gea_assumptions"
  
  if (verbose) {
    print(results)
  }
  
  invisible(results)
}

#' Test genotype-environment independence
#'
#' Tests whether genotypes are independent of environment (no G-E correlation
#' due to structure alone).
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with test results
#' @keywords internal
test_ge_independence <- function(holobiont, env_variable) {
  
  snps <- holobiont$host_genotypes
  env <- holobiont$environment[[env_variable]]
  
  if (is.null(snps) || is.null(env)) {
    return(list(status = "Insufficient data"))
  }
  
  # Test correlation between genetic PCs and environment
  # This detects if environment is confounded with population structure
  
  # PCA on genotypes
  pca_result <- fast_pca(snps, n_pcs = 5)
  
  # Test correlation with environment
  pc_env_cors <- sapply(1:ncol(pca_result$pcs), function(i) {
    cor.test(pca_result$pcs[, i], env)$estimate
  })
  
  # Test significance
  pc_env_pvals <- sapply(1:ncol(pca_result$pcs), function(i) {
    cor.test(pca_result$pcs[, i], env)$p.value
  })
  
  max_cor <- max(abs(pc_env_cors))
  significant <- any(pc_env_pvals < 0.05)
  
  list(
    pc_correlations = pc_env_cors,
    pc_pvalues = pc_env_pvals,
    max_correlation = max_cor,
    significant_confounding = significant,
    warning = significant && max_cor > 0.5,
    message = if (significant && max_cor > 0.5) {
      "WARNING: Strong correlation between genetic structure and environment detected"
    } else if (significant) {
      "Moderate correlation between genetic structure and environment"
    } else {
      "No significant genotype-environment confounding detected"
    }
  )
}

#' Fast PCA implementation
#'
#' @param snps SNP matrix
#' @param n_pcs Number of PCs to compute
#' @return List with PCs and variance explained
#' @keywords internal
fast_pca <- function(snps, n_pcs = 5) {
  
  # Remove SNPs with no variation
  snp_var <- apply(snps, 2, var, na.rm = TRUE)
  snps_filtered <- snps[, snp_var > 0 & !is.na(snp_var), drop = FALSE]
  
  # Impute missing with column means
  col_means <- colMeans(snps_filtered, na.rm = TRUE)
  for (j in 1:ncol(snps_filtered)) {
    missing <- is.na(snps_filtered[, j])
    if (any(missing)) {
      snps_filtered[missing, j] <- col_means[j]
    }
  }
  
  # Center and scale
  snps_scaled <- scale(snps_filtered)
  
  # SVD (faster than eigen for tall matrices)
  n_pcs <- min(n_pcs, nrow(snps_scaled) - 1, ncol(snps_scaled))
  
  if (ncol(snps_scaled) > 1000) {
    # Use randomized SVD for large matrices
    svd_result <- tryCatch({
      # If available, use irlba package
      if (requireNamespace("irlba", quietly = TRUE)) {
        irlba::irlba(snps_scaled, nv = n_pcs)
      } else {
        svd(snps_scaled)
      }
    }, error = function(e) {
      svd(snps_scaled)
    })
  } else {
    svd_result <- svd(snps_scaled)
  }
  
  # Extract PCs
  pcs <- svd_result$u[, 1:n_pcs, drop = FALSE]
  
  # Variance explained
  var_explained <- svd_result$d[1:n_pcs]^2 / sum(svd_result$d^2)
  
  list(
    pcs = pcs,
    variance_explained = var_explained
  )
}

#' Assess population structure confounding
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with confounding assessment
#' @keywords internal
assess_structure_confounding <- function(holobiont, env_variable) {
  
  if (is.null(holobiont$population)) {
    return(list(status = "No population structure information"))
  }
  
  env <- holobiont$environment[[env_variable]]
  pops <- holobiont$population
  
  # Test if environment varies between populations
  if (length(unique(pops)) > 1) {
    # ANOVA or Kruskal-Wallis
    if (is.numeric(env)) {
      test <- kruskal.test(env ~ factor(pops))
      
      # Calculate effect size (eta-squared)
      pop_means <- tapply(env, pops, mean, na.rm = TRUE)
      between_var <- var(pop_means, na.rm = TRUE)
      total_var <- var(env, na.rm = TRUE)
      eta_sq <- between_var / total_var
      
      list(
        test = "Kruskal-Wallis",
        statistic = test$statistic,
        p_value = test$p.value,
        effect_size = eta_sq,
        confounded = test$p.value < 0.05,
        message = if (test$p.value < 0.05) {
          sprintf("Environment varies significantly between populations (p = %.4f, eta^2 = %.3f)",
                 test$p.value, eta_sq)
        } else {
          "Environment does not vary significantly between populations"
        }
      )
    }
  } else {
    list(status = "Only one population")
  }
}

#' Estimate signal-to-noise ratio
#'
#' Estimates the ratio of genetic signal to environmental noise in the data.
#' Higher values indicate better power for detection.
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with SNR estimates
#' @keywords internal
estimate_signal_noise <- function(holobiont, env_variable) {
  
  snps <- holobiont$host_genotypes
  env <- holobiont$environment[[env_variable]]
  
  if (is.null(snps) || is.null(env)) {
    return(list(status = "Insufficient data"))
  }
  
  # Sample subset of SNPs for efficiency
  if (ncol(snps) > 1000) {
    sample_snps <- sample(ncol(snps), 1000)
    snps_sample <- snps[, sample_snps]
  } else {
    snps_sample <- snps
  }
  
  # Calculate correlation of each SNP with environment
  correlations <- apply(snps_sample, 2, function(snp) {
    valid <- !is.na(snp) & !is.na(env)
    if (sum(valid) < 10) return(NA)
    cor(snp[valid], env[valid], method = "spearman")
  })
  
  # Signal: magnitude of correlations
  signal <- median(abs(correlations), na.rm = TRUE)
  
  # Noise: expected correlation under null
  # For independent data, expected |r| ~ 1/sqrt(n)
  n <- length(env)
  expected_null <- 1 / sqrt(n)
  
  snr <- signal / expected_null
  
  list(
    median_abs_correlation = signal,
    expected_null = expected_null,
    signal_to_noise_ratio = snr,
    adequate_snr = snr > 1.5,
    message = if (snr > 2) {
      "Strong signal-to-noise ratio"
    } else if (snr > 1.5) {
      "Adequate signal-to-noise ratio"
    } else {
      "Low signal-to-noise ratio - consider increasing sample size"
    }
  )
}

#' Detect outliers in holobiont data
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with outlier information
#' @keywords internal
detect_outliers <- function(holobiont, env_variable) {
  
  env <- holobiont$environment[[env_variable]]
  
  if (is.null(env)) {
    return(list(status = "No environmental data"))
  }
  
  # Environmental outliers (3 SD rule)
  env_mean <- mean(env, na.rm = TRUE)
  env_sd <- sd(env, na.rm = TRUE)
  env_outliers <- abs(env - env_mean) > 3 * env_sd
  
  results <- list(
    n_env_outliers = sum(env_outliers, na.rm = TRUE),
    env_outlier_ids = which(env_outliers)
  )
  
  # Genetic outliers (heterozygosity)
  if (!is.null(holobiont$host_genotypes)) {
    het <- rowMeans(holobiont$host_genotypes == 1, na.rm = TRUE)
    het_mean <- mean(het, na.rm = TRUE)
    het_sd <- sd(het, na.rm = TRUE)
    het_outliers <- abs(het - het_mean) > 3 * het_sd
    
    results$n_het_outliers <- sum(het_outliers, na.rm = TRUE)
    results$het_outlier_ids <- which(het_outliers)
  }
  
  results$message <- if (results$n_env_outliers > 0) {
    sprintf("%d environmental outliers detected", results$n_env_outliers)
  } else {
    "No outliers detected"
  }
  
  return(results)
}

#' Assess environmental gradient quality
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with gradient quality metrics
#' @keywords internal
assess_gradient_quality <- function(holobiont, env_variable) {
  
  env <- holobiont$environment[[env_variable]]
  
  if (is.null(env)) {
    return(list(status = "No environmental data"))
  }
  
  # Range
  env_range <- diff(range(env, na.rm = TRUE))
  
  # Coefficient of variation
  cv <- sd(env, na.rm = TRUE) / mean(env, na.rm = TRUE)
  
  # Number of unique values (resolution)
  n_unique <- length(unique(env))
  resolution <- n_unique / length(env)
  
  # Test for uniformity vs gradient
  # A good gradient should have correlation with position
  position <- seq_along(env)
  gradient_test <- cor.test(position, env, method = "spearman")
  
  list(
    range = env_range,
    cv = cv,
    n_unique_values = n_unique,
    resolution = resolution,
    gradient_strength = abs(gradient_test$estimate),
    gradient_p_value = gradient_test$p.value,
    adequate_range = cv > 0.1,
    adequate_resolution = resolution > 0.5,
    message = if (cv < 0.1) {
      "WARNING: Low environmental variation"
    } else if (resolution < 0.3) {
      "WARNING: Low environmental resolution (few unique values)"
    } else {
      "Environmental gradient quality is adequate"
    }
  )
}

#' Check RDA-specific assumptions
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with RDA-specific checks
#' @keywords internal
check_rda_assumptions <- function(holobiont, env_variable) {
  
  list(
    method = "RDA",
    linearity = "RDA assumes linear relationships",
    normality = "RDA is robust to non-normality",
    multicollinearity = "Check VIF if using multiple environmental variables",
    recommendation = "RDA is suitable for most scenarios"
  )
}

#' Check LFMM-specific assumptions
#'
#' @param holobiont Holobiont data object
#' @param env_variable Environmental variable name
#' @return List with LFMM-specific checks
#' @keywords internal
check_lfmm_assumptions <- function(holobiont, env_variable) {
  
  # Estimate number of latent factors needed
  if (!is.null(holobiont$population)) {
    n_pops <- length(unique(holobiont$population))
    recommended_k <- max(1, n_pops - 1)
  } else {
    recommended_k <- "Unknown - run structure analysis"
  }
  
  list(
    method = "LFMM",
    latent_factors = "LFMM requires specifying number of latent factors (K)",
    recommended_k = recommended_k,
    structure_correction = "LFMM explicitly models population structure",
    recommendation = "Use PCA or STRUCTURE to estimate K first"
  )
}

#' Assess overall suitability for GEA
#'
#' @param assumption_results List of assumption check results
#' @return Overall assessment
#' @keywords internal
assess_overall_suitability <- function(assumption_results) {
  
  flags <- character()
  
  # Check each component
  if (!is.null(assumption_results$ge_independence$warning)) {
    if (assumption_results$ge_independence$warning) {
      flags <- c(flags, "Strong G-E confounding")
    }
  }
  
  if (!is.null(assumption_results$structure_confounding$confounded)) {
    if (assumption_results$structure_confounding$confounded) {
      flags <- c(flags, "Environment confounded with population structure")
    }
  }
  
  if (!is.null(assumption_results$signal_noise$adequate_snr)) {
    if (!assumption_results$signal_noise$adequate_snr) {
      flags <- c(flags, "Low signal-to-noise ratio")
    }
  }
  
  if (!is.null(assumption_results$gradient_quality$adequate_range)) {
    if (!assumption_results$gradient_quality$adequate_range) {
      flags <- c(flags, "Insufficient environmental variation")
    }
  }
  
  list(
    suitable = length(flags) == 0,
    flags = flags,
    recommendation = if (length(flags) == 0) {
      "Data appear suitable for GEA analysis"
    } else {
      paste("Potential issues detected:", paste(flags, collapse = "; "))
    }
  )
}

#' Print method for gea_assumptions
#'
#' @param x A gea_assumptions object
#' @param ... Additional arguments (ignored)
#' @export
print.gea_assumptions <- function(x, ...) {
  
  cat("GEA Assumption Checks\n")
  cat("=====================\n\n")
  
  # G-E independence
  if (!is.null(x$ge_independence$message)) {
    cat("Genotype-Environment Independence:\n")
    cat(" ", x$ge_independence$message, "\n")
    if (!is.null(x$ge_independence$max_correlation)) {
      cat(sprintf("  Max PC-environment correlation: %.3f\n", 
                 x$ge_independence$max_correlation))
    }
    cat("\n")
  }
  
  # Structure confounding
  if (!is.null(x$structure_confounding$message)) {
    cat("Population Structure:\n")
    cat(" ", x$structure_confounding$message, "\n\n")
  }
  
  # Signal-to-noise
  if (!is.null(x$signal_noise$message)) {
    cat("Signal-to-Noise Ratio:\n")
    cat(" ", x$signal_noise$message, "\n")
    if (!is.null(x$signal_noise$signal_to_noise_ratio)) {
      cat(sprintf("  SNR: %.2f\n", x$signal_noise$signal_to_noise_ratio))
    }
    cat("\n")
  }
  
  # Outliers
  if (!is.null(x$outliers$message)) {
    cat("Outliers:\n")
    cat(" ", x$outliers$message, "\n\n")
  }
  
  # Gradient quality
  if (!is.null(x$gradient_quality$message)) {
    cat("Environmental Gradient:\n")
    cat(" ", x$gradient_quality$message, "\n")
    if (!is.null(x$gradient_quality$cv)) {
      cat(sprintf("  Coefficient of variation: %.3f\n", x$gradient_quality$cv))
    }
    cat("\n")
  }
  
  # Overall
  cat("Overall Assessment:\n")
  cat(" ", x$overall_suitable$recommendation, "\n")
  
  if (length(x$overall_suitable$flags) > 0) {
    cat("\nFlags:\n")
    for (flag in x$overall_suitable$flags) {
      cat("  -", flag, "\n")
    }
  }
  
  invisible(x)
}

#' Compare power between different scenarios
#'
#' Compare statistical power across different experimental designs or parameter settings.
#' Useful for deciding between alternative study designs.
#'
#' @param ... Named gea_power_analysis objects to compare
#' @return Comparison plot and summary table
#'
#' @examples
#' \dontrun{
#' power1 <- simulate_gea_power(effect_size = 0.2)
#' power2 <- simulate_gea_power(effect_size = 0.3)
#' power3 <- simulate_gea_power(effect_size = 0.4)
#' 
#' compare_power_scenarios(
#'   "Small effect" = power1,
#'   "Medium effect" = power2,
#'   "Large effect" = power3
#' )
#' }
#'
#' @export
compare_power_scenarios <- function(...) {
  
  scenarios <- list(...)
  
  if (length(scenarios) < 2) {
    stop("Need at least 2 scenarios to compare")
  }
  
  # Check all inputs are gea_power_analysis objects
  valid <- sapply(scenarios, inherits, "gea_power_analysis")
  if (!all(valid)) {
    stop("All inputs must be gea_power_analysis objects")
  }
  
  # Get scenario names
  scenario_names <- names(scenarios)
  if (is.null(scenario_names)) {
    scenario_names <- paste("Scenario", seq_along(scenarios))
  }
  
  # Combine data
  combined_data <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
    df <- scenarios[[i]]
    df$scenario <- scenario_names[i]
    df
  }))
  
  # Plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(combined_data, 
                        ggplot2::aes(x = n_individuals, y = power, 
                                    color = scenario, linetype = scenario)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", 
                         color = "gray40") +
      ggplot2::labs(
        title = "Power Comparison Across Scenarios",
        x = "Sample size (n individuals)",
        y = "Statistical power",
        color = "Scenario",
        linetype = "Scenario"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold")
      )
    
    print(p)
  }
  
  # Summary table
  summary_table <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
    result <- scenarios[[i]]
    required_n <- tryCatch(
      calculate_required_n(result, target_power = 0.8),
      error = function(e) NA
    )
    
    params <- attr(result, "parameters")
    
    data.frame(
      Scenario = scenario_names[i],
      Effect_size = params$effect_size,
      N_causal = params$n_causal,
      Required_n_80pct = required_n,
      Max_power_achieved = max(result$power, na.rm = TRUE)
    )
  }))
  
  cat("\nScenario Comparison Summary:\n")
  print(summary_table, row.names = FALSE)
  
  invisible(list(
    plot = if (exists("p")) p else NULL,
    summary = summary_table,
    data = combined_data
  ))
}