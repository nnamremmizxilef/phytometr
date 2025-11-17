#' Diagnostic functions for holobiont data and GEA analyses
#'
#' @name diagnostics
#' @keywords internal
NULL

#' Diagnose holobiont data quality
#'
#' Comprehensive quality check for holobiont datasets including genetic data,
#' environmental variables, and holobiont responses.
#'
#' @param holobiont A holobiont_data object
#' @param verbose Logical. Print detailed diagnostics? (default: TRUE)
#' @return List of diagnostic results with class 'holobiont_diagnostics'
#'
#' @examples
#' \dontrun{
#' holo <- create_holobiont_data(100)
#' holo <- simulate_host_genetics(holo, params)
#' diag <- diagnose_holobiont_data(holo)
#' print(diag)
#' plot(diag)
#' }
#'
#' @export
diagnose_holobiont_data <- function(holobiont, verbose = TRUE) {
  
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  
  diagnostics <- list()
  
  # Sample size
  diagnostics$n_individuals <- length(holobiont$host_id)
  
  # Genetic data diagnostics
  if (!is.null(holobiont$host_genotypes)) {
    diagnostics$genetic <- diagnose_genetic_data(holobiont$host_genotypes)
  } else {
    diagnostics$genetic <- list(status = "No genetic data")
  }
  
  # Environmental data diagnostics
  diagnostics$environment <- diagnose_environmental_data(holobiont$environment)
  
  # Holobiont response diagnostics
  diagnostics$holobiont_response <- diagnose_holobiont_response(holobiont)
  
  # Population structure diagnostics
  if (!is.null(holobiont$population)) {
    diagnostics$population <- diagnose_population_structure(holobiont)
  }
  
  # Overall data completeness
  diagnostics$completeness <- calculate_data_completeness(holobiont)
  
  # Sample size recommendations
  diagnostics$sample_size_assessment <- assess_sample_size(holobiont)
  
  class(diagnostics) <- "holobiont_diagnostics"
  
  if (verbose) {
    print(diagnostics)
  }
  
  invisible(diagnostics)
}

#' Diagnose genetic data quality
#'
#' @param snps SNP matrix
#' @return List of genetic diagnostics
#' @keywords internal
diagnose_genetic_data <- function(snps) {
  
  diag <- list()
  
  diag$n_loci <- ncol(snps)
  diag$n_individuals <- nrow(snps)
  
  # Missing data
  diag$missing_rate <- sum(is.na(snps)) / length(snps)
  diag$missing_per_locus <- colMeans(is.na(snps))
  diag$missing_per_individual <- rowMeans(is.na(snps))
  
  # Allele frequencies
  afs <- colMeans(snps, na.rm = TRUE) / 2
  diag$maf <- pmin(afs, 1 - afs)
  diag$mean_maf <- mean(diag$maf, na.rm = TRUE)
  diag$n_rare_variants <- sum(diag$maf < 0.05, na.rm = TRUE)
  
  # Heterozygosity
  het <- rowMeans(snps == 1, na.rm = TRUE)
  diag$heterozygosity <- het
  diag$mean_heterozygosity <- mean(het, na.rm = TRUE)
  diag$sd_heterozygosity <- sd(het, na.rm = TRUE)
  
  # Linkage disequilibrium (sample)
  if (ncol(snps) > 1 && ncol(snps) <= 1000) {
    diag$ld_estimate <- estimate_ld_decay(snps)
  } else if (ncol(snps) > 1000) {
    # Sample 1000 SNPs for LD calculation
    sample_snps <- sample(ncol(snps), 1000)
    diag$ld_estimate <- estimate_ld_decay(snps[, sample_snps])
  }
  
  # Genotype quality flags
  diag$flags <- list()
  if (diag$missing_rate > 0.2) {
    diag$flags$high_missing <- "WARNING: >20% missing data"
  }
  if (diag$n_rare_variants / diag$n_loci > 0.5) {
    diag$flags$many_rare <- "WARNING: >50% variants with MAF < 0.05"
  }
  if (any(diag$heterozygosity < 0.01 | diag$heterozygosity > 0.8)) {
    diag$flags$extreme_het <- "WARNING: Individuals with extreme heterozygosity detected"
  }
  
  return(diag)
}

#' Estimate linkage disequilibrium decay
#'
#' @param snps SNP matrix
#' @param max_pairs Maximum number of SNP pairs to sample (default: 10000)
#' @return List with LD statistics
#' @keywords internal
estimate_ld_decay <- function(snps, max_pairs = 10000) {
  
  n_snps <- ncol(snps)
  
  # Sample SNP pairs if too many
  if (n_snps * (n_snps - 1) / 2 > max_pairs) {
    pairs <- sample_snp_pairs(n_snps, max_pairs)
  } else {
    pairs <- which(upper.tri(matrix(0, n_snps, n_snps)), arr.ind = TRUE)
  }
  
  # Calculate r^2 for sampled pairs
  r2_values <- apply(pairs, 1, function(pair) {
    calculate_r2(snps[, pair[1]], snps[, pair[2]])
  })
  
  # Calculate distance (assuming uniform spacing for simulated data)
  distances <- abs(pairs[, 1] - pairs[, 2])
  
  # Summary statistics
  list(
    mean_r2 = mean(r2_values, na.rm = TRUE),
    median_r2 = median(r2_values, na.rm = TRUE),
    max_r2 = max(r2_values, na.rm = TRUE),
    n_pairs_high_ld = sum(r2_values > 0.8, na.rm = TRUE),
    ld_decay = estimate_decay_rate(distances, r2_values)
  )
}

#' Calculate r^2 between two SNPs
#'
#' @param snp1 First SNP vector
#' @param snp2 Second SNP vector
#' @return r^2 value
#' @keywords internal
calculate_r2 <- function(snp1, snp2) {
  valid <- !is.na(snp1) & !is.na(snp2)
  if (sum(valid) < 10) return(NA)
  
  cor_val <- cor(snp1[valid], snp2[valid])
  return(cor_val^2)
}

#' Sample SNP pairs for LD calculation
#'
#' @param n_snps Number of SNPs
#' @param n_pairs Number of pairs to sample
#' @return Matrix of SNP pair indices
#' @keywords internal
sample_snp_pairs <- function(n_snps, n_pairs) {
  idx1 <- sample(n_snps, n_pairs, replace = TRUE)
  idx2 <- sample(n_snps, n_pairs, replace = TRUE)
  
  # Ensure different SNPs
  same <- idx1 == idx2
  while (any(same)) {
    idx2[same] <- sample(n_snps, sum(same), replace = TRUE)
    same <- idx1 == idx2
  }
  
  cbind(idx1, idx2)
}

#' Estimate LD decay rate
#'
#' @param distances Vector of distances between SNP pairs
#' @param r2_values Vector of r^2 values
#' @return Decay rate estimate
#' @keywords internal
estimate_decay_rate <- function(distances, r2_values) {
  valid <- !is.na(r2_values) & distances > 0
  
  if (sum(valid) < 10) return(NA)
  
  # Fit exponential decay: r^2 ~ exp(-distance * decay_rate)
  tryCatch({
    fit <- nls(r2_values[valid] ~ a * exp(-b * distances[valid]),
              start = list(a = 0.5, b = 0.01),
              control = nls.control(maxiter = 100, warnOnly = TRUE))
    coef(fit)["b"]
  }, error = function(e) {
    NA
  })
}

#' Diagnose environmental data
#'
#' @param environment List of environmental variables
#' @return List of environmental diagnostics
#' @keywords internal
diagnose_environmental_data <- function(environment) {
  
  diag <- list()
  
  # Count available variables
  available <- !sapply(environment, is.null)
  diag$n_variables <- sum(available)
  diag$available_variables <- names(environment)[available]
  
  if (diag$n_variables == 0) {
    diag$status <- "No environmental data"
    return(diag)
  }
  
  # Analyze each variable
  diag$variables <- lapply(environment[available], function(var) {
    
    var_diag <- list()
    
    if (is.numeric(var)) {
      var_diag$type <- "continuous"
      var_diag$mean <- mean(var, na.rm = TRUE)
      var_diag$sd <- sd(var, na.rm = TRUE)
      var_diag$range <- range(var, na.rm = TRUE)
      var_diag$missing <- sum(is.na(var))
      var_diag$cv <- var_diag$sd / var_diag$mean  # Coefficient of variation
      
      # Check for gradient
      var_diag$has_gradient <- assess_gradient_strength(var)
      
      # Check for spatial autocorrelation (if enough variation)
      if (var_diag$sd > 0) {
        var_diag$spatial_structure <- estimate_spatial_autocorrelation(var)
      }
      
    } else {
      var_diag$type <- "categorical"
      var_diag$categories <- unique(var)
      var_diag$n_categories <- length(var_diag$categories)
      var_diag$frequencies <- table(var)
    }
    
    var_diag
  })
  
  # Check for collinearity between variables
  if (diag$n_variables > 1) {
    diag$collinearity <- assess_env_collinearity(environment[available])
  }
  
  return(diag)
}

#' Assess strength of environmental gradient
#'
#' @param env_var Environmental variable vector
#' @return List with gradient assessment
#' @keywords internal
assess_gradient_strength <- function(env_var) {
  
  valid <- !is.na(env_var)
  if (sum(valid) < 10) return(list(present = FALSE))
  
  # Test for monotonic trend (Spearman with position)
  position <- seq_along(env_var)[valid]
  cor_test <- cor.test(position, env_var[valid], method = "spearman")
  
  # Test for spatial autocorrelation
  autocor <- calculate_autocorrelation(env_var[valid])
  
  list(
    present = abs(cor_test$estimate) > 0.3 || autocor > 0.3,
    strength = abs(cor_test$estimate),
    autocorrelation = autocor,
    p_value = cor_test$p.value
  )
}

#' Calculate autocorrelation in environmental data
#'
#' @param x Numeric vector
#' @return Lag-1 autocorrelation
#' @keywords internal
calculate_autocorrelation <- function(x) {
  if (length(x) < 3) return(0)
  
  # Lag-1 autocorrelation
  cor(x[-length(x)], x[-1], use = "complete.obs")
}

#' Estimate spatial autocorrelation
#'
#' @param env_var Environmental variable
#' @return Spatial autocorrelation estimate
#' @keywords internal
estimate_spatial_autocorrelation <- function(env_var) {
  
  # Moran's I approximation using sequential ordering
  n <- length(env_var)
  valid <- !is.na(env_var)
  
  if (sum(valid) < 10) return(NA)
  
  # Calculate Moran's I with nearest neighbor weights
  x <- env_var[valid]
  x_mean <- mean(x)
  
  numerator <- 0
  denominator <- sum((x - x_mean)^2)
  
  for (i in seq_along(x)[-length(x)]) {
    numerator <- numerator + (x[i] - x_mean) * (x[i+1] - x_mean)
  }
  
  morans_i <- (length(x) - 1) * numerator / denominator
  
  return(morans_i)
}

#' Assess environmental collinearity
#'
#' @param env_list List of environmental variables
#' @return Matrix of correlations and VIF values
#' @keywords internal
assess_env_collinearity <- function(env_list) {
  
  # Extract numeric variables
  numeric_vars <- env_list[sapply(env_list, is.numeric)]
  
  if (length(numeric_vars) < 2) {
    return(list(status = "Need at least 2 numeric variables"))
  }
  
  # Create data frame
  env_df <- as.data.frame(numeric_vars)
  
  # Correlation matrix
  cor_matrix <- cor(env_df, use = "pairwise.complete.obs")
  
  # Calculate VIF (Variance Inflation Factor)
  vif_values <- sapply(names(env_df), function(var) {
    if (ncol(env_df) == 2) {
      # Simple case: VIF = 1/(1-r^2)
      other_var <- setdiff(names(env_df), var)
      r2 <- cor_matrix[var, other_var]^2
      return(1 / (1 - r2))
    } else {
      # Multiple regression
      formula_str <- paste(var, "~", paste(setdiff(names(env_df), var), collapse = " + "))
      tryCatch({
        model <- lm(as.formula(formula_str), data = env_df)
        r2 <- summary(model)$r.squared
        1 / (1 - r2)
      }, error = function(e) NA)
    }
  })
  
  list(
    correlation_matrix = cor_matrix,
    vif = vif_values,
    high_collinearity = any(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.7, na.rm = TRUE),
    high_vif = any(vif_values > 10, na.rm = TRUE)
  )
}

#' Diagnose holobiont response data
#'
#' @param holobiont Holobiont data object
#' @return List of response diagnostics
#' @keywords internal
diagnose_holobiont_response <- function(holobiont) {
  
  diag <- list()
  
  # Check microbiome data
  if (!is.null(holobiont$microbiome$target_taxon)) {
    diag$abundance <- list(
      mean = mean(holobiont$microbiome$target_taxon, na.rm = TRUE),
      sd = sd(holobiont$microbiome$target_taxon, na.rm = TRUE),
      range = range(holobiont$microbiome$target_taxon, na.rm = TRUE),
      n_zeros = sum(holobiont$microbiome$target_taxon == 0, na.rm = TRUE),
      type = "taxon abundance"
    )
  }
  
  if (!is.null(holobiont$microbiome$alpha_diversity)) {
    diag$diversity <- list(
      mean = mean(holobiont$microbiome$alpha_diversity, na.rm = TRUE),
      sd = sd(holobiont$microbiome$alpha_diversity, na.rm = TRUE),
      range = range(holobiont$microbiome$alpha_diversity, na.rm = TRUE),
      type = "alpha diversity"
    )
  }
  
  if (!is.null(holobiont$microbiome$pc1)) {
    diag$composition <- list(
      mean = mean(holobiont$microbiome$pc1, na.rm = TRUE),
      sd = sd(holobiont$microbiome$pc1, na.rm = TRUE),
      range = range(holobiont$microbiome$pc1, na.rm = TRUE),
      type = "community composition PC1"
    )
  }
  
  if (!is.null(holobiont$microbiome$fungal_community)) {
    diag$community <- diagnose_community_matrix(holobiont$microbiome$fungal_community)
  }
  
  # Check phenotype data
  pheno_available <- !sapply(holobiont$phenotypes, is.null)
  if (any(pheno_available)) {
    diag$phenotypes <- names(holobiont$phenotypes)[pheno_available]
  }
  
  if (length(diag) == 0) {
    diag$status <- "No response data"
  }
  
  return(diag)
}

#' Diagnose community matrix
#'
#' @param community_matrix Community matrix (samples x taxa)
#' @return List of community diagnostics
#' @keywords internal
diagnose_community_matrix <- function(community_matrix) {
  
  diag <- list()
  
  diag$n_samples <- nrow(community_matrix)
  diag$n_taxa <- ncol(community_matrix)
  
  # Prevalence
  diag$taxa_prevalence <- colMeans(community_matrix > 0)
  diag$mean_prevalence <- mean(diag$taxa_prevalence)
  diag$n_rare_taxa <- sum(diag$taxa_prevalence < 0.1)
  
  # Abundance distribution
  diag$mean_abundance_per_taxon <- colMeans(community_matrix)
  diag$evenness <- calculate_evenness(community_matrix)
  
  # Check for compositional data
  row_sums <- rowSums(community_matrix)
  diag$is_compositional <- all(abs(row_sums - 1) < 1e-6)
  
  return(diag)
}

#' Calculate community evenness
#'
#' @param community_matrix Community matrix
#' @return Vector of evenness values (Pielou's J)
#' @keywords internal
calculate_evenness <- function(community_matrix) {
  
  apply(community_matrix, 1, function(row) {
    present <- row > 0
    if (sum(present) <= 1) return(NA)
    
    p <- row[present] / sum(row[present])
    H <- -sum(p * log(p))  # Shannon diversity
    S <- sum(present)      # Richness
    J <- H / log(S)        # Pielou's evenness
    
    return(J)
  })
}

#' Diagnose population structure
#'
#' @param holobiont Holobiont data object
#' @return List of population structure diagnostics
#' @keywords internal
diagnose_population_structure <- function(holobiont) {
  
  diag <- list()
  
  if (is.null(holobiont$population)) {
    return(list(status = "No population structure"))
  }
  
  diag$n_populations <- length(unique(holobiont$population))
  diag$population_sizes <- table(holobiont$population)
  diag$balanced <- sd(as.numeric(diag$population_sizes)) / 
                   mean(as.numeric(diag$population_sizes)) < 0.2
  
  # Calculate Fst if genetic data available
  if (!is.null(holobiont$host_genotypes)) {
    diag$fst <- calculate_fst(holobiont$host_genotypes, holobiont$population)
  }
  
  return(diag)
}

#' Calculate Fst between populations
#'
#' @param snps SNP matrix
#' @param populations Population assignments
#' @return Weighted average Fst
#' @keywords internal
calculate_fst <- function(snps, populations) {
  
  pops <- unique(populations)
  
  if (length(pops) < 2) return(NA)
  
  # Calculate allele frequencies per population
  pop_freqs <- lapply(pops, function(pop) {
    pop_snps <- snps[populations == pop, , drop = FALSE]
    colMeans(pop_snps, na.rm = TRUE) / 2
  })
  
  # Overall allele frequencies
  overall_freq <- colMeans(snps, na.rm = TRUE) / 2
  
  # Calculate Fst per locus (Weir & Cockerham 1984 approximation)
  fst_per_locus <- sapply(seq_along(overall_freq), function(i) {
    p <- overall_freq[i]
    if (is.na(p) || p == 0 || p == 1) return(NA)
    
    # Variance in allele frequencies among populations
    pop_p <- sapply(pop_freqs, function(pf) pf[i])
    var_p <- var(pop_p, na.rm = TRUE)
    
    # Expected heterozygosity
    He <- 2 * p * (1 - p)
    
    # Fst
    if (He > 0) {
      var_p / (p * (1 - p))
    } else {
      NA
    }
  })
  
  # Weighted mean Fst
  mean(fst_per_locus, na.rm = TRUE)
}

#' Calculate data completeness
#'
#' @param holobiont Holobiont data object
#' @return Data completeness score (0-1)
#' @keywords internal
calculate_data_completeness <- function(holobiont) {
  
  scores <- numeric()
  
  # Genetic data
  scores["genetic"] <- ifelse(is.null(holobiont$host_genotypes), 0, 1)
  
  # Environmental data
  n_env <- sum(!sapply(holobiont$environment, is.null))
  scores["environment"] <- min(n_env / 3, 1)  # Expect at least 3 variables
  
  # Holobiont response
  n_response <- sum(!sapply(holobiont$microbiome, is.null)) +
                sum(!sapply(holobiont$phenotypes, is.null))
  scores["response"] <- min(n_response / 2, 1)  # Expect at least 2 types
  
  list(
    overall = mean(scores),
    components = scores
  )
}

#' Assess sample size adequacy
#'
#' @param holobiont Holobiont data object
#' @return List with sample size assessment
#' @keywords internal
assess_sample_size <- function(holobiont) {
  
  n <- length(holobiont$host_id)
  
  assessment <- list(
    current_n = n,
    adequate_for_basic_gea = n >= 100,
    adequate_for_gxe = n >= 200,
    adequate_for_rare_variants = n >= 500
  )
  
  # Context-specific recommendations
  if (!is.null(holobiont$host_genotypes)) {
    n_loci <- ncol(holobiont$host_genotypes)
    assessment$snp_to_sample_ratio <- n_loci / n
    assessment$multiple_testing_burden <- n_loci > n * 10
  }
  
  if (!is.null(holobiont$population)) {
    n_pops <- length(unique(holobiont$population))
    assessment$samples_per_population <- n / n_pops
    assessment$adequate_per_population <- assessment$samples_per_population >= 20
  }
  
  return(assessment)
}

#' Print method for holobiont_diagnostics
#'
#' @param x A holobiont_diagnostics object
#' @param ... Additional arguments (ignored)
#' @export
print.holobiont_diagnostics <- function(x, ...) {
  
  cat("Holobiont Data Diagnostics\n")
  cat("==========================\n\n")
  
  cat("Sample size:", x$n_individuals, "individuals\n\n")
  
  # Genetic data
  if (!is.null(x$genetic$n_loci)) {
    cat("Genetic Data:\n")
    cat(sprintf("  Number of loci: %d\n", x$genetic$n_loci))
    cat(sprintf("  Missing data rate: %.2f%%\n", x$genetic$missing_rate * 100))
    cat(sprintf("  Mean MAF: %.3f\n", x$genetic$mean_maf))
    cat(sprintf("  Rare variants (MAF < 0.05): %d (%.1f%%)\n", 
               x$genetic$n_rare_variants,
               100 * x$genetic$n_rare_variants / x$genetic$n_loci))
    cat(sprintf("  Mean heterozygosity: %.3f (SD: %.3f)\n", 
               x$genetic$mean_heterozygosity,
               x$genetic$sd_heterozygosity))
    
    if (!is.null(x$genetic$ld_estimate)) {
      cat(sprintf("  Mean LD (r^2): %.3f\n", x$genetic$ld_estimate$mean_r2))
      cat(sprintf("  SNP pairs with high LD (r^2 > 0.8): %d\n", 
                 x$genetic$ld_estimate$n_pairs_high_ld))
    }
    
    if (length(x$genetic$flags) > 0) {
      cat("\n  Quality flags:\n")
      for (flag in x$genetic$flags) {
        cat("   ", flag, "\n")
      }
    }
    cat("\n")
  }
  
  # Environmental data
  if (x$environment$n_variables > 0) {
    cat("Environmental Data:\n")
    cat(sprintf("  Number of variables: %d\n", x$environment$n_variables))
    cat("  Variables:", paste(x$environment$available_variables, collapse = ", "), "\n")
    
    for (var_name in names(x$environment$variables)) {
      var <- x$environment$variables[[var_name]]
      cat(sprintf("\n  %s (%s):\n", var_name, var$type))
      
      if (var$type == "continuous") {
        cat(sprintf("    Mean: %.2f (SD: %.2f)\n", var$mean, var$sd))
        cat(sprintf("    Range: [%.2f, %.2f]\n", var$range[1], var$range[2]))
        cat(sprintf("    CV: %.2f\n", var$cv))
        
        if (!is.null(var$has_gradient$present)) {
          if (var$has_gradient$present) {
            cat(sprintf("    Gradient detected: strength = %.2f, p = %.4f\n", 
                       var$has_gradient$strength, var$has_gradient$p_value))
          }
        }
      } else {
        cat(sprintf("    Categories: %d\n", var$n_categories))
      }
    }
    
    if (!is.null(x$environment$collinearity)) {
      cat("\n  Collinearity:\n")
      if (x$environment$collinearity$high_collinearity) {
        cat("    WARNING: High correlation detected between variables\n")
      }
      if (x$environment$collinearity$high_vif) {
        cat("    WARNING: High VIF detected (>10)\n")
      }
      cat("    VIF values:\n")
      for (i in seq_along(x$environment$collinearity$vif)) {
        cat(sprintf("      %s: %.2f\n", 
                   names(x$environment$collinearity$vif)[i],
                   x$environment$collinearity$vif[i]))
      }
    }
    cat("\n")
  }
  
  # Holobiont response
  if (length(x$holobiont_response) > 1 || 
      (length(x$holobiont_response) == 1 && is.null(x$holobiont_response$status))) {
    cat("Holobiont Response Data:\n")
    
    for (response_type in names(x$holobiont_response)) {
      if (response_type == "status") next
      
      resp <- x$holobiont_response[[response_type]]
      if (!is.null(resp$type)) {
        cat(sprintf("  %s:\n", resp$type))
        cat(sprintf("    Mean: %.2f (SD: %.2f)\n", resp$mean, resp$sd))
        cat(sprintf("    Range: [%.2f, %.2f]\n", resp$range[1], resp$range[2]))
        
        if (!is.null(resp$n_zeros)) {
          cat(sprintf("    Zero inflation: %d (%.1f%%)\n", 
                     resp$n_zeros, 100 * resp$n_zeros / x$n_individuals))
        }
      }
    }
    
    if (!is.null(x$holobiont_response$community)) {
      cat("\n  Community matrix:\n")
      cat(sprintf("    Taxa: %d\n", x$holobiont_response$community$n_taxa))
      cat(sprintf("    Mean prevalence: %.2f\n", 
                 x$holobiont_response$community$mean_prevalence))
      cat(sprintf("    Rare taxa (prev < 0.1): %d\n", 
                 x$holobiont_response$community$n_rare_taxa))
      cat(sprintf("    Compositional: %s\n", 
                 x$holobiont_response$community$is_compositional))
    }
    cat("\n")
  }
  
  # Population structure
  if (!is.null(x$population)) {
    cat("Population Structure:\n")
    cat(sprintf("  Number of populations: %d\n", x$population$n_populations))
    cat(sprintf("  Population sizes: %s\n", 
               paste(x$population$population_sizes, collapse = ", ")))
    cat(sprintf("  Balanced design: %s\n", x$population$balanced))
    if (!is.null(x$population$fst)) {
      cat(sprintf("  Mean Fst: %.4f\n", x$population$fst))
    }
    cat("\n")
  }
  
  # Data completeness
  cat("Data Completeness:\n")
  cat(sprintf("  Overall: %.0f%%\n", x$completeness$overall * 100))
  cat("  Components:\n")
  for (comp in names(x$completeness$components)) {
    cat(sprintf("    %s: %.0f%%\n", comp, x$completeness$components[comp] * 100))
  }
  cat("\n")
  
  # Sample size assessment
  cat("Sample Size Assessment:\n")
  cat(sprintf("  Current n: %d\n", x$sample_size_assessment$current_n))
  cat(sprintf("  Adequate for basic GEA: %s\n", 
             x$sample_size_assessment$adequate_for_basic_gea))
  cat(sprintf("  Adequate for GxE: %s\n", 
             x$sample_size_assessment$adequate_for_gxe))
  
  if (!is.null(x$sample_size_assessment$multiple_testing_burden)) {
    if (x$sample_size_assessment$multiple_testing_burden) {
      cat("  WARNING: High multiple testing burden (loci >> samples)\n")
    }
  }
  
  if (!is.null(x$sample_size_assessment$adequate_per_population)) {
    if (!x$sample_size_assessment$adequate_per_population) {
      cat(sprintf("  WARNING: Low samples per population (%.1f)\n",
                 x$sample_size_assessment$samples_per_population))
    }
  }
  
  invisible(x)
}

#' Plot holobiont diagnostics
#'
#' @param x A holobiont_diagnostics object
#' @param which Character vector. Which plots to generate: "genetic", "environment", 
#'   "response", "all" (default: "all")
#' @param ... Additional arguments passed to plotting functions
#' @return List of ggplot objects (invisible)
#' @export
plot.holobiont_diagnostics <- function(x, which = "all", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting diagnostics")
  }
  
  plots <- list()
  
  # Determine which plots to make
  if ("all" %in% which) {
    which <- c("genetic", "environment", "response")
  }
  
  # Genetic data plots
  if ("genetic" %in% which && !is.null(x$genetic$n_loci)) {
    plots$maf_distribution <- plot_maf_distribution(x$genetic)
    plots$heterozygosity <- plot_heterozygosity(x$genetic)
    plots$missing_data <- plot_missing_data(x$genetic)
  }
  
  # Environmental data plots
  if ("environment" %in% which && x$environment$n_variables > 0) {
    plots$environment <- plot_environmental_diagnostics(x$environment)
  }
  
  # Response data plots
  if ("response" %in% which && length(x$holobiont_response) > 1) {
    plots$response <- plot_response_diagnostics(x$holobiont_response)
  }
  
  # Display plots
  for (plot_name in names(plots)) {
    print(plots[[plot_name]])
    if (plot_name != names(plots)[length(plots)]) {
      readline(prompt = "Press [enter] for next plot")
    }
  }
  
  invisible(plots)
}

#' Plot MAF distribution
#'
#' @param genetic_diag Genetic diagnostics list
#' @return ggplot object
#' @keywords internal
plot_maf_distribution <- function(genetic_diag) {
  
  df <- data.frame(maf = genetic_diag$maf)
  
  ggplot2::ggplot(df, ggplot2::aes(x = maf)) +
    ggplot2::geom_histogram(bins = 50, fill = "#2c7fb8", color = "white") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", 
                       color = "#d7301f", linewidth = 1) +
    ggplot2::annotate("text", x = 0.05, y = Inf, 
                     label = "MAF = 0.05", hjust = -0.1, vjust = 2) +
    ggplot2::labs(
      title = "Minor Allele Frequency Distribution",
      x = "Minor Allele Frequency",
      y = "Count"
    ) +
    ggplot2::theme_minimal()
}

#' Plot heterozygosity
#'
#' @param genetic_diag Genetic diagnostics list
#' @return ggplot object
#' @keywords internal
plot_heterozygosity <- function(genetic_diag) {
  
  df <- data.frame(
    individual = seq_along(genetic_diag$heterozygosity),
    heterozygosity = genetic_diag$heterozygosity
  )
  
  mean_het <- mean(genetic_diag$heterozygosity, na.rm = TRUE)
  sd_het <- sd(genetic_diag$heterozygosity, na.rm = TRUE)
  
  ggplot2::ggplot(df, ggplot2::aes(x = individual, y = heterozygosity)) +
    ggplot2::geom_point(alpha = 0.5, color = "#2c7fb8") +
    ggplot2::geom_hline(yintercept = mean_het, linetype = "solid", color = "black") +
    ggplot2::geom_hline(yintercept = mean_het + 3*sd_het, 
                       linetype = "dashed", color = "#d7301f") +
    ggplot2::geom_hline(yintercept = mean_het - 3*sd_het, 
                       linetype = "dashed", color = "#d7301f") +
    ggplot2::labs(
      title = "Individual Heterozygosity",
      subtitle = sprintf("Mean: %.3f (SD: %.3f)", mean_het, sd_het),
      x = "Individual",
      y = "Heterozygosity"
    ) +
    ggplot2::theme_minimal()
}

#' Plot missing data
#'
#' @param genetic_diag Genetic diagnostics list
#' @return ggplot object
#' @keywords internal
plot_missing_data <- function(genetic_diag) {
  
  df <- data.frame(
    individual = seq_along(genetic_diag$missing_per_individual),
    missing_rate = genetic_diag$missing_per_individual * 100
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = individual, y = missing_rate)) +
    ggplot2::geom_col(fill = "#2c7fb8") +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", 
                       color = "#d7301f", linewidth = 1) +
    ggplot2::labs(
      title = "Missing Data per Individual",
      x = "Individual",
      y = "Missing data (%)"
    ) +
    ggplot2::theme_minimal()
}

#' Plot environmental diagnostics
#'
#' @param env_diag Environmental diagnostics list
#' @return ggplot object
#' @keywords internal
plot_environmental_diagnostics <- function(env_diag) {
  
  # Extract continuous variables for correlation plot
  if (!is.null(env_diag$collinearity)) {
    cor_matrix <- env_diag$collinearity$correlation_matrix
    
    # Reshape for plotting
    cor_df <- reshape_cor_matrix(cor_matrix)
    
    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                   midpoint = 0, limits = c(-1, 1),
                                   name = "Correlation") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 3) +
      ggplot2::labs(
        title = "Environmental Variable Correlations",
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  return(NULL)
}

#' Reshape correlation matrix for plotting
#'
#' @param cor_matrix Correlation matrix
#' @return Data frame for ggplot
#' @keywords internal
reshape_cor_matrix <- function(cor_matrix) {
  cor_df <- as.data.frame(as.table(cor_matrix))
  names(cor_df) <- c("Var1", "Var2", "value")
  cor_df
}

#' Plot response diagnostics
#'
#' @param response_diag Response diagnostics list
#' @return ggplot object or NULL
#' @keywords internal
plot_response_diagnostics <- function(response_diag) {
  # Placeholder for response-specific plots
  # Could add distribution plots, Q-Q plots, etc.
  return(NULL)
}