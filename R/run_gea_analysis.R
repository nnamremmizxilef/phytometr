#' Run GEA analysis on holobiont data
#'
#' @param holobiont Holobiont data object with genotypes, environment, and response
#' @param params Analysis parameters including gea_method, alpha, correction_method
#' @return List with power, fdr, n_detected, p_values, and detected_loci indices
#' @keywords internal
run_gea_analysis <- function(holobiont, params) {
  
  # Extract data
  genotypes <- holobiont$host_genotypes
  
  # Get environment variable
  if (params$env_type == "continuous") {
    environment <- holobiont$environment$env1
  } else {
    environment <- holobiont$environment[[1]]
  }
  
  # Get response variable
  if (params$response_variable == "abundance") {
    response <- holobiont$microbiome$response
  } else if (params$response_variable == "diversity") {
    response <- holobiont$microbiome$diversity
  } else if (params$response_variable == "composition") {
    response <- holobiont$microbiome$pc1
  } else if (params$response_variable == "phenotype") {
    response <- holobiont$phenotype$trait1
  }
  
  # Check dimensions
  n_ind <- nrow(genotypes)
  n_snps <- ncol(genotypes)
  
  if (length(environment) != n_ind) {
    stop("Environment and genotype dimensions don't match")
  }
  if (length(response) != n_ind) {
    stop("Response and genotype dimensions don't match")
  }
  
  # Run GEA method
  if (params$gea_method == "rda") {
    result <- run_rda_gea(genotypes, environment, response, params)
  } else if (params$gea_method == "lfmm") {
    result <- run_lfmm_gea(genotypes, environment, response, params)
  } else {
    # Default to RDA
    result <- run_rda_gea(genotypes, environment, response, params)
  }
  
  return(result)
}

#' Run RDA-based GEA
#' @keywords internal
run_rda_gea <- function(genotypes, environment, response, params) {
  
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' required for RDA. Install with: install.packages('vegan')")
  }
  
  # Prepare data
  # RDA: response ~ environment, conditioned on genotypes
  # We want to test genotype-environment associations affecting response
  
  # Create SNP data frame
  snp_df <- as.data.frame(genotypes)
  colnames(snp_df) <- paste0("SNP", seq_len(ncol(genotypes)))
  
  # Add environment and response
  snp_df$environment <- scale(environment)
  snp_df$response <- scale(response)
  
  # Test each SNP for association with response, accounting for environment
  p_values <- numeric(ncol(genotypes))
  
  for (i in seq_len(ncol(genotypes))) {
    
    # Model: response ~ SNP + environment
    tryCatch({
      model <- lm(response ~ genotypes[, i] + environment)
      # P-value for SNP coefficient
      p_values[i] <- summary(model)$coefficients[2, 4]
    }, error = function(e) {
      p_values[i] <- 1  # Assign non-significant if error
    })
  }
  
  # Apply multiple testing correction
  if (params$correction_method == "bonferroni") {
    p_adjusted <- p.adjust(p_values, method = "bonferroni")
  } else if (params$correction_method == "fdr") {
    p_adjusted <- p.adjust(p_values, method = "fdr")
  } else {
    p_adjusted <- p_values
  }
  
  # Detect significant SNPs
  detected <- which(p_adjusted < params$alpha)
  
  # Calculate power (if we know causal loci)
  if (!is.null(holobiont$causal_loci)) {
    true_positives <- sum(detected %in% holobiont$causal_loci)
    power <- true_positives / length(holobiont$causal_loci)
    
    false_positives <- sum(!(detected %in% holobiont$causal_loci))
    fdr <- if (length(detected) > 0) {
      false_positives / length(detected)
    } else {
      0
    }
  } else {
    power <- NA
    fdr <- NA
  }
  
  result <- list(
    p_values = p_values,
    p_adjusted = p_adjusted,
    detected_loci = detected,
    n_detected = length(detected),
    power = power,
    fdr = fdr
  )
  
  return(result)
}

#' Run LFMM-based GEA
#' @keywords internal
run_lfmm_gea <- function(genotypes, environment, response, params) {
  
  # Check if LEA package available
  if (!requireNamespace("LEA", quietly = TRUE)) {
    warning("Package 'LEA' not available. Falling back to RDA.")
    return(run_rda_gea(genotypes, environment, response, params))
  }
  
  # For now, fall back to RDA
  # TODO: Implement LFMM properly
  warning("LFMM not yet fully implemented. Using RDA.")
  return(run_rda_gea(genotypes, environment, response, params))
}