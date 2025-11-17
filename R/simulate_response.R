#' Simulate holobiont response variable
#'
#' @param holobiont A holobiont_data object
#' @param params List of simulation parameters
#' @return Updated holobiont_data object with response data and causal_loci
#' @keywords internal
simulate_response <- function(holobiont, params) {
  
  snps <- holobiont$host_genotypes
  env <- holobiont$environment$main
  n <- nrow(snps)
  n_snps <- ncol(snps)
  
  # Check if we have enough SNPs
  if (n_snps < params$n_causal) {
    stop(sprintf("Not enough SNPs (%d) for requested causal loci (%d). Increase n_loci or decrease maf_threshold.",
                n_snps, params$n_causal))
  }
  
  # Select causal loci
  causal_loci <- sample(n_snps, params$n_causal)
  holobiont$causal_loci <- causal_loci
  
  # Generate effect sizes
  effects <- generate_effect_sizes(params)
  
  # Calculate genetic component
  genetic_value <- calculate_genetic_value(snps, causal_loci, effects, params)
  
  # Calculate environmental component
  env_value <- calculate_environmental_effect(env, params)
  
  # Calculate GxE interaction if specified
  gxe_value <- if (params$gxe_interaction) {
    genetic_value * env_value * params$gxe_strength
  } else {
    0
  }
  
  # Calculate residual variance
  total_explained <- params$holobiont_heritability + params$env_effect_on_holobiont
  if (params$gxe_interaction) {
    total_explained <- total_explained + (params$gxe_strength^2)
  }
  
  residual_var <- max(0, 1 - total_explained)
  noise <- rnorm(n, sd = sqrt(residual_var))
  
  # Combined response
  response <- genetic_value + env_value + gxe_value + noise
  
  # Format based on response type
  holobiont <- format_response_variable(holobiont, response, params)
  
  return(holobiont)
}

#' Generate effect sizes for causal loci
#'
#' @param params Simulation parameters
#' @return Numeric vector of effect sizes
#' @keywords internal
generate_effect_sizes <- function(params) {
  
  if (params$effect_size_distribution == "equal") {
    effects <- rep(params$effect_size, params$n_causal)
    
  } else if (params$effect_size_distribution == "exponential") {
    # Few large effects, many small effects
    effects <- rexp(params$n_causal, rate = 1/params$effect_size)
    # Standardize to maintain overall effect size
    effects <- effects / sum(effects) * params$effect_size * params$n_causal
    
  } else if (params$effect_size_distribution == "normal") {
    effects <- abs(rnorm(params$n_causal, 
                        mean = params$effect_size, 
                        sd = params$effect_size * 0.3))
  }
  
  # Random signs
  effects <- effects * sample(c(-1, 1), params$n_causal, replace = TRUE)
  
  return(effects)
}

#' Calculate genetic value from genotypes
#'
#' @param snps SNP matrix
#' @param causal_loci Integer vector of causal locus positions
#' @param effects Numeric vector of effect sizes
#' @param params Simulation parameters
#' @return Numeric vector of genetic values
#' @keywords internal
calculate_genetic_value <- function(snps, causal_loci, effects, params) {
  
  causal_snps <- snps[, causal_loci, drop = FALSE]
  
  # Apply dominance if not additive
  if (params$dominance != 0.5) {
    # Convert to 0/1/2 coding if needed
    causal_snps <- apply(causal_snps, 2, function(x) {
      # Assuming input is 0/1/2
      ifelse(x == 1, params$dominance * 2, x)
    })
  }
  
  # Calculate breeding values
  genetic_value <- as.vector(causal_snps %*% effects)
  
  # Standardize to desired heritability
  genetic_value <- scale(genetic_value)[,1] * sqrt(params$holobiont_heritability)
  
  return(genetic_value)
}

#' Calculate environmental effect
#'
#' @param env Environmental variable
#' @param params Simulation parameters
#' @return Numeric vector of environmental effects
#' @keywords internal
calculate_environmental_effect <- function(env, params) {
  
  # Handle categorical environment
  if (is.character(env) || is.factor(env)) {
    # Convert to numeric codes
    env <- as.numeric(factor(env))
  }
  
  env_value <- scale(env)[,1] * sqrt(params$env_effect_on_holobiont)
  
  return(env_value)
}

#' Format response variable based on type
#'
#' @param holobiont Holobiont data object
#' @param response Raw response values
#' @param params Simulation parameters
#' @return Updated holobiont object
#' @keywords internal
format_response_variable <- function(holobiont, response, params) {
  
  if (params$response_variable == "abundance") {
    # Single taxon abundance (0-1)
    response <- plogis(response)
    
    # Add measurement error
    if (params$measurement_error > 0) {
      response <- response + rnorm(length(response), sd = params$measurement_error)
      response <- pmax(0, pmin(1, response))
    }
    
    holobiont$microbiome$target_taxon <- response
    
  } else if (params$response_variable == "diversity") {
    # Alpha diversity (count data)
    response <- params$baseline_diversity + response * 20
    response <- pmax(0, round(response))
    
    # Add measurement error
    if (params$measurement_error > 0) {
      response <- response + rnorm(length(response), 
                                   sd = params$measurement_error * params$baseline_diversity)
      response <- pmax(0, round(response))
    }
    
    holobiont$microbiome$alpha_diversity <- response
    
  } else if (params$response_variable == "composition") {
    # Beta diversity PC1 (continuous)
    
    # Add measurement error
    if (params$measurement_error > 0) {
      response <- response + rnorm(length(response), sd = params$measurement_error)
    }
    
    holobiont$microbiome$pc1 <- response
    
    # Optionally simulate full community matrix
    if (params$microbiome_structure == "compositional") {
      holobiont$microbiome$fungal_community <- 
        simulate_community_matrix(length(response), params$n_taxa, response, params)
    }
    
  } else if (params$response_variable == "phenotype") {
    
    # Add measurement error
    if (params$measurement_error > 0) {
      response <- response + rnorm(length(response), sd = params$measurement_error)
    }
    
    holobiont$phenotypes$trait <- response
  }
  
  return(holobiont)
}

#' Simulate full community matrix
#'
#' @param n Number of individuals
#' @param n_taxa Number of taxa
#' @param pc1_values PC1 values driving community composition
#' @param params Simulation parameters
#' @return Community matrix (individuals x taxa)
#' @keywords internal
simulate_community_matrix <- function(n, n_taxa, pc1_values, params) {
  
  # Generate base community with random structure
  comm <- matrix(rexp(n * n_taxa), nrow = n, ncol = n_taxa)
  
  # Make subset of taxa respond to PC1
  n_responsive <- round(n_taxa * 0.2)
  responsive_taxa <- sample(n_taxa, n_responsive)
  
  for (taxon in responsive_taxa) {
    effect <- rnorm(1, sd = 0.5)
    comm[, taxon] <- comm[, taxon] * exp(effect * scale(pc1_values)[,1])
  }
  
  # Convert to relative abundance (compositional data)
  comm <- comm / rowSums(comm)
  
  return(comm)
}