#' Simulate environmental data
#'
#' @param holobiont A holobiont_data object
#' @param params List of simulation parameters
#' @return Updated holobiont_data object with environment data
#' @keywords internal
simulate_environment_data <- function(holobiont, params) {
  
  n <- length(holobiont$host_id)
  
  if (params$env_type == "continuous") {
    
    env <- generate_continuous_environment(n, holobiont, params)
    holobiont$environment$main <- env
    
  } else if (params$env_type == "categorical") {
    
    if (is.null(params$env_categories)) {
      params$env_categories <- c("low", "medium", "high")
    }
    holobiont$environment$main <- sample(
      params$env_categories, 
      n, 
      replace = TRUE
    )
    
  } else if (params$env_type == "multivariate") {
    
    env_data <- generate_multivariate_environment(n, params)
    
    for (i in seq_len(params$env_n_variables)) {
      holobiont$environment[[paste0("env_", i)]] <- env_data[, i]
    }
    
    # Set first variable as main
    holobiont$environment$main <- env_data[, 1]
  }
  
  return(holobiont)
}

#' Generate continuous environmental gradient
#'
#' @param n Number of individuals
#' @param holobiont Holobiont data object
#' @param params Simulation parameters
#' @return Numeric vector of environmental values
#' @keywords internal
generate_continuous_environment <- function(n, holobiont, params) {
  
  if (params$gradient_type == "linear") {
    # Linear gradient (e.g., latitude)
    env <- seq(params$env_range[1], params$env_range[2], length.out = n)
    
  } else if (params$gradient_type == "nonlinear") {
    # Nonlinear gradient (quadratic)
    x <- seq(0, 1, length.out = n)
    env <- params$env_range[1] + 
          diff(params$env_range) * x^2
    
  } else if (params$gradient_type == "patchy") {
    # Patchy environment with discrete patches
    if (!is.null(holobiont$population)) {
      # Use existing population structure
      n_patches <- length(unique(holobiont$population))
      patch_means <- seq(params$env_range[1], params$env_range[2], 
                        length.out = n_patches)
      env <- patch_means[holobiont$population]
    } else {
      n_patches <- params$n_populations
      patch_means <- seq(params$env_range[1], params$env_range[2], 
                        length.out = n_patches)
      env <- rep(patch_means, each = ceiling(n / n_patches))[seq_len(n)]
    }
    
  } else if (params$gradient_type == "random") {
    # Random environmental variation
    env <- runif(n, params$env_range[1], params$env_range[2])
  }
  
  # Add noise
  if (params$gradient_noise > 0) {
    noise_sd <- diff(params$env_range) * params$gradient_noise
    env <- env + rnorm(n, sd = noise_sd)
  }
  
  return(env)
}

#' Generate multivariate environmental data
#'
#' @param n Number of individuals
#' @param params Simulation parameters
#' @return Matrix of environmental variables (n x n_variables)
#' @keywords internal
generate_multivariate_environment <- function(n, params) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multivariate environment simulation")
  }
  
  # Build correlation matrix
  cor_matrix <- matrix(params$env_correlation, 
                      params$env_n_variables, 
                      params$env_n_variables)
  diag(cor_matrix) <- 1
  
  # Convert to covariance matrix
  env_sd <- diff(params$env_range) / 4  # Approximate SD
  cov_matrix <- cor_matrix * env_sd^2
  
  # Generate correlated variables
  env_matrix <- MASS::mvrnorm(
    n = n,
    mu = rep(mean(params$env_range), params$env_n_variables),
    Sigma = cov_matrix
  )
  
  return(env_matrix)
}