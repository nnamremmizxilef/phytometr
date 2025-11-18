#' Simulate a holobiont response driven by one or more environmental variables
#'
#' This function simulates a single univariate holobiont response
#' (e.g. microbial diversity, community PC1, trait value) such that
#' a specified proportion of its variance is explained by one or
#' several environmental variables.
#'
#' Typical use cases:
#' \itemize{
#'   \item "Temperature explains ~30% of variation in root bacterial diversity"
#'   \item "Temperature explains 20% and water availability 10% of fungal richness"
#' }
#' Given an environmental design and literature-based R^2 values,
#' the function returns a simulated response vector with that
#' approximate variance structure.
#'
#' @param env Numeric vector, matrix, or data.frame of environmental
#'   variables. Each column is a predictor (e.g. temperature, water).
#'   \itemize{
#'     \item A numeric vector means a single predictor.
#'     \item A matrix/data.frame means multiple predictors.
#'   }
#' @param r2_env Either:
#'   \itemize{
#'     \item a single numeric in [0, 1]: total proportion of response
#'           variance explained by all environmental variables together
#'           (shared equally across columns), or
#'     \item a numeric vector in [0, 1] of length equal to \code{ncol(env)},
#'           optionally named to match \code{colnames(env)}, giving the
#'           variance contribution of each environmental variable.
#'   }
#'   For example, if \code{env} has columns \code{temp} and \code{water},
#'   you can use \code{r2_env = c(temp = 0.3, water = 0.1)}.
#' @param total_var Numeric > 0. Total variance of the simulated response.
#'   Default is 1. The environmental and residual components are scaled
#'   so that their variances sum to \code{total_var}.
#' @param seed Optional integer. If provided, used to set the random seed
#'   for reproducibility.
#'
#' @return Numeric vector of length \code{length(env)} (if \code{env} is
#'   a vector) or \code{nrow(env)} (if \code{env} is a matrix/data.frame),
#'   representing the simulated holobiont response.
#'
#' @examples
#' # Single environmental variable: temperature
#' set.seed(1)
#' temp <- seq(10, 20, length.out = 100)
#' resp <- simulate_env_response(temp, r2_env = 0.3)
#' summary(lm(resp ~ temp))$r.squared  # roughly ~0.3
#'
#' # Two variables: temperature and water
#' set.seed(1)
#' n <- 120
#' temp  <- seq(15, 25, length.out = n)
#' water <- seq(0.3, 0.9, length.out = n)
#' env <- cbind(temp = temp, water = water)
#'
#' resp2 <- simulate_env_response(
#'   env    = env,
#'   r2_env = c(temp = 0.2, water = 0.1)  # ~30% total explained by env
#' )
#' summary(lm(resp2 ~ temp + water))$r.squared
#'
#' @export
simulate_env_response <- function(
  env,
  r2_env = 0.3,
  total_var = 1,
  seed = NULL
) {
  # --- Standardise env input to a numeric matrix ---
  if (is.vector(env) && !is.list(env)) {
    # Single environmental variable
    env <- as.numeric(env)
    n <- length(env)
    env_mat <- matrix(env, ncol = 1)
    colnames(env_mat) <- "env1"
  } else {
    # Multiple environmental variables
    env_mat <- as.matrix(env)
    n <- nrow(env_mat)
  }

  if (!is.numeric(env_mat)) {
    stop("'env' must be numeric (vector, matrix, or data.frame).")
  }
  if (n < 2L) {
    stop("'env' must have at least 2 observations.")
  }

  p <- ncol(env_mat)  # number of environmental predictors

  # --- Interpret r2_env ---
  if (!is.numeric(r2_env) || any(r2_env < 0) || any(r2_env > 1)) {
    stop("'r2_env' must be numeric and in [0, 1].")
  }

  if (length(r2_env) == 1L) {
    # User gave total R^2 for all env variables together:
    # split it equally across the p variables
    r2_vec <- rep(r2_env / p, p)
  } else if (length(r2_env) == p) {
    # User gave one R^2 per env variable
    r2_vec <- as.numeric(r2_env)
  } else {
    stop("If 'env' has ", p, " columns, 'r2_env' must have length 1 or ", p, ".")
  }

  if (!is.numeric(total_var) || length(total_var) != 1L || total_var <= 0) {
    stop("'total_var' must be a single positive number.")
  }

  # Total R^2 contributed by environment variables
  total_r2_env <- sum(r2_vec)
  if (total_r2_env > 1 + 1e-8) {
    stop("Sum of 'r2_env' components cannot exceed 1.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Build environmental contribution to the response ---
  env_effect <- numeric(n)  # start with all zeros

  for (j in seq_len(p)) {
    env_j <- env_mat[, j]

    # Center the variable so intercept = 0 on average
    env_j_centered <- env_j - mean(env_j, na.rm = TRUE)
    var_env_j <- stats::var(env_j_centered, na.rm = TRUE)

    # If this environmental variable is constant, skip it
    if (var_env_j == 0) {
      next
    }

    # Desired variance contribution from this variable
    var_env_component_j <- r2_vec[j] * total_var

    if (var_env_component_j > 0) {
      # Find slope beta_j such that:
      #   var(beta_j * env_j_centered) = var_env_component_j
      beta_j <- sqrt(var_env_component_j / var_env_j)
      env_effect <- env_effect + beta_j * env_j_centered
    }
  }

  # --- Add residual noise to reach total_var ---
  var_env_total <- stats::var(env_effect)

  # Numerical safety: clamp between 0 and total_var
  var_env_total <- max(min(var_env_total, total_var), 0)

  var_residual <- total_var - var_env_total
  if (var_residual < 0) {
    var_residual <- 0
  }

  residual <- stats::rnorm(n, mean = 0, sd = sqrt(var_residual))

  # Final response = env-driven signal + noise
  response <- env_effect + residual

  return(response)
}
