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
simulate_holobiont_components <- function(components, env_data) {

  if (missing(components)) {
    stop("'components' must be provided and must be a named list.")
  }

  if (!is.list(components) || is.null(names(components))) {
    stop("'components' must be a *named* list.")
  }

  if (missing(env_data)) {
    stop("'env_data' must be provided and must contain environmental variables.")
  }

  env_data <- as.data.frame(env_data)
  n <- nrow(env_data)

  results <- data.frame(id = seq_len(n))

  # Loop over holobiont components
  for (comp in names(components)) {

    spec <- components[[comp]]

    if (!("r2_env" %in% names(spec))) {
      stop(sprintf("Component '%s' must contain 'r2_env'.", comp))
    }

    r2_vec <- spec$r2_env

    # Check that env predictors exist
    if (!all(names(r2_vec) %in% colnames(env_data))) {
      stop(sprintf("Environmental variables for '%s' not found in env_data", comp))
    }

    # Convert R² to beta coefficients: beta = sqrt(R²)
    betas <- sqrt(r2_vec)

    signal <-
      as.numeric(as.matrix(env_data[, names(r2_vec), drop = FALSE]) %*% betas)

    # Add environmental noise scaled so the resulting R² matches
    noise_sd <- sd(signal) * (1 - mean(r2_vec))
    noise <- rnorm(n, mean = 0, sd = noise_sd)

    results[[comp]] <- signal + noise
  }

  return(results)
}
