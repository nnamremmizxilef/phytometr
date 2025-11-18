#' Simulate multiple holobiont components driven by environment
#'
#' This function simulates several holobiont response variables
#' (e.g. root bacteria, root fungi, herbivores, "whole holobiont")
#' from one or more environmental predictors. Each component can
#' have its own literature-based effect sizes (R² values) for
#' the environmental variables.
#'
#' Under the hood, it calls \code{simulate_env_response()} once
#' for each component and returns all simulated responses in a
#' single data frame.
#'
#' @param env Numeric vector, matrix, or data.frame of environmental
#'   variables. Each column is a predictor (e.g. temperature, water).
#'   This is passed directly to \code{simulate_env_response()}.
#' @param components Named list describing each holobiont component.
#'   Each element should be a list with at least:
#'   \itemize{
#'     \item \code{r2_env}: either a single numeric (total R² for this
#'           component) or a numeric vector per environmental variable,
#'           as in \code{simulate_env_response()}.
#'     \item \code{total_var}: optional, total variance of the component
#'           (default 1 if not provided).
#'   }
#'
#'   Example:
#'   \preformatted{
#'   components <- list(
#'     root_bacteria = list(
#'       r2_env    = c(temp = 0.30, water = 0.10),
#'       total_var = 1
#'     ),
#'     root_fungi = list(
#'       r2_env    = c(temp = 0.20),
#'       total_var = 1
#'     ),
#'     herbivores = list(
#'       r2_env    = c(temp = 0.05, water = 0.25)
#'     )
#'   )
#'   }
#'
#' @param seed Optional integer. If provided, used as a base seed.
#'   Each component will internally use \code{seed + i} (for component
#'   i) to remain reproducible but independent.
#'
#' @return A data.frame with one column per component and as many
#'   rows as there are observations in \code{env}.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' temp  <- seq(15, 25, length.out = n)
#' water <- seq(0.3, 0.9, length.out = n)
#' env   <- cbind(temp = temp, water = water)
#'
#' components <- list(
#'   root_bacteria = list(
#'     r2_env    = c(temp = 0.30, water = 0.10),
#'     total_var = 1
#'   ),
#'   root_fungi = list(
#'     r2_env    = c(temp = 0.20, water = 0.05),
#'     total_var = 1
#'   )
#' )
#'
#' holo <- simulate_holobiont_components(env, components, seed = 123)
#' str(holo)
#' summary(lm(root_bacteria ~ temp + water, data = cbind(holo, env)))
#'
#' @export
simulate_holobiont_components <- function(
  env,
  components,
  seed = NULL
) {
  # --- Basic checks on components ---
  if (!is.list(components) || is.null(names(components))) {
    stop("'components' must be a *named* list.")
  }

  comp_names <- names(components)
  if (any(comp_names == "")) {
    stop("All entries in 'components' must have a non-empty name.")
  }

  # Determine number of observations (n) from env
  if (is.vector(env) && !is.list(env)) {
    n <- length(env)
  } else {
    env <- as.matrix(env)
    n <- nrow(env)
  }

  # Prepare output data.frame
  result <- as.data.frame(matrix(NA_real_, nrow = n, ncol = length(components)))
  colnames(result) <- comp_names

  # Optional: set a base seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Loop over components and simulate each one ---
  for (i in seq_along(components)) {
    comp_name <- comp_names[i]
    comp_spec <- components[[i]]

    # Extract r2_env and total_var, with defaults
    if (is.null(comp_spec$r2_env)) {
      stop("Component '", comp_name, "' must have an 'r2_env' element.")
    }
    r2_env    <- comp_spec$r2_env
    total_var <- if (!is.null(comp_spec$total_var)) comp_spec$total_var else 1

    # For reproducibility: different seed for each component
    comp_seed <- if (!is.null(seed)) seed + i else NULL

    # Call lower-level simulator for this component
    result[[comp_name]] <- simulate_env_response(
      env       = env,
      r2_env    = r2_env,
      total_var = total_var,
      seed      = comp_seed
    )
  }

  return(result)
}
