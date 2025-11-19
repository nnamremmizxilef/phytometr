#' Simulate holobiont components at the group level
#'
#' This function generates *group-level* holobiont component values
#' (e.g. mean root fungal diversity per site) from environmental
#' data and a component specification.
#'
#' @param components A named list of component specifications. Each
#'   entry must be a list with at least:
#'   \itemize{
#'     \item \code{r2_env}: named numeric vector of non-negative values
#'           giving the proportion of variance explained by each
#'           environmental variable (R² per variable). Names must match
#'           columns in \code{env_data}.
#'     \item \code{sign_env} (optional): named numeric vector of the
#'           same length as \code{r2_env}, with values typically in
#'           \code{-1, 0, 1}. This encodes *direction* of the effect:
#'           negative values yield negative correlations with the
#'           corresponding environmental variable. If omitted, all
#'           signs default to +1.
#'   }
#' @param env_data A data.frame with one row per group (e.g. site or
#'   condition). Must contain the environmental variables referenced
#'   in \code{r2_env}. If a column \code{group_id} is present, it is
#'   carried along in the output; otherwise groups are numbered 1..N.
#' @param total_var Numeric. Target total variance for each component
#'   (after scaling). Default is 1. Internally, values are scaled to
#'   have approximately unit variance, so the absolute scale is less
#'   important than the *relative contributions* of environment and
#'   noise.
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item \code{id}: group index
#'     \item \code{group_id}: copied from \code{env_data$group_id} if present
#'     \item one column per component (e.g. \code{root_fungi_diversity})
#'   }
#'
#' @details
#' The environmental contribution for a component is constructed as:
#' \deqn{
#'   \text{signal} = \sum_j \text{sign}_j \sqrt{R^2_j} \cdot \tilde{E}_j,
#' }
#' where \eqn{\tilde{E}_j} is the standardized environmental variable
#' and \eqn{R^2_j \ge 0} is the variance contribution of variable \eqn{j}.
#' The sign vector \code{sign_env} controls whether the correlation
#' with a given environmental variable is positive or negative.
#'
#' The remaining variance is filled with Gaussian noise. The output
#' columns are finally scaled to have approximately unit variance.
#'
#' @export
simulate_holobiont_components <- function(
  components,
  env_data,
  total_var = 1
) {
  if (!is.list(components) || length(components) == 0L) {
    stop("'components' must be a non-empty list.")
  }
  if (!is.data.frame(env_data)) {
    stop("'env_data' must be a data.frame.")
  }

  n_groups <- nrow(env_data)
  if (n_groups < 1L) {
    stop("'env_data' must have at least one row.")
  }

  # group id handling
  if ("group_id" %in% colnames(env_data)) {
    group_id <- env_data$group_id
  } else {
    group_id <- seq_len(n_groups)
  }

  out <- data.frame(
    id       = seq_len(n_groups),
    group_id = group_id
  )

  # available env variables for modelling
  env_vars_all <- setdiff(colnames(env_data), "group_id")

  for (comp_name in names(components)) {
    spec <- components[[comp_name]]

    if (is.null(spec$r2_env)) {
      stop("Component '", comp_name, "' must contain 'r2_env'.")
    }

    r2_env <- spec$r2_env

    if (any(is.na(r2_env))) {
      stop("Component '", comp_name, "': 'r2_env' contains NA.")
    }
    if (any(r2_env < 0)) {
      stop("Component '", comp_name,
           "': all entries of 'r2_env' must be >= 0. ",
           "Use 'sign_env' to encode negative effects.")
    }

    # match environmental variables
    env_names <- names(r2_env)
    if (!all(env_names %in% env_vars_all)) {
      missing_vars <- env_names[!env_names %in% env_vars_all]
      stop("Component '", comp_name,
           "': environmental variables not found in env_data: ",
           paste(missing_vars, collapse = ", "))
    }

    # sign vector
    if (!is.null(spec$sign_env)) {
      sign_env <- spec$sign_env
      if (!all(names(sign_env) %in% env_names) ||
          length(sign_env) != length(r2_env)) {
        stop("Component '", comp_name, "': 'sign_env' must have the same ",
             "names and length as 'r2_env'.")
      }
      sign_env <- sign_env[env_names]
    } else {
      # default: all positive relationships
      sign_env <- rep(1, length(r2_env))
      names(sign_env) <- env_names
    }

    total_r2 <- sum(r2_env)

    if (total_r2 <= 0) {
      warning("Component '", comp_name, "' has total R2 <= 0; using pure noise.")
      signal <- rep(0, n_groups)
      noise_sd <- sqrt(total_var)
    } else if (total_r2 >= 1) {
      warning("Component '", comp_name,
              "' has total R2 >= 1; capping at 0.99 for simulation.")
      total_r2 <- 0.99
      r2_env <- r2_env * (total_r2 / sum(r2_env))
      # sign_env unchanged
    }

    if (total_r2 > 0) {
      # standardize env vars
      env_mat <- as.matrix(env_data[, env_names, drop = FALSE])
      env_scaled <- scale(env_mat)

      # weights = sign * sqrt(R2_j)
      weights <- sign_env * sqrt(r2_env)
      signal <- as.numeric(env_scaled %*% weights)

      # environment explains total_r2 of the variance
      noise_sd <- sqrt(total_var * (1 - total_r2))
    }

    noise <- stats::rnorm(n_groups, mean = 0, sd = noise_sd)

    comp_vals <- signal + noise

    # rescale to approx unit variance for convenience
    comp_vals <- scale(comp_vals)[, 1]

    out[[comp_name]] <- comp_vals
  }

  out
}
