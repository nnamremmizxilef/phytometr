#' Simulate holobiont components from environmental predictors
#'
#' @param components A named list. Each entry defines one holobiont component
#'   and must contain an element \code{r2_env}, which is a named numeric vector
#'   giving the (approximate) variance contribution of each environmental
#'   variable.
#' @param env_data A data frame or matrix with one row per individual (or group)
#'   and one column per environmental variable. If a column called
#'   \code{group_id} is present, it will be carried through to the output so
#'   that components can be grouped in plots.
#'
#' @return A data frame with:
#'   \itemize{
#'     \item \code{id} – row index (1..n),
#'     \item optional \code{group_id} (if present in \code{env_data}),
#'     \item one column per holobiont component.
#'   }
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
  n        <- nrow(env_data)

  has_group <- "group_id" %in% names(env_data)

  # we will only use non-group env vars for simulation
  env_vars <- setdiff(colnames(env_data), "group_id")

  res <- data.frame(id = seq_len(n))

  for (comp_name in names(components)) {

    spec <- components[[comp_name]]

    if (!("r2_env" %in% names(spec))) {
      stop(sprintf("Component '%s' must contain 'r2_env'.", comp_name))
    }

    r2_vec <- spec$r2_env

    # Check variables present in env_data
    if (!all(names(r2_vec) %in% env_vars)) {
      missing_vars <- setdiff(names(r2_vec), env_vars)
      stop(
        sprintf(
          "For component '%s', env variables not found in env_data: %s",
          comp_name,
          paste(missing_vars, collapse = ", ")
        )
      )
    }

    # Extract and standardise relevant env variables
    env_sub <- scale(env_data[, names(r2_vec), drop = FALSE])

    # Total R² explained by all env variables together
    total_r2 <- sum(r2_vec)
    if (total_r2 <= 0) {
      warning(
        sprintf(
          "Component '%s' has non-positive total R² (sum(r2_env) = %g); using noise only.",
          comp_name, total_r2
        )
      )
      res[[comp_name]] <- stats::rnorm(n)
      next
    }
    if (total_r2 >= 1) {
      warning(
        sprintf(
          "Component '%s' has total R² >= 1 (sum(r2_env) = %g); truncating to 0.99.",
          comp_name, total_r2
        )
      )
      total_r2 <- 0.99
    }

    # Relative weights across env variables
    weights <- r2_vec / sum(r2_vec)

    # Latent environmental score (mean 0, var ~1)
    env_score <- drop(as.matrix(env_sub) %*% weights)

    # Independent noise
    epsilon <- stats::rnorm(n)

    # Construct response: cor(env_score, Y)^2 ≈ total_r2
    y_latent <- sqrt(total_r2) * env_score + sqrt(1 - total_r2) * epsilon

    res[[comp_name]] <- y_latent
  }

  # Attach group_id if present
  if (has_group) {
    res$group_id <- env_data$group_id
    # reorder columns: id, group_id, components...
    res <- res[, c("id", "group_id", setdiff(names(res), c("id", "group_id")))]
  }

  return(res)
}
