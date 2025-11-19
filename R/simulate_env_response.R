#' Simulate holobiont responses at the individual level
#'
#' This function simulates individual-level holobiont component values
#' (e.g. diversity, abundance, herbivore load) from:
#' \itemize{
#'   \item environmental data per individual,
#'   \item a host model (clone or natural),
#'   \item and a component specification with R² and signs.
#' }
#'
#' @param env_data A data.frame with one row per individual. Must contain
#'   environmental variables referenced in each component's \code{r2_env}.
#'   A column \code{group_id} is optional (used for plotting/design).
#' @param components A named list in the same format as for
#'   \code{simulate_holobiont_components()}, i.e. each element is a list
#'   with at least \code{r2_env} and optional \code{sign_env}.
#' @param host_model An object of class \code{"host_model"} created by
#'   \code{define_host_model()}. Only \code{type = "clone"} and its
#'   \code{var_host} are currently used; natural host-variance partitions
#'   can be added later.
#' @param total_var Numeric. Total target variance per component (default 1).
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item \code{id}: individual index
#'     \item \code{group_id}: if present in \code{env_data}
#'     \item one column per holobiont component
#'   }
#'
#' @export
simulate_env_response <- function(
  env_data,
  components,
  host_model,
  total_var = 1
) {
  if (!is.data.frame(env_data)) {
    stop("'env_data' must be a data.frame.")
  }
  if (!is.list(components) || length(components) == 0L) {
    stop("'components' must be a non-empty list.")
  }
  if (!inherits(host_model, "host_model")) {
    stop("'host_model' must be created by define_host_model().")
  }

  n_ind <- nrow(env_data)
  if (n_ind < 1L) {
    stop("'env_data' must have at least one row (individual).")
  }

  # ids & group ids
  if ("id" %in% colnames(env_data)) {
    id <- env_data$id
  } else {
    id <- seq_len(n_ind)
  }

  if ("group_id" %in% colnames(env_data)) {
    group_id <- env_data$group_id
  } else {
    group_id <- NA_integer_
  }

  out <- data.frame(
    id       = id,
    group_id = group_id
  )

  env_vars_all <- setdiff(colnames(env_data), c("id", "group_id"))

  # host variance fraction (Phase 1: clone => ~0)
  var_host <- if (host_model$type == "clone") {
    0
  } else {
    # For now, take the var_host slot directly; you can refine later
    host_model$var_host %||% 0
  }

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
      sign_env <- rep(1, length(r2_env))
      names(sign_env) <- env_names
    }

    total_r2_env <- sum(r2_env)

    if (total_r2_env < 0 || total_r2_env + var_host > 1) {
      warning("Component '", comp_name,
              "': env R2 + host variance > 1; rescaling env R2.")
      # simple rescale to fit into (1 - var_host)
      if (total_r2_env > 0) {
        r2_env <- r2_env * ((1 - var_host) / total_r2_env)
        total_r2_env <- sum(r2_env)
      }
    }

    # ENVIRONMENTAL PART
    if (total_r2_env <= 0) {
      signal_env <- rep(0, n_ind)
    } else {
      env_mat <- as.matrix(env_data[, env_names, drop = FALSE])
      env_scaled <- scale(env_mat)
      weights <- sign_env * sqrt(r2_env)
      signal_env <- as.numeric(env_scaled %*% weights)
      var_env_target <- total_var * total_r2_env
      if (stats::var(signal_env) > 0) {
        signal_env <- signal_env / sqrt(stats::var(signal_env)) *
          sqrt(var_env_target)
      }
    }

    # HOST PART (placeholder for future Phase 2; clone => zero)
    if (var_host > 0) {
      host_effect <- stats::rnorm(n_ind, sd = sqrt(total_var * var_host))
    } else {
      host_effect <- rep(0, n_ind)
    }

    # RESIDUAL PART
    var_resid <- total_var - (total_var * total_r2_env) - (total_var * var_host)
    if (var_resid < 0) var_resid <- 0
    noise <- stats::rnorm(n_ind, sd = sqrt(var_resid))

    comp_vals <- signal_env + host_effect + noise

    # final scaling to sd ~1 for convenience
    if (stats::var(comp_vals) > 0) {
      comp_vals <- scale(comp_vals)[, 1]
    }

    out[[comp_name]] <- comp_vals
  }

  out
}
