#' Define a multivariate environmental space (several variables)
#'
#' This function defines a set of potential sites/conditions along
#' multiple environmental axes (e.g. temperature, moisture, pH),
#' without fixing the number of individuals. It calls
#' \code{define_env_gradient()} once for each variable and combines
#' the results.
#'
#' For each environmental variable, you specify:
#' \itemize{
#'   \item a feasible range (\code{range}),
#'   \item a gradient shape (\code{shape}),
#'   \item how precisely you can hit the target mean (\code{target_sd}),
#'   \item expected within-group variation (\code{within_group_sd}).
#' }
#'
#' All variables share the same number of groups (\code{n_groups}) and
#' context (\code{"field"} or \code{"controlled"}).
#'
#' @param context Character, either \code{"field"} or \code{"controlled"}.
#'   Describes whether this is a field gradient or a controlled experiment.
#' @param variables Named list specifying each environmental variable.
#'   Each element must be a list with at least a \code{range} element,
#'   and optionally \code{shape}, \code{target_sd}, and
#'   \code{within_group_sd}. Example:
#'   \preformatted{
#'   variables <- list(
#'     temperature = list(
#'       range          = c(15, 25),
#'       shape          = "linear",
#'       target_sd      = 1,
#'       within_group_sd = 0.5
#'     ),
#'     moisture = list(
#'       range          = c(0.2, 0.6),
#'       shape          = "random",
#'       target_sd      = 0.05,
#'       within_group_sd = 0.1
#'     )
#'   )
#'   }
#' @param n_groups Integer. Number of potential sites/conditions along
#'   each environmental gradient.
#' @param seed Optional integer. If provided, used as a base seed.
#'   Each variable will internally use \code{seed + i} for reproducibility.
#'
#' @return A list with class \code{"env_space"} containing:
#'   \itemize{
#'     \item \code{context}: field or controlled.
#'     \item \code{n_groups}: number of groups.
#'     \item \code{variables}: character vector of variable names.
#'     \item \code{planned}: data.frame with columns \code{group_id} and
#'           one column per variable giving the planned group means.
#'     \item \code{realised}: data.frame with columns \code{group_id} and
#'           one column per variable giving realised group means.
#'     \item \code{within_group_sd}: named numeric vector of within-group SDs
#'           for each variable.
#'   }
#'
#'   You can obtain a matrix of realised group means for use in
#'   holobiont simulations via:
#'   \code{as.matrix(env_space$realised[,-1])}.
#'
#' @examples
#' # Define multivariate field environment: temperature + moisture
#' vars <- list(
#'   temperature = list(
#'     range          = c(15, 25),
#'     shape          = "linear",
#'     target_sd      = 1,
#'     within_group_sd = 0.5
#'   ),
#'   moisture = list(
#'     range          = c(0.2, 0.6),
#'     shape          = "random",
#'     target_sd      = 0.05,
#'     within_group_sd = 0.1
#'   )
#' )
#'
#' env_space <- define_env_space(
#'   context  = "field",
#'   variables = vars,
#'   n_groups = 6,
#'   seed     = 123
#' )
#'
#' env_space$planned
#' env_space$realised
#'
#' # Use realised group means (e.g. as env for holobiont simulation)
#' realised_env_mat <- as.matrix(env_space$realised[,-1])
#'
#' @export
define_env_space <- function(
  context   = c("field", "controlled"),
  variables,
  n_groups  = 5L,
  seed      = NULL
) {
  context <- match.arg(context)

  # Check variables list
  if (!is.list(variables) || is.null(names(variables))) {
    stop("'variables' must be a *named* list.")
  }
  var_names <- names(variables)
  if (any(var_names == "")) {
    stop("All entries in 'variables' must have non-empty names.")
  }

  if (!is.numeric(n_groups) || length(n_groups) != 1L || n_groups < 1) {
    stop("'n_groups' must be a positive integer.")
  }
  n_groups <- as.integer(n_groups)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Containers for outputs
  planned_df  <- data.frame(group_id = seq_len(n_groups))
  realised_df <- data.frame(group_id = seq_len(n_groups))
  within_sd_vec <- numeric(length(variables))
  names(within_sd_vec) <- var_names

  # Loop over variables and call define_env_gradient() for each
  for (i in seq_along(variables)) {
    v_name <- var_names[i]
    v_spec <- variables[[i]]

    if (is.null(v_spec$range) || length(v_spec$range) != 2L) {
      stop("Variable '", v_name, "' must have a 'range' element of length 2.")
    }

    v_shape           <- if (!is.null(v_spec$shape)) v_spec$shape else "linear"
    v_target_sd       <- if (!is.null(v_spec$target_sd)) v_spec$target_sd else 0
    v_within_group_sd <- if (!is.null(v_spec$within_group_sd)) v_spec$within_group_sd else 0

    # Different seed per variable if base seed is provided
    v_seed <- if (!is.null(seed)) seed + i else NULL

    grad <- define_env_gradient(
      context         = context,
      env_range       = v_spec$range,
      n_groups        = n_groups,
      shape           = v_shape,
      target_sd       = v_target_sd,
      within_group_sd = v_within_group_sd,
      seed            = v_seed
    )

    planned_df[[v_name]]  <- grad$planned_env
    realised_df[[v_name]] <- grad$realised_env
    within_sd_vec[v_name] <- v_within_group_sd
  }

  out <- list(
    context         = context,
    n_groups        = n_groups,
    variables       = var_names,
    planned         = planned_df,
    realised        = realised_df,
    within_group_sd = within_sd_vec
  )
  class(out) <- "env_space"

  return(out)
}
