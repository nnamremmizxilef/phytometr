#' Simulate individual-level environmental values from an env_space
#'
#' Given an \code{env_space} object (created by \code{define_env_space()}),
#' this function simulates individual-level environmental values by
#' drawing around the realised group means with the specified
#' within-group standard deviations.
#'
#' This is the bridge between:
#' \itemize{
#'   \item environmental design (\code{env_space}: which sites/conditions
#'         are possible), and
#'   \item individual-level simulations (power analysis, holobiont responses).
#' }
#'
#' @param env_space Object of class \code{"env_space"} created by
#'   \code{define_env_space()}.
#' @param n_per_group Integer or integer vector. Number of individuals
#'   per selected group. If a single integer, it is recycled for all
#'   groups. If a vector, its length must equal the number of groups.
#' @param groups Optional vector of group IDs to use (a subset of
#'   \code{env_space$planned$group_id}). If \code{NULL} (default), all
#'   groups in \code{env_space} are used.
#' @param seed Optional integer. If provided, used to set the random
#'   seed for reproducibility.
#'
#' @return A data.frame with one row per simulated individual and
#'   columns:
#'   \itemize{
#'     \item \code{group_id}: site/condition index.
#'     \item one column per environmental variable (same names as in
#'           \code{env_space$variables}).
#'   }
#'
#' @examples
#' \dontrun{
#' env_space <- define_env_space(
#'   context   = "field",
#'   variables = list(
#'     temperature = list(
#'       range           = c(15, 25),
#'       shape           = "linear",
#'       target_sd       = 1,
#'       within_group_sd = 0.5
#'     ),
#'     moisture = list(
#'       range           = c(0.2, 0.6),
#'       shape           = "random",
#'       target_sd       = 0.05,
#'       within_group_sd = 0.1
#'     )
#'   ),
#'   n_groups = 6,
#'   seed     = 123
#' )
#'
#' ind_env <- simulate_individual_env(env_space, n_per_group = 10)
#' head(ind_env)
#' }
#'
#' @export
simulate_individual_env <- function(
  env_space,
  n_per_group,
  groups = NULL,
  seed = NULL
) {
  if (!inherits(env_space, "env_space")) {
    stop("'env_space' must be an object created by define_env_space().")
  }

  # Select groups
  all_groups <- env_space$planned$group_id
  if (is.null(groups)) {
    groups <- all_groups
  } else {
    if (!all(groups %in% all_groups)) {
      stop("Some 'groups' are not present in env_space$planned$group_id.")
    }
  }
  groups <- as.integer(groups)
  n_groups <- length(groups)

  # Handle n_per_group
  if (length(n_per_group) == 1L) {
    n_per_group <- rep(as.integer(n_per_group), n_groups)
  } else {
    if (length(n_per_group) != n_groups) {
      stop("'n_per_group' must be a single integer or have length equal to length(groups).")
    }
    n_per_group <- as.integer(n_per_group)
  }
  if (any(n_per_group < 1)) {
    stop("All entries in 'n_per_group' must be >= 1.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  var_names <- env_space$variables
  realised  <- env_space$realised
  within_sd <- env_space$within_group_sd

  # Prepare output container
  total_n <- sum(n_per_group)
  out <- data.frame(
    group_id = integer(total_n),
    matrix(NA_real_, nrow = total_n, ncol = length(var_names))
  )
  colnames(out)[-1] <- var_names

  row_index <- 1L

  for (i in seq_along(groups)) {
    g      <- groups[i]
    ng     <- n_per_group[i]
    g_row  <- realised[realised$group_id == g, , drop = FALSE]

    if (nrow(g_row) != 1L) {
      stop("Internal error: multiple rows for group_id ", g, " in env_space$realised.")
    }

    # For each variable, simulate NG values around realised mean
    for (v in var_names) {
      mu <- g_row[[v]]
      sd <- within_sd[[v]]
      if (is.null(sd) || is.na(sd)) {
        sd <- 0
      }
      out[row_index:(row_index + ng - 1), v] <- if (sd > 0) {
        stats::rnorm(ng, mean = mu, sd = sd)
      } else {
        rep(mu, ng)
      }
    }

    out$group_id[row_index:(row_index + ng - 1)] <- g
    row_index <- row_index + ng
  }

  return(out)
}
