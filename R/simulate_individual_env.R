#' Simulate individual-level environmental values
#'
#' Given an `env_space` object (from `define_env_space()`,
#' `define_env_space_grid()`, or `define_env_locations()`), generate
#' individual-level environmental conditions by adding within-group
#' noise around the group means.
#'
#' @param env_space An object of class `env_space`.
#' @param n_per_group Integer. Either a single number (same N for all
#'   groups) or a vector of length `env_space$n_groups` giving the
#'   number of individuals per group.
#'
#' @return A data frame with columns
#'   `group_id` and one column per environmental variable.
#' @export
simulate_individual_env <- function(env_space, n_per_group) {

  if (!inherits(env_space, "env_space")) {
    stop("`env_space` must be an object of class 'env_space'.")
  }

  realised  <- env_space$realised
  variables <- env_space$variables
  n_groups  <- env_space$n_groups
  within_sd <- env_space$within_group_sd

  if (!("group_id" %in% names(realised))) {
    stop("`env_space$realised` must contain a 'group_id' column.")
  }

  # allow single n_per_group or per-group vector
  if (length(n_per_group) == 1L) {
    n_per_group <- rep(n_per_group, n_groups)
  } else if (length(n_per_group) != n_groups) {
    stop("`n_per_group` must be length 1 or length env_space$n_groups (",
         n_groups, ").")
  }

  # sanity checks
  if (!all(variables %in% names(realised))) {
    missing <- setdiff(variables, names(realised))
    stop("The following variables listed in env_space$variables ",
         "are missing from env_space$realised: ",
         paste(missing, collapse = ", "))
  }

  if (length(within_sd) == 1L) {
    within_sd <- rep(within_sd, length(variables))
    names(within_sd) <- variables
  }

  if (!all(variables %in% names(within_sd))) {
    missing <- setdiff(variables, names(within_sd))
    stop("within_group_sd must have entries for all variables. Missing: ",
         paste(missing, collapse = ", "))
  }

  # we keep group labels as they are (numeric, character, factor)
  group_labels <- realised$group_id

  if (length(group_labels) != n_groups) {
    stop("Internal error: length(group_labels) != env_space$n_groups.")
  }

  out_list <- vector("list", n_groups)

  for (i in seq_len(n_groups)) {
    g_label <- group_labels[i]

    # exactly one row per group_id in realised
    row_i <- realised[i, , drop = FALSE]

    n_i <- n_per_group[i]

    # simulate env variables for individuals in group i
    sim_env <- lapply(variables, function(v) {
      mu <- row_i[[v]]
      sd <- within_sd[[v]]
      stats::rnorm(n_i, mean = mu, sd = sd)
    })

    sim_df <- as.data.frame(sim_env)
    names(sim_df) <- variables

    sim_df$group_id <- rep(g_label, n_i)
    # put group_id first
    sim_df <- sim_df[, c("group_id", variables), drop = FALSE]

    out_list[[i]] <- sim_df
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
