#' Define factorial environmental space (grid of combinations)
#'
#' This function creates an `env_space` object where each group
#' corresponds to a combination of levels across environmental
#' variables (e.g. low/high temperature × low/high moisture).
#'
#' @param variables A named list, one element per environmental
#'   variable, each containing at least `range = c(min, max)`.
#'   Example:
#'   `list(temperature = list(range = c(15, 30)),
#'         moisture    = list(range = c(20, 80)))`
#'
#'   Optionally, each variable list can also contain
#'   `n_levels` (integer) and `within_group_sd` (numeric).
#'
#' @param n_levels Either a single integer (same number of levels
#'   for all variables) or a named integer vector with one value per
#'   variable (names must match `variables`).
#'   Ignored for variables that have their own `n_levels` entry.
#'
#' @param context Character, `"field"` or `"controlled"`.
#'
#' @param default_within_sd Numeric. Default within-group SD used
#'   if neither the variable nor `within_group_sd` specifies a value.
#'
#' @return An `env_space` object.
#' @export
define_env_space_grid <- function(variables,
                                  n_levels = 2,
                                  context = c("field", "controlled"),
                                  default_within_sd = 0) {

  context   <- match.arg(context)
  var_names <- names(variables)

  if (is.null(var_names) || any(var_names == "")) {
    stop("`variables` must be a *named* list, one element per variable.")
  }

  # helper: get levels per variable
  get_levels_for_var <- function(vname) {
    v <- variables[[vname]]

    if (is.null(v$range) || length(v$range) != 2) {
      stop("Variable '", vname, "' must have a numeric `range = c(min, max)`.")
    }

    # priority: variable-specific n_levels, then global n_levels
    if (!is.null(v$n_levels)) {
      n <- as.integer(v$n_levels)
    } else if (length(n_levels) == 1L) {
      n <- as.integer(n_levels)
    } else {
      if (is.null(names(n_levels)) ||
          !(vname %in% names(n_levels))) {
        stop("`n_levels` must be a single integer or a named vector ",
             "containing '", vname, "'.")
      }
      n <- as.integer(n_levels[[vname]])
    }

    seq(from = v$range[1],
        to   = v$range[2],
        length.out = n)
  }

  # build levels list (for expand.grid)
  levels_list <- lapply(var_names, get_levels_for_var)
  names(levels_list) <- var_names

  # factorial combination of all levels
  grid <- do.call(expand.grid, c(levels_list, KEEP.OUT.ATTRS = FALSE))
  grid$group_id <- seq_len(nrow(grid))
  # reorder: group_id first
  grid <- grid[ , c("group_id", var_names), drop = FALSE]

  # within-group SD per variable
  within_sd <- sapply(var_names, function(vname) {
    v <- variables[[vname]]
    if (!is.null(v$within_group_sd)) {
      return(as.numeric(v$within_group_sd))
    }
    default_within_sd
  })
  names(within_sd) <- var_names

  res <- list(
    context         = context,
    n_groups        = nrow(grid),
    variables       = var_names,
    planned         = grid,
    realised        = grid,
    within_group_sd = within_sd
  )
  class(res) <- "env_space"
  res
}
