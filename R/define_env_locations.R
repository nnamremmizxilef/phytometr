#' Define environmental conditions from explicit locations
#'
#' This is for the case where you already know which sites or
#' experimental conditions are available (e.g. 3 field sites, 4
#' growth chambers) and their typical environmental values.
#'
#' The function converts those into an `env_space` object that can be
#' used by the rest of the phytometr workflow.
#'
#' @param locations A data frame with one row per site/condition.
#'   Must contain either a column `group_id` or any of
#'   `location`, `site`, `condition`. All other numeric columns are
#'   treated as environmental variables (e.g. temperature, moisture).
#' @param context Character, `"field"` or `"controlled"`.
#'   Only used as a label in the resulting object.
#' @param within_group_sd Either
#'   * a single numeric (same SD for all variables), or
#'   * a named numeric vector with one SD per variable
#'     (names must match the environmental column names).
#'   This describes how much individual-level variation you expect
#'   *within* each site/condition.
#'
#' @return An object of class `env_space`, with slots
#'   `context`, `n_groups`, `variables`, `planned`,
#'   `realised`, and `within_group_sd`.
#' @export
define_env_locations <- function(locations,
                                 context = c("field", "controlled"),
                                 within_group_sd = 0) {

  context  <- match.arg(context)
  loc_df   <- as.data.frame(locations)

  # ----- identify / create group_id ---------------------------------
  if (!("group_id" %in% names(loc_df))) {
    # try common alternatives
    possible_id <- intersect(
      c("location", "site", "condition"),
      names(loc_df)
    )

    if (length(possible_id) > 0) {
      loc_df$group_id <- as.character(loc_df[[possible_id[1]]])
    } else {
      loc_df$group_id <- seq_len(nrow(loc_df))
    }
  }

  # Move group_id to first column for nicer printing
  id_col   <- loc_df$group_id
  env_cols <- setdiff(names(loc_df), "group_id")

  # keep only numeric env variables
  numeric_env <- env_cols[sapply(loc_df[env_cols], is.numeric)]
  if (length(numeric_env) == 0L) {
    stop("No numeric environmental variables found in `locations`.")
  }

  planned <- data.frame(
    group_id = id_col,
    loc_df[ , numeric_env, drop = FALSE]
  )

  variables <- numeric_env

  # ----- within-group SD handling -----------------------------------
  if (length(within_group_sd) == 1L) {
    within_group_sd <- rep(within_group_sd, length(variables))
    names(within_group_sd) <- variables
  } else {
    within_group_sd <- as.numeric(within_group_sd)
    if (is.null(names(within_group_sd))) {
      stop("If `within_group_sd` has length > 1, it must be a *named* vector ",
           "with names matching the environmental variables.")
    }
    if (!all(variables %in% names(within_group_sd))) {
      stop("Names of `within_group_sd` must include all variables: ",
           paste(variables, collapse = ", "))
    }
    within_group_sd <- within_group_sd[variables]
  }

  res <- list(
    context        = context,
    n_groups       = nrow(planned),
    variables      = variables,
    planned        = planned,
    realised       = planned,        # initially the same
    within_group_sd = within_group_sd
  )
  class(res) <- "env_space"
  res
}
