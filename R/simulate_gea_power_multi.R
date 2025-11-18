#' Simulate GEA power analysis for multiple methods
#'
#' Runs \code{simulate_gea_power()} for several GEA methods and
#' returns all power curves in a single object.
#'
#' @param gea_methods Character vector of GEA methods:
#'   any of \code{c("rda", "lfmm", "baypass", "gemma", "all")}.
#'   If \code{"all"} is included, it is expanded to all implemented methods.
#' @inheritParams simulate_gea_power
#'
#' @return A named list of \code{"gea_power_analysis"} objects,
#'   with class \code{"gea_power_multi"} added.
#' @export
simulate_gea_power_multi <- function(
  gea_methods = c("rda", "lfmm"),
  ...
) {
  # Normalize method names
  gea_methods <- match.arg(
    gea_methods,
    choices = c("rda", "lfmm", "baypass", "gemma", "all"),
    several.ok = TRUE
  )

  # If "all" is requested, expand to all implemented methods
  if ("all" %in% gea_methods) {
    gea_methods <- c("rda", "lfmm", "baypass", "gemma")
    gea_methods <- unique(gea_methods)
  }

  # Run simulate_gea_power() for each method
  results_list <- lapply(gea_methods, function(m) {
    res <- simulate_gea_power(gea_method = m, ...)
    # make sure the method is stored in parameters for plotting
    params <- attr(res, "parameters")
    if (!is.null(params)) {
      params$gea_method <- m
      attr(res, "parameters") <- params
    }
    res
  })
  names(results_list) <- gea_methods

  class(results_list) <- c("gea_power_multi", "list")
  results_list
}
