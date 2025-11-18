#' Simulate holobiont responses to environmental conditions
#'
#' This is a high-level wrapper that takes:
#' \itemize{
#'   \item individual- or group-level environmental data,
#'   \item a specification of holobiont components and their
#'         environment-driven R² values (see \code{simulate_holobiont_components()}),
#'   \item an optional host model (clone vs natural population).
#' }
#'
#' For now, the host model is accepted but not used to add extra
#' host-driven variance (this can be extended in future phases).
#'
#' @param env_data A data frame or matrix with environmental values.
#'   Rows are individuals (or groups), columns are variables
#'   like "temperature", "moisture", etc.
#' @param components A named list describing holobiont components
#'   and their R² with environment (see
#'   \code{simulate_holobiont_components()}).
#' @param host_model Optional host model object created by
#'   \code{define_host_model()}. Currently included for API
#'   consistency; host-specific variance can be added in later
#'   development.
#'
#' @return A data frame with one row per row in \code{env_data}
#'   and one column per holobiont component (plus an \code{id}
#'   column).
#'
#' @examples
#' \dontrun{
#' env_data <- env_space$realised[, c("temperature", "moisture")]
#' components_spec <- list(
#'   root_bacteria = list(r2_env = c(temperature = 0.3)),
#'   root_fungi    = list(r2_env = c(temperature = 0.25, moisture = 0.10))
#' )
#' resp <- simulate_env_response(env_data, components_spec)
#' }
#'
#' @export
simulate_env_response <- function(env_data,
                                  components,
                                  host_model = NULL) {

  # Basic checks
  if (missing(env_data)) {
    stop("`env_data` must be provided and must contain environmental variables.")
  }
  if (missing(components)) {
    stop("`components` must be provided and must be a named list.")
  }

  env_data <- as.data.frame(env_data)

  # Use the lower-level engine to construct component responses
  resp <- simulate_holobiont_components(
    components = components,
    env_data   = env_data
  )

  # Placeholder: here we could incorporate host_model effects later
  # e.g. add random host variance for natural populations.

  return(resp)
}
