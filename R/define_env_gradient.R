#' Define a potential environmental gradient or set of conditions
#'
#' This function defines a set of *potential* sites (field) or
#' conditions (controlled experiments) along a single environmental
#' axis (e.g. temperature, moisture), without deciding how many
#' individuals will be sampled. It is meant as a first step to
#' describe what environmental range and structure are realistically
#' achievable.
#'
#' You specify:
#' \itemize{
#'   \item the feasible environmental range (e.g. 15--25 °C in the field,
#'         10--35 °C in growth chambers),
#'   \item the number of potential sites/conditions along that range,
#'   \item the shape of the gradient (linear, random, patchy),
#'   \item how precisely you can hit the targeted group means
#'         (\code{target_sd}),
#'   \item how much variation you expect within each group
#'         (\code{within_group_sd}).
#' }
#'
#' The output describes the *environment space* for a single
#' variable that later power functions can use to decide how many
#' sites/conditions and replicates per group are needed.
#'
#' @param context Character, either \code{"field"} or \code{"controlled"}.
#'   This is mainly descriptive, but typically:
#'   \itemize{
#'     \item \code{"field"} for natural gradients (limited range, less control),
#'     \item \code{"controlled"} for chambers/greenhouses (wider range, more control).
#'   }
#' @param env_range Numeric length-2 vector \code{c(min, max)} giving
#'   the feasible range of the environmental variable.
#' @param n_groups Integer. Number of potential sites/conditions along
#'   the gradient. This does \emph{not} fix the number of individuals;
#'   it only defines possible group mean values.
#' @param shape Character. Shape of the gradient of planned group means:
#'   \itemize{
#'     \item \code{"linear"}: evenly spaced from \code{env_range[1]} to
#'           \code{env_range[2]}.
#'     \item \code{"random"}: group means randomly scattered within the
#'           range, then sorted.
#'     \item \code{"patchy"}: group means clustered toward the extremes
#'           of the range (approximate patchiness).
#'   }
#' @param target_sd Numeric >= 0. Standard deviation of the deviation
#'   of realised group means from planned values. Larger values mean
#'   you cannot exactly hit the intended condition (e.g. inter-annual
#'   variability in the field). Use small values (e.g. 0.1) for
#'   controlled experiments.
#' @param within_group_sd Numeric >= 0. Expected standard deviation of
#'   environmental values within each site/condition (micro-variation).
#'   This is stored as metadata for later individual-level simulation.
#' @param seed Optional integer. If provided, used to set the random
#'   seed for reproducibility.
#'
#' @return A data.frame with one row per potential site/condition, and
#'   columns:
#'   \itemize{
#'     \item \code{group_id}: integer group index.
#'     \item \code{context}: field or controlled.
#'     \item \code{planned_env}: planned group mean value.
#'     \item \code{realised_env}: realised group mean value after
#'           adding \code{target_sd} noise and clipping to \code{env_range}.
#'     \item \code{within_group_sd}: expected within-group SD.
#'   }
#'
#' @examples
#' # Field gradient: 15-25 °C, 6 sites, moderate control
#' field_grad <- define_env_gradient(
#'   context         = "field",
#'   env_range       = c(15, 25),
#'   n_groups        = 6,
#'   shape           = "linear",
#'   target_sd       = 1,    # years can be off by ~1 °C
#'   within_group_sd = 0.5   # micro-variation within site
#' )
#' field_grad
#'
#' # Controlled experiment: 10-35 °C, 8 chambers, very precise
#' ctrl_grad <- define_env_gradient(
#'   context         = "controlled",
#'   env_range       = c(10, 35),
#'   n_groups        = 8,
#'   shape           = "linear",
#'   target_sd       = 0.1,  # you can hit target ~exactly
#'   within_group_sd = 0.2
#' )
#' ctrl_grad
#'
#' @export
define_env_gradient <- function(
  context         = c("field", "controlled"),
  env_range       = c(15, 25),
  n_groups        = 5L,
  shape           = c("linear", "random", "patchy"),
  target_sd       = 0,
  within_group_sd = 0,
  seed            = NULL
) {
  # Match arguments and basic checks
  context <- match.arg(context)
  shape   <- match.arg(shape)

  if (!is.numeric(env_range) || length(env_range) != 2L) {
    stop("'env_range' must be a numeric vector of length 2: c(min, max).")
  }
  env_min <- env_range[1]
  env_max <- env_range[2]
  if (env_min >= env_max) {
    stop("'env_range[1]' must be < 'env_range[2]'.")
  }

  if (!is.numeric(n_groups) || length(n_groups) != 1L || n_groups < 1) {
    stop("'n_groups' must be a positive integer.")
  }
  n_groups <- as.integer(n_groups)

  if (!is.numeric(target_sd) || length(target_sd) != 1L || target_sd < 0) {
    stop("'target_sd' must be a single number >= 0.")
  }
  if (!is.numeric(within_group_sd) || length(within_group_sd) != 1L || within_group_sd < 0) {
    stop("'within_group_sd' must be a single number >= 0.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Define planned group means
  if (shape == "linear") {
    planned_env <- seq(env_min, env_max, length.out = n_groups)
  } else if (shape == "random") {
    planned_env <- sort(stats::runif(n_groups, min = env_min, max = env_max))
  } else if (shape == "patchy") {
    # Rough "patchy" behaviour: half of groups near lower end, half near upper end
    n_low  <- floor(n_groups / 2)
    n_high <- n_groups - n_low
    range_quarter <- (env_max - env_min) / 4

    low_vals  <- stats::runif(n_low,
                              min = env_min,
                              max = env_min + range_quarter)
    high_vals <- stats::runif(n_high,
                              min = env_max - range_quarter,
                              max = env_max)
    planned_env <- sort(c(low_vals, high_vals))
  }

  # Realised group means (imperfect control)
  if (target_sd > 0) {
    realised_env <- planned_env + stats::rnorm(n_groups, mean = 0, sd = target_sd)
    realised_env <- pmax(pmin(realised_env, env_max), env_min)  # clip to range
  } else {
    realised_env <- planned_env
  }

  out <- data.frame(
    group_id        = seq_len(n_groups),
    context         = context,
    planned_env     = planned_env,
    realised_env    = realised_env,
    within_group_sd = within_group_sd
  )

  return(out)
}
