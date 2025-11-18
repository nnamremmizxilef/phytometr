#' Power analysis for holobiont responses to environment
#'
#' This function performs Monte Carlo power analysis for detecting
#' environmental effects on holobiont responses (whole holobiont or
#' components), given:
#' \itemize{
#'   \item a multivariate environmental space (\code{env_space}),
#'   \item a host model (clone vs natural; \code{host_model}),
#'   \item one or more holobiont components with literature-based
#'         environmental effect sizes (\code{components}),
#'   \item a sampling design: number of groups (sites/conditions)
#'         and individuals per group.
#' }
#'
#' For each design, the function:
#' \enumerate{
#'   \item simulates individual-level environment from \code{env_space},
#'   \item constructs a response for each component:
#'         \eqn{Y = Y_\mathrm{env} + Y_\mathrm{host} + \epsilon},
#'   \item fits a linear model \code{Y ~ env_1 + ... + env_p},
#'   \item tests the environmental block using an F-test,
#'   \item estimates power as the proportion of simulations with
#'         p-value < \code{alpha}.
#' }
#'
#' Host structure is taken into account via \code{host_model$var_host}:
#' \itemize{
#'   \item for \code{type = "clone"}, \code{var_host = 0} (no host variance),
#'   \item for \code{type = "natural"}, \code{var_host > 0} adds extra
#'         between-individual variation, reducing power for the same N.
#' }
#'
#' @param env_space Object of class \code{"env_space"} created by
#'   \code{define_env_space()}.
#' @param host_model Object of class \code{"host_model"} created by
#'   \code{define_host_model()}.
#' @param components Named list describing each holobiont component.
#'   Each element should be a list with at least:
#'   \itemize{
#'     \item \code{r2_env}: numeric vector between 0 and 1, either:
#'       \itemize{
#'         \item length 1: total R² for this component explained by
#'               \emph{all} environmental variables together, or
#'         \item length = number of environmental variables: per-variable
#'               R² contributions, whose sum is the total environmental
#'               R² for this component.
#'       }
#'   }
#'   The sum of \code{r2_env} plus \code{host_model$var_host} must be
#'   <= 1. The remaining variance is residual.
#' @param n_groups_use Integer vector. Candidate numbers of groups
#'   (sites/conditions) to use from \code{env_space}. For each value,
#'   groups are chosen evenly across the gradient.
#' @param n_per_group Integer vector. Candidate numbers of individuals
#'   per group. The function evaluates the Cartesian product of
#'   \code{n_groups_use} and \code{n_per_group}.
#' @param n_replicates Integer. Number of simulation replicates per
#'   design. Larger values give more precise power estimates but are
#'   slower.
#' @param alpha Numeric between 0 and 1. Significance threshold.
#' @param seed Optional integer. Base seed. Different designs and
#'   replicates will use offsets of this seed to keep simulations
#'   reproducible but distinct.
#' @param verbose Logical. If TRUE, print progress messages.
#'
#' @return A data.frame with one row per design and component, and
#'   columns:
#'   \itemize{
#'     \item \code{component}: name of the holobiont component.
#'     \item \code{n_groups}: number of groups used.
#'     \item \code{n_per_group}: individuals per group.
#'     \item \code{N_total}: total sample size.
#'     \item \code{env_r2_total}: total env R² for this component
#'           (sum of \code{r2_env}).
#'     \item \code{var_host}: host variance fraction.
#'     \item \code{power_env}: estimated power to detect an environmental
#'           effect (joint F-test of all env variables).
#'   }
#'
#' @examples
#' \dontrun{
#' # 1) Define environment
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
#' # 2) Host model: natural population with 20% host-driven variance
#' host_nat <- define_host_model(
#'   type              = "natural",
#'   pop_structure     = "panmictic",
#'   mutation_rate     = 1,
#'   recombination_rate = 0,
#'   var_host          = 0.2
#' )
#'
#' # 3) Holobiont components with env effect sizes
#' components <- list(
#'   root_bacteria = list(
#'     r2_env = c(temperature = 0.25, moisture = 0.05)  # 30% total
#'   ),
#'   root_fungi = list(
#'     r2_env = c(temperature = 0.15, moisture = 0.05)  # 20% total
#'   )
#' )
#'
#' # 4) Power analysis across designs
#' power_res <- simulate_holobiont_power(
#'   env_space    = env_space,
#'   host_model   = host_nat,
#'   components   = components,
#'   n_groups_use = c(3, 6),
#'   n_per_group  = c(10, 20),
#'   n_replicates = 200,
#'   alpha        = 0.05,
#'   seed         = 42,
#'   verbose      = TRUE
#' )
#'
#' power_res
#' }
#'
#' @export
simulate_holobiont_power <- function(
  env_space,
  host_model,
  components,
  n_groups_use,
  n_per_group,
  n_replicates = 200,
  alpha = 0.05,
  seed = NULL,
  verbose = TRUE
) {
  if (!inherits(env_space, "env_space")) {
    stop("'env_space' must be an object created by define_env_space().")
  }
  if (!inherits(host_model, "host_model")) {
    stop("'host_model' must be an object created by define_host_model().")
  }
  if (!is.list(components) || is.null(names(components))) {
    stop("'components' must be a *named* list.")
  }

  var_names <- env_space$variables
  p_env     <- length(var_names)

  # Check designs
  n_groups_use <- as.integer(n_groups_use)
  n_per_group  <- as.integer(n_per_group)

  if (any(n_groups_use < 1)) {
    stop("All values in 'n_groups_use' must be >= 1.")
  }
  if (any(n_per_group < 1)) {
    stop("All values in 'n_per_group' must be >= 1.")
  }
  if (any(n_groups_use > env_space$n_groups)) {
    stop("Cannot use more groups than available in env_space$n_groups.")
  }

  # Prepare design grid
  design_grid <- expand.grid(
    n_groups    = n_groups_use,
    n_per_group = n_per_group
  )
  design_grid$design_id <- seq_len(nrow(design_grid))

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Host variance fraction
  var_host <- host_model$var_host

  # Result container
  res_list <- list()

  # Loop over designs
  for (d in seq_len(nrow(design_grid))) {
    n_groups_d    <- design_grid$n_groups[d]
    n_per_group_d <- design_grid$n_per_group[d]
    N_total_d     <- n_groups_d * n_per_group_d

    # Choose groups roughly evenly across gradient
    all_groups <- env_space$planned$group_id
    idx <- round(seq(1, length(all_groups), length.out = n_groups_d))
    groups_d <- all_groups[idx]

    if (verbose) {
      message(
        "[Design ", d, "] n_groups = ", n_groups_d,
        ", n_per_group = ", n_per_group_d,
        " (N = ", N_total_d, ")"
      )
    }

    # Pre-allocate power counters for each component
    comp_names <- names(components)
    n_comp     <- length(components)
    hits_env   <- numeric(n_comp)
    total_env  <- rep(n_replicates, n_comp)

    # Loop over replicates
    for (r in seq_len(n_replicates)) {
      # Use different seed per replicate & design (if base seed was given)
      if (!is.null(seed)) {
        set.seed(seed + 1000L * d + r)
      }

      # 1) Simulate individual environment
      ind_env <- simulate_individual_env(
        env_space  = env_space,
        n_per_group = n_per_group_d,
        groups      = groups_d
      )
      env_mat <- as.matrix(ind_env[, var_names, drop = FALSE])

      # 2) For each component, simulate response and test env effect
      for (c_idx in seq_along(components)) {
        comp_name <- comp_names[c_idx]
        comp_spec <- components[[c_idx]]

        if (is.null(comp_spec$r2_env)) {
          stop("Component '", comp_name, "' must have an 'r2_env' element.")
        }

        r2_env_vec <- comp_spec$r2_env
        # Interpret r2_env: either single total or per-variable vector
        if (length(r2_env_vec) == 1L) {
          total_r2_env <- as.numeric(r2_env_vec)
          # Split equally across variables internally for env_effect
          r2_env_per_var <- rep(1 / p_env, p_env)
        } else {
          if (length(r2_env_vec) != p_env) {
            stop(
              "Component '", comp_name,
              "': length(r2_env) must be 1 or equal to number of env variables (", p_env, ")."
            )
          }
          # If named, ensure order matches var_names
          if (!is.null(names(r2_env_vec))) {
            r2_env_vec <- r2_env_vec[var_names]
          }
          total_r2_env <- sum(r2_env_vec)
          if (total_r2_env <= 0) {
            stop("Component '", comp_name, "': sum(r2_env) must be > 0.")
          }
          # Internal relative contributions (sum to 1)
          r2_env_per_var <- r2_env_vec / sum(r2_env_vec)
        }

        if (total_r2_env + var_host > 1) {
          stop(
            "For component '", comp_name, "': total_r2_env + var_host > 1. ",
            "Please reduce environmental R² or host variance."
          )
        }

        # 2a) Simulate *pure* env-driven signal (no residual)
        env_effect <- simulate_env_response(
          env       = env_mat,
          r2_env    = r2_env_per_var,
          total_var = total_r2_env
        )

        # 2b) Simulate host effect (clone vs natural)
        if (var_host > 0) {
          host_effect <- stats::rnorm(N_total_d, mean = 0, sd = sqrt(var_host))
        } else {
          host_effect <- rep(0, N_total_d)
        }

        # 2c) Residual noise to reach total variance 1
        var_resid <- 1 - total_r2_env - var_host
        if (var_resid < 0) var_resid <- 0
        resid_effect <- if (var_resid > 0) {
          stats::rnorm(N_total_d, mean = 0, sd = sqrt(var_resid))
        } else {
          rep(0, N_total_d)
        }

        response <- env_effect + host_effect + resid_effect

        # 3) Fit LM and test environmental block
        dat_lm <- data.frame(
          response = as.numeric(response),
          env_mat
        )

        fit_full <- stats::lm(response ~ ., data = dat_lm)
        fit_null <- stats::lm(response ~ 1, data = dat_lm)

        an <- stats::anova(fit_null, fit_full)
        p_env <- an$`Pr(>F)`[2L]

        if (!is.na(p_env) && p_env < alpha) {
          hits_env[c_idx] <- hits_env[c_idx] + 1
        }
      } # end component loop
    }   # end replicate loop

    # Summarise power for this design
    for (c_idx in seq_along(components)) {
      comp_name   <- comp_names[c_idx]
      comp_spec   <- components[[c_idx]]
      r2_env_vec  <- comp_spec$r2_env
      total_r2_env <- if (length(r2_env_vec) == 1L) {
        as.numeric(r2_env_vec)
      } else {
        sum(r2_env_vec)
      }

      res_list[[length(res_list) + 1L]] <- data.frame(
        component   = comp_name,
        n_groups    = n_groups_d,
        n_per_group = n_per_group_d,
        N_total     = N_total_d,
        env_r2_total = total_r2_env,
        var_host    = var_host,
        power_env   = hits_env[c_idx] / total_env[c_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  return(res)
}
