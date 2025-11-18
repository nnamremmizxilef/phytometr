#' Simulate power for detecting env–holobiont effects
#'
#' This function explores how design choices (number of groups/sites and
#' individuals per group) affect the power to detect environmental effects
#' on holobiont components.
#'
#' It uses three main ingredients:
#' \itemize{
#'   \item an \code{env_space} object from \code{define_env_space()},
#'   \item a \code{host_model} from \code{define_host_model()} (currently
#'         only used to distinguish clone vs natural; host variance is not
#'         yet simulated),
#'   \item a component specification list (the same structure you pass to
#'         \code{simulate_holobiont_components()}).
#' }
#'
#' For each combination of \code{n_groups_use} and \code{n_per_group}, the
#' function:
#' \enumerate{
#'   \item Simulates individual-level environmental values around the group
#'         means in \code{env_space}, using the within-group SDs stored there.
#'   \item Simulates holobiont component values from the specified R²
#'         relationships using \code{simulate_holobiont_components()}.
#'   \item Fits simple linear models of each component on all environmental
#'         variables, and records whether at least one environmental
#'         coefficient is significant at \code{alpha}.
#'   \item Estimates power for each component as the fraction of replicates
#'         with at least one significant environmental effect.
#' }
#'
#' @param env_space An object created by \code{define_env_space()}.
#' @param host_model Host model object from \code{define_host_model()}.
#'        (Currently included for API consistency; host variance is not yet
#'        simulated.)
#' @param components A named list describing holobiont components and their
#'        environment-driven R², exactly as used in
#'        \code{simulate_holobiont_components()}.
#' @param n_groups_use Integer vector with candidate numbers of groups/sites
#'        to sample from \code{env_space}.
#' @param n_per_group Integer vector with candidate numbers of individuals
#'        per group.
#' @param n_replicates Number of simulation replicates per design.
#' @param alpha Significance level for detecting any environmental effect.
#'
#' @return A data.frame with one row per design and columns:
#'   \itemize{
#'     \item \code{n_groups}, \code{n_per_group}, \code{N} – design size,
#'     \item \code{power_env} – mean power across components,
#'     \item one column per component with its individual power.
#'   }
#'   The object is given class \code{"holobiont_power"}.
#'
#' @export
simulate_holobiont_power <- function(env_space,
                                     host_model,
                                     components,
                                     n_groups_use,
                                     n_per_group,
                                     n_replicates = 100,
                                     alpha = 0.05) {

  # ---- Basic checks -----------------------------------------------------

  if (!inherits(env_space, "env_space")) {
    stop("`env_space` must be an 'env_space' object from define_env_space().")
  }
  if (!is.list(components) || is.null(names(components))) {
    stop("`components` must be a named list (the same structure used by simulate_holobiont_components()).")
  }
  if (missing(n_groups_use) || missing(n_per_group)) {
    stop("Please provide both `n_groups_use` and `n_per_group`.")
  }

  env_means <- env_space$realised
  env_vars  <- env_space$variables
  within_sd <- env_space$within_group_sd[env_vars]

  # internal helper: simulate individual-level env around group means
  simulate_ind_env <- function(n_groups, n_per_group) {

    if (n_groups > nrow(env_means)) {
      stop("Requested n_groups exceeds the number of groups in env_space$realised.")
    }

    chosen_ids <- sort(sample(env_means$group_id, n_groups))
    pieces <- lapply(chosen_ids, function(gid) {
      mu_row <- env_means[env_means$group_id == gid, env_vars, drop = FALSE]
      mu     <- as.numeric(mu_row[1, ])

      env_mat <- vapply(
        seq_along(env_vars),
        function(j) stats::rnorm(n_per_group, mean = mu[j], sd = within_sd[j]),
        FUN.VALUE = numeric(n_per_group)
      )

      df <- as.data.frame(env_mat)
      names(df) <- env_vars
      df$group_id <- gid
      df
    })

    out <- do.call(rbind, pieces)
    out <- out[, c("group_id", env_vars), drop = FALSE]
    rownames(out) <- NULL
    out
  }

  designs <- expand.grid(
    n_groups   = n_groups_use,
    n_per_group = n_per_group,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  comp_names <- names(components)
  results_list <- vector("list", nrow(designs))

  # ---- Loop over designs -----------------------------------------------

  for (i in seq_len(nrow(designs))) {

    ng  <- designs$n_groups[i]
    npg <- designs$n_per_group[i]
    N   <- ng * npg

    message(sprintf("[Design %d] n_groups = %d, n_per_group = %d (N = %d)",
                    i, ng, npg, N))

    # record successes per component
    succ <- matrix(0, nrow = n_replicates, ncol = length(comp_names))
    colnames(succ) <- comp_names

    for (rep in seq_len(n_replicates)) {

      # 1) simulate individual envs
      env_ind <- simulate_ind_env(ng, npg)

      # 2) simulate holobiont components
      comp_sim <- simulate_holobiont_components(
        components = components,
        env_data   = env_ind[, env_vars, drop = FALSE]
      )

      # merge env + components (drop comp_sim$id; same ordering)
      dat <- cbind(
        env_ind,
        comp_sim[, setdiff(names(comp_sim), "id"), drop = FALSE]
      )

      # 3) for each component, fit lm and see if any env effect is sig.
      for (j in seq_along(comp_names)) {
        cname <- comp_names[j]

        f <- stats::as.formula(
          paste(cname, "~", paste(env_vars, collapse = " + "))
        )

        fit <- stats::lm(f, data = dat)

        # drop intercept
        pvals <- summary(fit)$coefficients[-1, 4]

        succ[rep, j] <- as.integer(any(pvals < alpha))
      }
    }

    power_per_comp <- colMeans(succ)
    power_env      <- mean(power_per_comp)

    results_list[[i]] <- data.frame(
      n_groups    = ng,
      n_per_group = npg,
      N           = N,
      power_env   = power_env,
      t(power_per_comp),
      row.names   = NULL,
      check.names = FALSE
    )
  }

  out <- do.call(rbind, results_list)
  class(out) <- c("holobiont_power", class(out))
  attr(out, "components")      <- comp_names
  attr(out, "n_replicates")    <- n_replicates
  attr(out, "alpha")           <- alpha
  attr(out, "env_variables")   <- env_vars
  attr(out, "host_model_type") <- host_model$type

  out
}
