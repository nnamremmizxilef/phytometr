# ============================================================
# diagnose_holobiont_design()
# Full multi-panel diagnostic engine
# ============================================================

#' Diagnose a holobiont experimental design
#'
#' @description
#' Provides a full statistical diagnostic report for a given design:
#' environmental structure, predicted holobiont behaviour, variance
#' decomposition, PCA structure, and correlation matrices.
#'
#' @param env_space Output of define_env_space() or define_env_locations()
#' @param components List of holobiont components with r2_env specification
#' @param host_model Output of define_host_model()
#'
#' @return An object of class "phytomet_diagnostics" containing:
#' \itemize{
#'   \item summary — key warnings and metrics
#'   \item env_cor — correlation matrix among environmental variables
#'   \item comp_cor — correlation matrix among holobiont components
#'   \item env_comp_cor — correlation env → holobiont
#'   \item pca — PCA list with scores & loadings
#'   \item plots — a list of ggplot objects
#'
#' }
#'
#' @export
diagnose_holobiont_design <- function(env_space, components, host_model) {

  # ============================================================
  # SAFETY CHECKS
  # ============================================================

  if (!inherits(env_space, "env_space") &&
      !inherits(env_space, "env_locations") &&
      !inherits(env_space, "env_grid")) {
    stop("env_space must be created by define_env_space(), define_env_grid(), or define_env_locations().")
  }

  if (!is.list(components)) {
    stop("components must be a named list of holobiont components with r2_env values.")
  }

  # ============================================================
  # 1. COMPUTE ENVIRONMENTAL DIAGNOSTICS
  # ============================================================

  env_vals <- env_space$realised
  if ("group_id" %in% colnames(env_vals))
    env_vals <- env_vals[, !(colnames(env_vals) == "group_id")]

  env_cor <- stats::cor(env_vals)
  env_cor_df <- as.data.frame(env_cor)

  # Issue warnings for multicollinearity
  multicol_warnings <- names(which(abs(env_cor) > 0.8 & env_cor != 1))

  # ============================================================
  # 2. SIMULATE A REPRESENTATIVE HOLOBIONT DATASET
  # ============================================================

  # Simulate "environment per individual"
  env_indiv <- simulate_individual_env(env_space, n_per_group = 20)

  holobiont_sim <- simulate_holobiont_components(
    components = components,
    env_data   = env_indiv
  )

  comp_vals <- holobiont_sim[, names(holobiont_sim) != "id"]
  comp_cor <- stats::cor(comp_vals)

  # ============================================================
  # 3. ENV → HOLOBIONT CORRELATION MATRIX
  # ============================================================

  env_indiv_only <- env_indiv[, !(colnames(env_indiv) == "id")]
  env_comp_cor <- stats::cor(env_indiv_only[, -1], # remove group id
                      comp_vals[, -1])      # remove group id

  # ============================================================
  # 4. PCA DIAGNOSTICS
  # ============================================================

  pca <- stats::prcomp(comp_vals[, -1], scale. = TRUE)  # remove group id

  # ============================================================
  # 5. VARIANCE DECOMPOSITION EXPECTATIONS
  # ============================================================

  var_table <- lapply(names(components), function(cmp) {
    r2_env <- components[[cmp]]$r2_env
    r2_total <- sum(r2_env)
    data.frame(
      component = cmp,
      r2_env = r2_total,
      r2_host = host_model$var_host,
      r2_residual = 1 - r2_total - host_model$var_host
    )
  })

  var_table <- do.call(rbind, var_table)

  # ============================================================
  # 6. BUILD SUMMARY OF WARNINGS
  # ============================================================

  summary_list <- list()

  if (length(multicol_warnings) > 0) {
    summary_list$environment_multicollinearity <-
      paste("Strong environmental correlations detected among:",
            paste(multicol_warnings, collapse = ", "))
  }

  high_r2 <- which(var_table$r2_env > 0.8)
  if (length(high_r2) > 0) {
    summary_list$high_env_r2 <-
      paste("Components with extremely high r2_env:",
            paste(var_table$component[high_r2], collapse = ", "))
  }

  if (host_model$type == "natural" && host_model$var_host < 0.05) {
    summary_list$low_host_variance <-
      "Natural population selected but host variance is extremely low."
  }

  if (length(summary_list) == 0) {
    summary_list$all_good <- "No critical issues detected."
  }

  # ============================================================
  # 7. PRODUCE PLOTS USING YOUR EXISTING FUNCTIONS
  # ============================================================

  plot_env <- plot_env_space(env_space)

  plot_components <- plot_holobiont_components(
    components_df = holobiont_sim,
    env           = env_indiv[, c("temperature", "moisture")],
    color_by      = env_indiv$group_id
  )

  # ============================================================
  # RETURN OBJECT
  # ============================================================

  out <- list(
    summary        = summary_list,
    env_cor        = env_cor_df,
    comp_cor       = comp_cor,
    env_comp_cor   = env_comp_cor,
    pca            = pca,
    var_table      = var_table,
    plots          = list(
      env_space      = plot_env,
      pca            = plot_components$PCA,
      comp_corr      = plot_components$CorMatrix,
      env_effects    = plot_components$EnvEffects
    )
  )

  class(out) <- "phytomet_diagnostics"
  return(out)
}
