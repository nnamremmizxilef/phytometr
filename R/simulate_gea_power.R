#' Simulate GEA power analysis
#'
#' Performs power analysis for genotype-environment association (GEA) studies
#' by simulating genetic data, environmental gradients, and holobiont responses.
#'
#' @param n_individuals Numeric vector. Sample sizes to test (default: seq(20, 200, by = 20))
#' @param n_sites Integer. Number of sampling sites (optional)
#' @param individuals_per_site Integer. Individuals per site if n_sites specified
#' @param n_loci Integer. Number of genetic loci to simulate (default: 5000)
#' @param n_causal Integer. Number of causal loci (default: 10)
#' @param mutation_rate Numeric. Mutation rate for coalescent simulation (default: 10)
#' @param recombination_rate Numeric. Recombination rate (default: 1)
#' @param maf_threshold Numeric. Minor allele frequency filter threshold (default: 0.05)
#' @param pop_structure Character. Population structure type: "panmictic", "isolation_by_distance", 
#'   "discrete", or "admixture" (default: "panmictic")
#' @param n_populations Integer. Number of populations for structured models (default: 3)
#' @param migration_rate Numeric. Migration rate between populations (default: 0.01)
#' @param migration_symmetric Logical. Symmetric migration? (default: TRUE)
#' @param fst Numeric. Target Fst between populations (alternative to migration_rate)
#' @param bottleneck List. Bottleneck parameters: list(time = 0.5, reduction = 0.1)
#' @param expansion List. Expansion parameters: list(time = 0.3, growth_rate = 2)
#' @param env_type Character. Environment type: "continuous", "categorical", "multivariate" (default: "continuous")
#' @param env_range Numeric vector. Range for continuous environment (default: c(15, 25))
#' @param env_categories Character vector. Categories for categorical environment
#' @param env_n_variables Integer. Number of environmental variables if multivariate (default: 1)
#' @param env_correlation Numeric. Correlation between environmental variables (default: 0)
#' @param gradient_type Character. Gradient shape: "linear", "nonlinear", "patchy", "random" (default: "linear")
#' @param gradient_noise Numeric. Within-site environmental noise (default: 0.1)
#' @param effect_size Numeric. Strength of genotype-environment association (default: 0.3)
#' @param effect_size_distribution Character. Distribution of effect sizes: "equal", "exponential", "normal" (default: "equal")
#' @param selection_coefficient Numeric. Selection coefficient (alternative to effect_size)
#' @param dominance Numeric. Dominance coefficient: 0 (recessive), 0.5 (additive), 1 (dominant) (default: 0.5)
#' @param response_variable Character. Response type: "abundance", "diversity", "composition", "phenotype" (default: "abundance")
#' @param holobiont_heritability Numeric. Host genetic control of microbiome (default: 0.3)
#' @param env_effect_on_holobiont Numeric. Direct environmental effect on holobiont (default: 0.3)
#' @param gxe_interaction Logical. Include GxE interaction? (default: FALSE)
#' @param gxe_strength Numeric. Strength of GxE interaction if included (default: 0.2)
#' @param n_taxa Integer. Number of microbial taxa (default: 100)
#' @param microbiome_structure Character. Microbiome data type: "compositional", "presence_absence", "abundance" (default: "compositional")
#' @param baseline_diversity Numeric. Mean alpha diversity (default: 50)
#' @param dispersal_limitation Numeric. Spatial autocorrelation in microbiome (default: 0.1)
#' @param measurement_error Numeric. Technical measurement noise (default: 0.1)
#' @param biological_noise Numeric. Unexplained biological variation (default: 0.2)
#' @param missing_data_rate Numeric. Proportion of missing genotypes (default: 0)
#' @param n_replicates Integer. Number of simulation replicates per sample size (default: 100)
#' @param gea_method Character. GEA method(s): "rda", "lfmm", "baypass", "gemma", "all" (default: "rda")
#' @param correction_method Character. Multiple testing correction: "bonferroni", "fdr", "qvalue" (default: "fdr")
#' @param alpha Numeric. Significance threshold (default: 0.05)
#' @param n_latent_factors Integer. Number of latent factors for LFMM (default: NULL, auto-detect)
#' @param min_effect_to_detect Numeric. Minimum effect size of interest
#' @param target_power Numeric. Target statistical power (default: 0.8)
#' @param calculate_fdr Logical. Calculate false discovery rate? (default: TRUE)
#' @param parallel Logical. Use parallel processing? (default: TRUE)
#' @param n_cores Integer. Number of cores for parallel processing (default: detectCores() - 1)
#' @param seed Integer. Random seed for reproducibility
#' @param verbose Logical. Print progress messages? (default: TRUE)
#'
#' @return Object of class 'gea_power_analysis' containing:
#' \itemize{
#'   \item power: Mean statistical power for each sample size
#'   \item fdr: Mean false discovery rate
#'   \item n_detected: Mean number of loci detected
#'   \item parameters: List of simulation parameters used
#' }
#'
#' @examples
#' \dontrun{
#' # Basic power analysis
#' power_results <- simulate_gea_power(
#'   n_individuals = seq(50, 200, by = 25),
#'   n_causal = 10,
#'   effect_size = 0.3
#' )
#' plot(power_results)
#'
#' # Complex scenario with population structure
#' power_complex <- simulate_gea_power(
#'   n_individuals = seq(30, 150, by = 20),
#'   pop_structure = "isolation_by_distance",
#'   n_populations = 5,
#'   migration_rate = 0.02,
#'   gradient_type = "linear",
#'   response_variable = "composition",
#'   gxe_interaction = TRUE,
#'   n_replicates = 200
#' )
#' }
#'
#' @export
simulate_gea_power <- function(
  # Sample size parameters
  n_individuals = seq(20, 200, by = 20),
  n_sites = NULL,
  individuals_per_site = NULL,
  
  # Genetic parameters
  n_loci = 5000,
  n_causal = 10,
  mutation_rate = 10,
  recombination_rate = 1,
  maf_threshold = 0.05,
  
  # Population structure parameters
  pop_structure = c("panmictic", "isolation_by_distance", "discrete", "admixture"),
  n_populations = 3,
  migration_rate = 0.01,
  migration_symmetric = TRUE,
  fst = NULL,
  
  # Demographic history
  bottleneck = NULL,
  expansion = NULL,
  
  # Environmental parameters
  env_type = c("continuous", "categorical", "multivariate"),
  env_range = c(15, 25),
  env_categories = NULL,
  env_n_variables = 1,
  env_correlation = 0,
  
  # Environmental gradient shape
  gradient_type = c("linear", "nonlinear", "patchy", "random"),
  gradient_noise = 0.1,
  
  # Selection parameters
  effect_size = 0.3,
  effect_size_distribution = c("equal", "exponential", "normal"),
  selection_coefficient = NULL,
  dominance = 0.5,
  
  # Holobiont-specific parameters
  response_variable = c("abundance", "diversity", "composition", "phenotype"),
  holobiont_heritability = 0.3,
  env_effect_on_holobiont = 0.3,
  gxe_interaction = FALSE,
  gxe_strength = 0.2,
  
  # Microbiome parameters
  n_taxa = 100,
  microbiome_structure = c("compositional", "presence_absence", "abundance"),
  baseline_diversity = 50,
  dispersal_limitation = 0.1,
  
  # Noise parameters
  measurement_error = 0.1,
  biological_noise = 0.2,
  missing_data_rate = 0,
  
  # Statistical parameters
  n_replicates = 100,
  gea_method = c("rda", "lfmm", "baypass", "gemma", "all"),
  correction_method = c("bonferroni", "fdr", "qvalue"),
  alpha = 0.05,
  n_latent_factors = NULL,
  
  # Power calculation parameters
  min_effect_to_detect = NULL,
  target_power = 0.8,
  calculate_fdr = TRUE,
  
  # Computational parameters
  parallel = TRUE,
  n_cores = parallel::detectCores() - 1,
  seed = NULL,
  verbose = TRUE
) {
# Match arguments
  pop_structure <- match.arg(pop_structure)
  env_type <- match.arg(env_type)
  gradient_type <- match.arg(gradient_type)
  response_variable <- match.arg(response_variable)
  effect_size_distribution <- match.arg(effect_size_distribution)
  microbiome_structure <- match.arg(microbiome_structure)
  correction_method <- match.arg(correction_method)
  
  gea_method <- match.arg(gea_method)
  
  if (effect_size < 0) {
  stop("effect_size must be non-negative")
  }
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Validate parameters
  params <- validate_simulation_parameters(
    n_individuals = n_individuals,
    n_loci = n_loci,
    n_causal = n_causal,
    pop_structure = pop_structure,
    env_type = env_type,
    response_variable = response_variable,
    gea_method = gea_method,
    holobiont_heritability = holobiont_heritability,
    env_effect_on_holobiont = env_effect_on_holobiont,
    n_replicates = n_replicates
  )
  
  # Store all parameters for later use
  params$all <- list(
    n_loci = n_loci,
    n_causal = n_causal,
    mutation_rate = mutation_rate,
    recombination_rate = recombination_rate,
    maf_threshold = maf_threshold,
    pop_structure = pop_structure,
    n_populations = n_populations,
    migration_rate = migration_rate,
    migration_symmetric = migration_symmetric,
    bottleneck = bottleneck,
    expansion = expansion,
    env_type = env_type,
    env_range = env_range,
    env_categories = env_categories,
    env_n_variables = env_n_variables,
    env_correlation = env_correlation,
    gradient_type = gradient_type,
    gradient_noise = gradient_noise,
    effect_size = effect_size,
    effect_size_distribution = effect_size_distribution,
    dominance = dominance,
    response_variable = response_variable,
    holobiont_heritability = holobiont_heritability,
    env_effect_on_holobiont = env_effect_on_holobiont,
    gxe_interaction = gxe_interaction,
    gxe_strength = gxe_strength,
    n_taxa = n_taxa,
    microbiome_structure = microbiome_structure,
    baseline_diversity = baseline_diversity,
    measurement_error = measurement_error,
    biological_noise = biological_noise,
    missing_data_rate = missing_data_rate,
    gea_method = gea_method,
    correction_method = correction_method,
    alpha = alpha,
    n_latent_factors = n_latent_factors,
    target_power = target_power
  )
  
  # Setup parallel processing if requested
  if (parallel && n_replicates > 1) {
    if (verbose) message("Setting up parallel processing with ", n_cores, " cores")
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export necessary objects to cluster
    parallel::clusterExport(cl, varlist = c("params"), envir = environment())
    parallel::clusterEvalQ(cl, library(phytometr))
  }
  
  # Main simulation loop over sample sizes
  if (verbose) message("Starting power simulation...")
  
  power_results <- lapply(seq_along(n_individuals), function(i) {
    n <- n_individuals[i]
    
    if (verbose) {
      message(sprintf("Simulating sample size %d/%d: n = %d", 
                     i, length(n_individuals), n))
    }
    
    # Run replicates for this sample size
    replicates <- if (parallel && n_replicates > 1) {
      parallel::parLapply(cl, seq_len(n_replicates), function(rep) {
        run_single_replicate(n, params$all)
      })
    } else {
      lapply(seq_len(n_replicates), function(rep) {
        run_single_replicate(n, params$all)
      })
    }
    
    # Aggregate results across replicates
    aggregate_replicates(replicates, params$all)
  })
  
  # Format and return results
  results <- format_power_results(power_results, n_individuals, params$all)
  
  if (verbose) message("Power simulation complete!")
  
  return(results)
}
