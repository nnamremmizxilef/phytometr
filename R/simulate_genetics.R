#' Simulate host genetics with population structure
#'
#' @param holobiont A holobiont_data object
#' @param params List of simulation parameters
#' @return Updated holobiont_data object with host_genotypes
#' @keywords internal
simulate_host_genetics <- function(holobiont, params) {
  
  if (!requireNamespace("coala", quietly = TRUE)) {
    stop("Package 'coala' is required for genetic simulations. Please install it with: install.packages('coala')")
  }
  
  # Check if coala simulate function exists
  if (!exists("simulate", where = asNamespace("coala"), mode = "function")) {
    stop("coala package is installed but 'simulate' function not found. Try reinstalling coala.")
  }
  
  n <- length(holobiont$host_id)
  
  # Build coala model based on population structure
  if (params$pop_structure == "panmictic") {
    model <- coala::coal_model(
      sample_size = n,
      loci_number = params$n_loci
    ) +
    coala::feat_mutation(rate = params$mutation_rate) +
    coala::feat_recombination(rate = params$recombination_rate)
    
  } else if (params$pop_structure == "isolation_by_distance") {
    n_demes <- params$n_populations
    samples_per_deme <- rep(floor(n / n_demes), n_demes)
    samples_per_deme[1] <- samples_per_deme[1] + (n - sum(samples_per_deme))
    
    model <- coala::coal_model(
      sample_size = samples_per_deme,
      loci_number = params$n_loci
    ) +
    coala::feat_mutation(rate = params$mutation_rate) +
    coala::feat_migration(
      rate = params$migration_rate, 
      symmetric = params$migration_symmetric
    ) +
    coala::feat_recombination(rate = params$recombination_rate)
    
    # Store population assignments
    holobiont$population <- rep(seq_len(n_demes), times = samples_per_deme)
    
  } else if (params$pop_structure == "discrete") {
    n_pops <- params$n_populations
    samples_per_pop <- rep(floor(n / n_pops), n_pops)
    samples_per_pop[1] <- samples_per_pop[1] + (n - sum(samples_per_pop))
    
    # Lower migration for discrete populations
    mig_rate <- params$migration_rate * 0.1
    
    model <- coala::coal_model(
      sample_size = samples_per_pop,
      loci_number = params$n_loci
    ) +
    coala::feat_mutation(rate = params$mutation_rate) +
    coala::feat_migration(rate = mig_rate) +
    coala::feat_recombination(rate = params$recombination_rate)
    
    holobiont$population <- rep(seq_len(n_pops), times = samples_per_pop)
  }
  
  # Add demographic events if specified
  if (!is.null(params$bottleneck)) {
    model <- model + coala::feat_size_change(
      time = params$bottleneck$time,
      population = "all",
      size = params$bottleneck$reduction
    )
  }
  
  if (!is.null(params$expansion)) {
    model <- model + coala::feat_growth(
      rate = params$expansion$growth_rate,
      time = params$expansion$time
    )
  }
  
  # Simulate genetic data
  sim_data <- simulate(model)
  snp_matrix <- convert_to_snp_matrix(sim_data)
  
  # Filter by MAF
  if (params$maf_threshold > 0 && ncol(snp_matrix) > 0) {
    maf <- colMeans(snp_matrix, na.rm = TRUE) / 2
    maf <- pmin(maf, 1 - maf)
    snp_matrix <- snp_matrix[, maf >= params$maf_threshold, drop = FALSE]
  }
  
  # Add missing data if specified
  if (params$missing_data_rate > 0 && length(snp_matrix) > 0) {
    n_missing <- floor(length(snp_matrix) * params$missing_data_rate)
    if (n_missing > 0) {
      missing_idx <- sample(length(snp_matrix), n_missing)
      snp_matrix[missing_idx] <- NA
    }
  }
  
  # Check if we have enough SNPs
  if (ncol(snp_matrix) < params$n_causal) {
    warning(sprintf("Only %d SNPs generated after MAF filtering (need %d causal). Consider lowering maf_threshold or increasing mutation_rate.",
                   ncol(snp_matrix), params$n_causal))
  }
  
  holobiont$host_genotypes <- snp_matrix
  return(holobiont)
}

#' Convert coala simulation output to SNP matrix
#'
#' @param coala_sim Output from coala::simulate()
#' @return SNP matrix (individuals x loci)
#' @keywords internal
convert_to_snp_matrix <- function(coala_sim) {
  
  # Extract segregating sites from all loci
  snp_list <- lapply(coala_sim$seg_sites, function(locus) {
    if (nrow(locus) > 0) {
      # Remove position column and transpose
      snp_data <- as.matrix(locus[, -1, drop = FALSE])
      t(snp_data)
    } else {
      NULL
    }
  })
  
  # Combine loci
  valid_snps <- snp_list[!sapply(snp_list, is.null)]
  
  if (length(valid_snps) == 0) {
    # No SNPs generated - return empty matrix with correct dimensions
    return(matrix(0, nrow = nrow(coala_sim$seg_sites[[1]]), ncol = 0))
  }
  
  snp_matrix <- do.call(cbind, valid_snps)
  
  # Ensure numeric matrix
  storage.mode(snp_matrix) <- "numeric"
  
  return(snp_matrix)
}