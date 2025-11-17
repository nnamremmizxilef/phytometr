#' Run GEA analysis
#'
#' @param holobiont A holobiont_data object
#' @param params Simulation parameters
#' @return List with p-values for each SNP
#' @keywords internal
run_gea_analysis <- function(holobiont, params) {

  snps <- holobiont$host_genotypes
  env <- holobiont$environment$main

  # Get response variable
  response <- extract_response_variable(holobiont, params)

  # Handle categorical environment
  if (is.character(env) || is.factor(env)) {
    env <- as.numeric(factor(env))
  }

  # Run appropriate GEA method
  if (params$gea_method == "rda" || params$gea_method == "all") {
    pvals <- run_rda(snps, env, response)
  } else if (params$gea_method == "lfmm") {
    pvals <- run_lfmm(snps, env, params)
  } else {
    # Placeholder for other methods
    pvals <- run_simple_correlation(snps, env, response)
  }

  return(list(pvals = pvals, method = params$gea_method))
}

#' Extract response variable from holobiont object
#'
#' @param holobiont Holobiont data object
#' @param params Simulation parameters
#' @return Numeric vector of response values
#' @keywords internal
extract_response_variable <- function(holobiont, params) {

  if (params$response_variable == "abundance") {
    response <- holobiont$microbiome$target_taxon
  } else if (params$response_variable == "diversity") {
    response <- holobiont$microbiome$alpha_diversity
  } else if (params$response_variable == "composition") {
    response <- holobiont$microbiome$pc1
  } else if (params$response_variable == "phenotype") {
    response <- holobiont$phenotypes$trait
  }

  return(response)
}

#' Run RDA-based GEA
#'
#' @param snps SNP matrix
#' @param env Environmental variable
#' @param response Response variable (optional, used for validation)
#' @return Numeric vector of p-values
#' @keywords internal
run_rda <- function(snps, env, response = NULL) {

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for RDA analysis")
  }

  # Remove SNPs with no variation or too much missing data
  snp_var <- apply(snps, 2, var, na.rm = TRUE)
  snps_filtered <- snps[, snp_var > 0 & !is.na(snp_var), drop = FALSE]

  if (ncol(snps_filtered) == 0) {
    return(rep(1, ncol(snps)))
  }

  # Run RDA
  tryCatch({
    rda_result <- vegan::rda(snps_filtered ~ env)

    # Extract loadings for constrained axis
    loadings <- scores(rda_result, choices = 1, display = "species",
                       scaling = "species")

    # Convert loadings to z-scores and p-values
    z_scores <- loadings / sd(loadings, na.rm = TRUE)
    pvals_filtered <- 2 * pnorm(-abs(z_scores))

    # Map back to original SNP positions
    pvals <- rep(1, ncol(snps))
    pvals[snp_var > 0 & !is.na(snp_var)] <- pvals_filtered

    return(pvals)

  }, error = function(e) {
    warning("RDA failed: ", e$message)
    return(rep(1, ncol(snps)))
  })
}

#' Run LFMM-based GEA
#'
#' @param snps SNP matrix
#' @param env Environmental variable
#' @param params Simulation parameters
#' @return Numeric vector of p-values
#' @keywords internal
run_lfmm <- function(snps, env, params) {

  if (!requireNamespace("LEA", quietly = TRUE)) {
    warning("Package 'LEA' not available, using correlation-based test instead")
    return(run_simple_correlation(snps, env, NULL))
  }

  # Determine number of latent factors if not specified
  if (is.null(params$n_latent_factors)) {
    # Simple heuristic based on population structure
    if (!is.null(params$pop_structure) && params$pop_structure != "panmictic") {
      K <- params$n_populations
    } else {
      K <- 1
    }
  } else {
    K <- params$n_latent_factors
  }

  # This is a simplified placeholder - actual LEA integration would be more complex
  # For now, fall back to correlation
  return(run_simple_correlation(snps, env, NULL))
}

#' Run simple correlation-based test
#'
#' @param snps SNP matrix
#' @param env Environmental variable
#' @param response Response variable (not used, for compatibility)
#' @return Numeric vector of p-values
#' @keywords internal
run_simple_correlation <- function(snps, env, response = NULL) {

  pvals <- apply(snps, 2, function(snp) {
    # Remove missing data
    valid <- !is.na(snp) & !is.na(env)
    if (sum(valid) < 10) return(1)

    # Correlation test
    tryCatch({
      test <- cor.test(snp[valid], env[valid], method = "spearman")
      test$p.value
    }, error = function(e) {
      1
    })
  })

  return(pvals)
}
