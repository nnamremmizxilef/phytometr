#' Create holobiont data structure
#'
#' Creates a structured object to store host genetic data, environmental variables,
#' microbiome composition, and associated metadata for holobiont studies.
#'
#' @param n_individuals Integer. Number of host individuals to simulate
#' @param host_name Character. Scientific name of host species (default: "Quercus_robur")
#'
#' @return An object of class 'holobiont_data' containing:
#' \itemize{
#'   \item host_id: Character vector of unique host identifiers
#'   \item host_name: Scientific name of host
#'   \item host_genotypes: SNP matrix (individuals x loci)
#'   \item causal_loci: Integer vector of causal locus positions
#'   \item environment: List of environmental variables
#'   \item microbiome: List of microbiome data
#'   \item herbivores: List of herbivore data
#'   \item phenotypes: List of phenotypic measurements
#' }
#'
#' @examples
#' # Create data structure for 100 oak trees
#' oak_data <- create_holobiont_data(n_individuals = 100, host_name = "Quercus_robur")
#' print(oak_data)
#'
#' @export
create_holobiont_data <- function(n_individuals, host_name = "Quercus_robur") {
  
  if (!is.numeric(n_individuals) || n_individuals < 1) {
    stop("n_individuals must be a positive integer")
  }
  
  holobiont <- list(
    host_id = paste0("host_", seq_len(n_individuals)),
    host_name = host_name,
    
    # Genetic data
    host_genotypes = NULL,
    causal_loci = NULL,
    
    # Environmental data
    environment = list(
      temperature = NULL,
      precipitation = NULL,
      soil_moisture = NULL,
      latitude = NULL,
      longitude = NULL,
      elevation = NULL,
      site_id = NULL
    ),
    
    # Holobiont components
    microbiome = list(
      fungi = NULL,
      bacteria = NULL,
      fungal_diversity = NULL,
      fungal_pcs = NULL,
      fungal_community = NULL
    ),
    
    # Biotic interactions
    herbivores = list(
      herbivore_load = NULL,
      herbivore_diversity = NULL,
      specific_taxa = NULL
    ),
    
    # Phenotypes
    phenotypes = list(
      growth_rate = NULL,
      drought_tolerance = NULL,
      leaf_chemistry = NULL,
      fitness = NULL
    ),
    
    # Metadata
    sampling_date = NULL,
    site = NULL,
    experimental_treatment = NULL,
    
    # Population structure
    population = NULL,
    spatial_coordinates = NULL
  )
  
  class(holobiont) <- "holobiont_data"
  return(holobiont)
}

#' Print method for holobiont_data
#' @param x A holobiont_data object
#' @param ... Additional arguments (ignored)
#' @export
print.holobiont_data <- function(x, ...) {
  cat("Holobiont dataset\n")
  cat("================\n")
  cat("Host species:", x$host_name, "\n")
  cat("N individuals:", length(x$host_id), "\n")
  cat("\nGenetic data:\n")
  cat("  N loci:", ifelse(is.null(x$host_genotypes), 0, ncol(x$host_genotypes)), "\n")
  cat("  N causal loci:", ifelse(is.null(x$causal_loci), 0, length(x$causal_loci)), "\n")
  cat("\nEnvironmental variables:", 
      sum(!sapply(x$environment, is.null)), "/", 
      length(x$environment), "\n")
  cat("Microbiome data:", 
      ifelse(is.null(x$microbiome$fungi), "No", "Yes"), "\n")
  cat("Phenotype data:",
      sum(!sapply(x$phenotypes, is.null)), "traits\n")
  invisible(x)
}

#' Get host genotypes from holobiont data
#' @param holobiont A holobiont_data object
#' @return SNP matrix
#' @export
get_host_genotypes <- function(holobiont) {
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  holobiont$host_genotypes
}

#' Get environment data from holobiont data
#' @param holobiont A holobiont_data object
#' @param variable Character. Specific environmental variable to retrieve (optional)
#' @return List of environmental data or specific variable vector
#' @export
get_environment <- function(holobiont, variable = NULL) {
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  if (is.null(variable)) {
    holobiont$environment
  } else {
    holobiont$environment[[variable]]
  }
}

#' Get microbiome data from holobiont data
#' @param holobiont A holobiont_data object
#' @param type Character. Type of microbiome data ("fungi", "bacteria", etc.)
#' @return Microbiome data
#' @export
get_microbiome <- function(holobiont, type = "fungi") {
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  holobiont$microbiome[[type]]
}

#' Set host genotypes in holobiont data
#' @param holobiont A holobiont_data object
#' @param snp_matrix Matrix of SNP genotypes
#' @return Updated holobiont_data object
#' @export
set_host_genotypes <- function(holobiont, snp_matrix) {
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  if (nrow(snp_matrix) != length(holobiont$host_id)) {
    stop("SNP matrix rows must match number of individuals")
  }
  holobiont$host_genotypes <- snp_matrix
  holobiont
}

#' Set environment data in holobiont data
#' @param holobiont A holobiont_data object
#' @param variable Character. Name of environmental variable
#' @param values Numeric vector of values
#' @return Updated holobiont_data object
#' @export
set_environment <- function(holobiont, variable, values) {
  if (!inherits(holobiont, "holobiont_data")) {
    stop("Input must be a holobiont_data object")
  }
  if (length(values) != length(holobiont$host_id)) {
    stop("Values length must match number of individuals")
  }
  holobiont$environment[[variable]] <- values
  holobiont
}