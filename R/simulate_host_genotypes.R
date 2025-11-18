#' Simulate host genotypes from a host model
#'
#' Given a \code{host_model} object (created with
#' \code{define_host_model()}), this function simulates host SNP
#' genotypes for a specified number of individuals and loci.
#'
#' For \code{type = "clone"}, all individuals are treated as having
#' effectively the same genotype (no true host variance). A trivial
#' genotype matrix is returned for compatibility.
#'
#' For \code{type = "natural"}, the \pkg{coala} package is used to
#' simulate SNPs under a coalescent model. At the moment, only a
#' panmictic population is implemented; more complex structures can
#' be added later.
#'
#' @param host_model Object of class \code{"host_model"} created by
#'   \code{define_host_model()}.
#' @param n_individuals Integer. Number of host individuals (diploid)
#'   to simulate.
#' @param n_loci Integer. Number of independent SNP loci to simulate.
#' @param maf_threshold Optional numeric between 0 and 0.5. If provided,
#'   SNPs with minor allele frequency below this threshold are
#'   filtered out.
#' @param seed Optional integer. If provided, used to set the random
#'   seed for reproducibility.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{genotypes}: numeric matrix of dimension
#'           \code{n_individuals x n_loci_filtered} with entries 0, 1, 2
#'           corresponding to minor-allele counts per individual.
#'     \item \code{host_model}: the input \code{host_model} (attached
#'           for convenience).
#'   }
#'
#' @details
#' For natural populations, this function currently implements a
#' simple panmictic coalescent model using \pkg{coala}. This is
#' intentionally conservative: instead of guessing complex migration
#' matrices, we start with a clean panmictic baseline that can be
#' refined as the project evolves.
#'
#' Technical note: the coalescent model is specified with
#' \code{ploidy = 2} so that the resulting segregating sites can
#' be interpreted as diploid genotypes (0/1/2). This assumes the
#' usual biallelic SNP setup.
#'
#' @examples
#' \dontrun{
#' library(coala)
#'
#' # Simple panmictic host
#' host_nat <- define_host_model(
#'   type              = "natural",
#'   pop_structure     = "panmictic",
#'   Ne                = 1e5,
#'   mutation_rate     = 1,
#'   recombination_rate = 0
#' )
#'
#' host_geno <- simulate_host_genotypes(
#'   host_model    = host_nat,
#'   n_individuals = 50,
#'   n_loci        = 500,
#'   maf_threshold = 0.05,
#'   seed          = 123
#' )
#'
#' dim(host_geno$genotypes)  # 50 x (<= 500) after MAF filtering
#' }
#'
#' @export
simulate_host_genotypes <- function(
  host_model,
  n_individuals,
  n_loci,
  maf_threshold = NULL,
  seed = NULL
) {
  if (!inherits(host_model, "host_model")) {
    stop("'host_model' must be an object created by define_host_model().")
  }

  if (!is.numeric(n_individuals) || length(n_individuals) != 1L ||
      n_individuals < 1) {
    stop("'n_individuals' must be a positive integer.")
  }
  n_individuals <- as.integer(n_individuals)

  if (!is.numeric(n_loci) || length(n_loci) != 1L || n_loci < 1) {
    stop("'n_loci' must be a positive integer.")
  }
  n_loci <- as.integer(n_loci)

  if (!is.null(maf_threshold)) {
    if (!is.numeric(maf_threshold) || length(maf_threshold) != 1L ||
        maf_threshold < 0 || maf_threshold > 0.5) {
      stop("'maf_threshold' must be in [0, 0.5] if provided.")
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ------------------------------------------------------------
  # Clone system: trivial genotype matrix (no host variance)
  # ------------------------------------------------------------
  if (host_model$type == "clone") {
    geno <- matrix(0L, nrow = n_individuals, ncol = n_loci)
    colnames(geno) <- paste0("locus_", seq_len(n_loci))
    rownames(geno) <- paste0("ind_", seq_len(n_individuals))

    return(list(
      genotypes  = geno,
      host_model = host_model
    ))
  }

  # ------------------------------------------------------------
  # Natural host: panmictic coalescent model via coala
  # ------------------------------------------------------------
  if (!requireNamespace("coala", quietly = TRUE)) {
    stop("Package 'coala' is required to simulate natural host genotypes.")
  }

  if (!identical(host_model$pop_structure, "panmictic")) {
    stop(
      "Currently, simulate_host_genotypes() only implements ",
      "pop_structure = 'panmictic' for natural hosts. ",
      "Support for structured populations can be added later."
    )
  }

  # Coalescent model: diploid individuals, n_loci independent loci
  model <- coala::coal_model(
    sample_size = n_individuals,
    loci_number = n_loci,
    ploidy      = 2
  ) +
    coala::feat_mutation(rate = host_model$mutation_rate) +
    coala::feat_recombination(rate = host_model$recombination_rate) +
    coala::sumstat_seg_sites()

  sim <- coala::simulate(model)
  seg_sites <- sim$seg_sites

  # seg_sites is a list of "seg_sites" objects, one per locus
  geno_list <- lapply(seg_sites, function(ss) {
    if (length(ss) == 0L) {
      return(NULL)
    }
    # as.matrix() should give a matrix: individuals x SNPs (0/1/2)
    mat <- as.matrix(ss)
    # Ensure dimensions: rows = individuals, columns = SNPs
    if (nrow(mat) != n_individuals) {
      stop("Unexpected number of rows in coala output for seg_sites.")
    }
    storage.mode(mat) <- "integer"
    mat
  })

  # Drop loci with no SNPs
  geno_list <- Filter(Negate(is.null), geno_list)

  if (length(geno_list) == 0L) {
    warning("No segregating sites simulated; returning an all-zero matrix.")
    geno <- matrix(0L, nrow = n_individuals, ncol = n_loci)
  } else {
    # Bind all SNPs across loci into one genotype matrix
    geno <- do.call(cbind, geno_list)
  }

  # Apply MAF filter if requested
  if (!is.null(maf_threshold) && ncol(geno) > 0) {
    freqs <- colMeans(geno) / 2  # allele frequency assuming 0/1/2
    maf   <- pmin(freqs, 1 - freqs)
    keep  <- maf >= maf_threshold
    if (any(!keep)) {
      geno <- geno[, keep, drop = FALSE]
    }
  }

  # Add dimnames
  if (ncol(geno) > 0) {
    colnames(geno) <- paste0("locus_", seq_len(ncol(geno)))
  }
  rownames(geno) <- paste0("ind_", seq_len(n_individuals))

  return(list(
    genotypes  = geno,
    host_model = host_model
  ))
}
