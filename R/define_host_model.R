#' Define the host genetic model for simulations
#'
#' This function specifies how host genotypes should be simulated:
#' either as a clonal system (no or negligible genotypic variation),
#' or as a natural population simulated under a coalescent model
#' using the \pkg{coala} package.
#'
#' It does **not** generate any genotypes yet. Instead, it stores
#' the assumptions and parameters needed to simulate host SNP data
#' later, when a specific sample size (number of individuals and
#' loci) is chosen.
#'
#' @param type Character. Either:
#'   \itemize{
#'     \item \code{"clone"}: clonal system; effectively one genotype,
#'           all individuals genetically identical for our purposes.
#'     \item \code{"natural"}: natural population, genotypes simulated
#'           with a coalescent model via \pkg{coala}.
#'   }
#' @param pop_structure Character (for \code{type = "natural"}). High-level
#'   description of population structure, currently one of:
#'   \itemize{
#'     \item \code{"panmictic"}: single well-mixed population.
#'     \item \code{"discrete"}: (reserved for future use, not yet implemented).
#'     \item \code{"isolation_by_distance"}: (reserved for future use).
#'   }
#'   At the moment, only \code{"panmictic"} is implemented for the
#'   coalescent simulation; other structures will give an error to
#'   avoid pretending it works.
#' @param Ne Numeric. Effective population size. Stored mainly for
#'   documentation; coalescent scaling is currently driven by
#'   \code{mutation_rate} and \code{recombination_rate}.
#' @param mutation_rate Numeric. Mutation rate parameter passed to
#'   \pkg{coala} (often denoted \eqn{\theta}). Adjust based on expected
#'   SNP density.
#' @param recombination_rate Numeric. Recombination rate parameter
#'   for \pkg{coala}. Set to 0 to ignore recombination.
#' @param migration_rate Numeric. Placeholder for future structured
#'   models. Currently not used for \code{"panmictic"}.
#' @param fst Optional numeric. Target Fst between populations. Stored
#'   as metadata for now; the mapping from Fst to coalescent parameters
#'   can be implemented later.
#' @param var_host Numeric in [0, 1]. Proportion of variance in the
#'   holobiont/phenotypic response that is expected to be explained
#'   by host genotype (based on literature or pilot data). This is
#'   used later when you simulate the holobiont response.
#' @param notes Optional character string for free-text comments
#'   (e.g. literature references, assumptions about Ne or Fst).
#'
#' @return A list of class \code{"host_model"} with elements:
#'   \itemize{
#'     \item \code{type}: "clone" or "natural".
#'     \item \code{pop_structure}: population structure (if natural).
#'     \item \code{Ne}, \code{mutation_rate}, \code{recombination_rate},
#'           \code{migration_rate}, \code{fst}: coalescent parameters.
#'     \item \code{var_host}: expected variance proportion explained by host.
#'     \item \code{notes}: user comments.
#'   }
#'
#' @examples
#' # 1) Clone system: one genotype, no host variance
#' host_clone <- define_host_model(
#'   type     = "clone",
#'   var_host = 0,
#'   notes    = "Clonal oak line; host variance assumed negligible."
#' )
#' host_clone
#'
#' # 2) Natural panmictic population: coalescent host
#' host_nat <- define_host_model(
#'   type              = "natural",
#'   pop_structure     = "panmictic",
#'   Ne                = 1e5,
#'   mutation_rate     = 1,
#'   recombination_rate = 0,
#'   var_host          = 0.15,
#'   notes             = "Panmictic host; h2 ~ 0.15 based on literature."
#' )
#' host_nat
#'
#' @export
define_host_model <- function(
  type              = c("clone", "natural"),
  pop_structure     = c("panmictic", "discrete", "isolation_by_distance"),
  Ne                = 1e5,
  mutation_rate     = 1,
  recombination_rate = 0,
  migration_rate    = 0.01,
  fst               = NULL,
  var_host          = 0,
  notes             = NULL
) {
  type <- match.arg(type)

  if (!is.numeric(var_host) || length(var_host) != 1L ||
      var_host < 0 || var_host > 1) {
    stop("'var_host' must be a single number between 0 and 1.")
  }

  if (type == "clone") {
    pop_structure <- NA_character_
  } else {
    pop_structure <- match.arg(pop_structure)
  }

  model <- list(
    type              = type,
    pop_structure     = pop_structure,
    Ne                = Ne,
    mutation_rate     = mutation_rate,
    recombination_rate = recombination_rate,
    migration_rate    = migration_rate,
    fst               = fst,
    var_host          = var_host,
    notes             = notes
  )

  class(model) <- "host_model"
  return(model)
}
