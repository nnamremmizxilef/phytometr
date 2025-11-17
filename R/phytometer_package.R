#' phytometr: Power Analysis for Genotype-Environment Association Studies in Holobiont Systems
#'
#' A comprehensive framework for experimental design and power analysis in 
#' genotype-environment association (GEA) studies, with specialized support
#' for holobiont systems (host-microbiome interactions). The package simulates
#' realistic population genetic data using coalescent theory, environmental
#' gradients, and holobiont responses to determine required sample sizes for
#' detecting genotype-environment associations.
#'
#' @section Main functions:
#' 
#' **Data Structure:**
#' \itemize{
#'   \item \code{\link{create_holobiont_data}} - Create integrated holobiont data structure
#'   \item \code{\link{get_host_genotypes}}, \code{\link{get_environment}}, \code{\link{get_microbiome}} - Data accessors
#'   \item \code{\link{set_host_genotypes}}, \code{\link{set_environment}} - Data setters
#' }
#' 
#' **Power Analysis:**
#' \itemize{
#'   \item \code{\link{simulate_gea_power}} - Simulate power for GEA studies
#'   \item \code{\link{calculate_required_n}} - Calculate required sample size for target power
#'   \item \code{\link{compare_power_scenarios}} - Compare power across study designs
#' }
#' 
#' **Diagnostics:**
#' \itemize{
#'   \item \code{\link{diagnose_holobiont_data}} - Comprehensive data quality assessment
#'   \item \code{\link{assess_gea_assumptions}} - Check GEA analysis assumptions
#' }
#' 
#' **Visualization:**
#' \itemize{
#'   \item \code{plot()} methods for power results and diagnostics
#'   \item \code{\link{plot_power_fdr}} - Plot power and false discovery rate
#'   \item \code{\link{plot_n_detected}} - Plot number of loci detected
#' }
#'
#' @section Workflow:
#' 
#' **1. Basic power analysis:**
#' ```
#' # Determine sample size needed for 80% power
#' power_results <- simulate_gea_power(
#'   n_individuals = seq(50, 200, by = 25),
#'   n_loci = 5000,
#'   n_causal = 10,
#'   effect_size = 0.3,
#'   response_variable = "composition"
#' )
#' plot(power_results)
#' calculate_required_n(power_results, target_power = 0.8)
#' ```
#' 
#' **2. Complex scenario with population structure:**
#' ```
#' power_complex <- simulate_gea_power(
#'   n_individuals = seq(30, 150, by = 20),
#'   pop_structure = "isolation_by_distance",
#'   n_populations = 5,
#'   gradient_type = "linear",
#'   gxe_interaction = TRUE,
#'   n_replicates = 200
#' )
#' ```
#' 
#' **3. Diagnose real data:**
#' ```
#' # Create holobiont data structure
#' my_data <- create_holobiont_data(n_individuals = 100)
#' my_data <- set_host_genotypes(my_data, snp_matrix)
#' my_data <- set_environment(my_data, "temperature", temp_values)
#' 
#' # Run diagnostics
#' diag <- diagnose_holobiont_data(my_data)
#' plot(diag)
#' 
#' # Check GEA assumptions
#' assumptions <- assess_gea_assumptions(my_data, env_variable = "temperature")
#' ```
#'
#' @section Key Features:
#' \itemize{
#'   \item **Realistic simulations**: Coalescent-based population genetics via coala
#'   \item **Flexible design**: 50+ customizable parameters for any GEA scenario
#'   \item **Holobiont focus**: Explicitly models host-microbiome-environment interactions
#'   \item **Population structure**: Supports panmictic, IBD, discrete populations, admixture
#'   \item **Multiple response types**: Abundance, diversity, composition, phenotypes
#'   \item **GxE interactions**: Model genotype-by-environment interactions
#'   \item **Parallel processing**: Multi-core support for fast simulations
#'   \item **Comprehensive diagnostics**: Unique quality checks not available elsewhere
#'   \item **Publication-ready plots**: Built-in ggplot2 visualizations
#' }
#'
#' @section Unique Diagnostic Features:
#' \itemize{
#'   \item G-E independence testing via genetic PC-environment correlations
#'   \item Signal-to-noise ratio estimation for power assessment
#'   \item Linkage disequilibrium decay rate estimation
#'   \item Environmental gradient strength and spatial autocorrelation
#'   \item Multi-variable collinearity with VIF calculation
#'   \item Holobiont community diagnostics (evenness, prevalence, composition)
#'   \item Population structure confounding detection
#'   \item Integrated assessment combining genetic, environmental, and response data
#' }
#'
#' @section Package dependencies:
#' 
#' **Required packages:**
#' \itemize{
#'   \item coala - Coalescent simulation for realistic population genetics
#'   \item vegan - RDA (Redundancy Analysis) for GEA
#'   \item ggplot2 - Publication-quality visualizations
#'   \item MASS - Multivariate normal distributions
#'   \item parallel - Multi-core processing
#' }
#' 
#' **Suggested packages (enhanced functionality):**
#' \itemize{
#'   \item LEA - LFMM (Latent Factor Mixed Models) analysis
#'   \item irlba - Fast PCA via randomized SVD for large genetic datasets
#' }
#'
#' @section Use Cases:
#' \itemize{
#'   \item **Experimental planning**: Determine required sample size before fieldwork
#'   \item **Grant proposals**: Justify sampling design with quantitative power analysis
#'   \item **Method comparison**: Compare statistical power across GEA methods (RDA, LFMM, etc.)
#'   \item **Data quality assessment**: Evaluate real datasets before expensive analyses
#'   \item **Sensitivity analysis**: Test how effect size, gradient strength affect power
#'   \item **Holobiont studies**: Optimize designs for host-microbiome association studies
#'   \item **Education**: Understand factors affecting GEA power and assumptions
#' }
#'
#' @section Application Examples:
#' 
#' **Ecological genomics:**
#' - Oak-ectomycorrhizal fungi associations along drought gradients
#' - Plant-pathogen interactions across temperature gradients
#' - Host-microbiome responses to environmental stress
#' 
#' **Landscape genomics:**
#' - Local adaptation to climate variables
#' - Elevation or latitudinal clines
#' - Isolation-by-environment patterns
#' 
#' **Microbiome studies:**
#' - Host genetic control of microbiome composition
#' - Environmental filtering of microbial communities
#' - Tripartite interactions (host × environment → microbiome)
#'
#' @section Citation:
#' To cite phytometr in publications, use:
#' 
#'   Felix Zimmermann (2025). phytometr: Power Analysis for Genotype-Environment
#'   Association Studies in Holobiont Systems. R package version 0.1.0.9000.
#'   https://github.com/nnamremmizxilef/phytometr
#'
#' @section Getting Help:
#' \itemize{
#'   \item Package documentation: \code{help(package = "phytometr")}
#'   \item Function help: \code{?simulate_gea_power}
#'   \item Vignettes: \code{browseVignettes("phytometr")}
#'   \item GitHub issues: \url{https://github.com/nnamremmizxilef/phytometr/issues}
#' }
#'
#' @docType package
#' @name phytometr-package
#' @aliases phytometr
#' 
#' @importFrom stats aov coef cor cor.test kruskal.test lm median nls nls.control
#' @importFrom stats p.adjust plogis pnorm qlogis rbinom rexp rnorm runif sd var
#' @importFrom stats as.formula
#' @importFrom parallel clusterEvalQ clusterExport detectCores makeCluster parLapply stopCluster
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom graphics plot
#' @importFrom grDevices dev.off pdf
#' @importFrom coala coal_model feat_mutation feat_recombination feat_migration feat_size_change feat_growth
#' 
#' @references
#' Forester BR, Lasky JR, Wagner HH, Urban DL (2018). 
#' Comparing methods for detecting multilocus adaptation with multivariate 
#' genotype-environment associations. Molecular Ecology, 27: 2215-2233.
#' 
#' Rellstab C, Gugerli F, Eckert AJ, Hancock AM, Holderegger R (2015).
#' A practical guide to environmental association analysis in landscape genomics.
#' Molecular Ecology, 24: 4348-4370.
#' 
#' Capblancq T, Forester BR (2021). Redundancy analysis: A Swiss Army Knife 
#' for landscape genomics. Methods in Ecology and Evolution, 12: 2298-2309.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL