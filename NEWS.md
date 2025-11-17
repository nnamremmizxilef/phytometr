# phytometr 0.1.0.9000

## Major Refactoring - GEA Power Analysis Framework

### Breaking Changes
* Removed mock data generation functions (`envirodata()`, `phenotydata()`, `genotydata()`)
* Replaced with comprehensive holobiont simulation and power analysis framework
* New focus: genotype-environment association (GEA) study design and power analysis

### New Core Data Structure
* Added `create_holobiont_data()` - Unified data structure for holobiont studies
  - Integrates host genetics, environment, microbiome, herbivores, and phenotypes
  - Single object tracks all components with consistent host IDs
  - S3 class with print, accessor, and setter methods
* Added `get_host_genotypes()`, `get_environment()`, `get_microbiome()` - Data accessors
* Added `set_host_genotypes()`, `set_environment()` - Data setters

### Simulation Framework
* Added `simulate_gea_power()` - Comprehensive power analysis for GEA studies
  - Tests multiple sample sizes in parallel
  - Highly customizable with 50+ parameters
  - Supports panmictic, IBD, discrete, and admixture population structures
  - Multiple response types: abundance, diversity, composition, phenotype
  - GxE interaction modeling
  - Returns `gea_power_analysis` class with built-in plotting

### Genetic Simulation
* Added `simulate_host_genetics()` - Population genetic simulation via coala
  - Coalescent-based neutral evolution
  - Population structure (panmictic, IBD, discrete populations)
  - Demographic history (bottlenecks, expansions)
  - Migration and gene flow
  - Recombination and mutation
  - MAF filtering and missing data simulation
* Added `convert_to_snp_matrix()` - Convert coala output to SNP matrices

### Environmental Simulation
* Added `simulate_environment_data()` - Environmental gradient generation
  - Continuous, categorical, and multivariate environments
  - Gradient types: linear, nonlinear, patchy, random
  - Spatial autocorrelation
  - Environmental noise
  - Multi-variable correlations
* Added `generate_continuous_environment()` - Continuous gradient builder
* Added `generate_multivariate_environment()` - Correlated environmental variables

### Holobiont Response Simulation
* Added `simulate_response()` - Holobiont response to genetics and environment
  - Host genetic control of microbiome (heritability)
  - Direct environmental effects
  - GxE interactions
  - Multiple response types (abundance, diversity, composition, phenotype)
* Added `generate_effect_sizes()` - Effect size distributions (equal, exponential, normal)
* Added `calculate_genetic_value()` - Breeding value calculation with dominance
* Added `simulate_community_matrix()` - Full microbiome community simulation

### GEA Analysis Methods
* Added `run_gea_analysis()` - GEA method implementations
  - RDA (Redundancy Analysis)
  - LFMM (Latent Factor Mixed Models) - framework ready
  - Correlation-based tests
  - Multiple testing correction (Bonferroni, FDR, q-value)
* Added `run_rda()` - RDA implementation using vegan
* Added `run_simple_correlation()` - Spearman correlation tests

### Power Calculation
* Added `run_single_replicate()` - Single simulation replicate
* Added `calc_power_metrics()` - Calculate power, FDR, true/false positives
* Added `aggregate_replicates()` - Aggregate across simulation replicates
* Added `format_power_results()` - Format results as `gea_power_analysis` object
* Added `calculate_required_n()` - Calculate sample size for target power

### Comprehensive Diagnostics
* Added `diagnose_holobiont_data()` - Complete data quality assessment
  - Genetic data: missing data, MAF, heterozygosity, LD decay
  - Environment: gradient strength, spatial autocorrelation, collinearity (VIF)
  - Holobiont response: abundance, diversity, composition metrics
  - Population structure: Fst calculation
  - Returns `holobiont_diagnostics` class with plotting
* Added `diagnose_genetic_data()` - Genetic data quality checks
  - Missing data rates per locus and individual
  - Allele frequency spectrum
  - Heterozygosity outlier detection
  - Linkage disequilibrium estimation
* Added `estimate_ld_decay()` - LD decay rate estimation
* Added `calculate_r2()` - Pairwise LD calculation
* Added `diagnose_environmental_data()` - Environmental data assessment
  - Gradient detection and strength
  - Spatial autocorrelation (Moran's I)
  - Multi-variable collinearity (VIF)
* Added `assess_gradient_strength()` - Monotonic trend testing
* Added `assess_env_collinearity()` - VIF calculation for environmental variables
* Added `diagnose_holobiont_response()` - Response variable diagnostics
  - Distribution statistics
  - Zero-inflation detection
  - Community matrix diagnostics (prevalence, evenness)
* Added `diagnose_community_matrix()` - Microbiome community quality
  - Taxa prevalence and abundance
  - Pielou's evenness
  - Compositional data checks
* Added `diagnose_population_structure()` - Population structure assessment
  - Population size balance
  - Between-population Fst
* Added `calculate_fst()` - Weir & Cockerham Fst estimation

### GEA-Specific Diagnostics
* Added `assess_gea_assumptions()` - Comprehensive GEA assumption testing
  - Genotype-environment independence via PCA
  - Population structure confounding detection
  - Signal-to-noise ratio estimation
  - Outlier detection (genetic and environmental)
  - Gradient quality assessment
  - Method-specific checks (RDA, LFMM)
  - Returns `gea_assumptions` class with recommendations
* Added `test_ge_independence()` - Test G-E correlation via genetic PCs
* Added `fast_pca()` - Efficient PCA for genetic data (with SVD)
* Added `assess_structure_confounding()` - Environmental variation across populations
* Added `estimate_signal_noise()` - SNR for power estimation
* Added `detect_outliers()` - Statistical outlier detection (3-SD rule)
* Added `assess_gradient_quality()` - Environmental gradient metrics
* Added `check_rda_assumptions()` - RDA-specific assumption checks
* Added `check_lfmm_assumptions()` - LFMM-specific checks with K estimation
* Added `compare_power_scenarios()` - Compare power across study designs

### Visualization
* Added `plot.gea_power_analysis()` - Power curve plotting
  - Sample size vs. power
  - Error bars for uncertainty
  - Target power line
  - Publication-ready ggplot2 graphics
* Added `plot_power_fdr()` - Dual plot of power and FDR
* Added `plot_n_detected()` - Number of loci detected vs. sample size
* Added `plot.holobiont_diagnostics()` - Diagnostic plots
  - MAF distribution
  - Heterozygosity per individual
  - Missing data patterns
  - Environmental correlations
* Added `plot_maf_distribution()` - Allele frequency spectrum
* Added `plot_heterozygosity()` - Individual heterozygosity with outliers
* Added `plot_missing_data()` - Missing data visualization
* Added `plot_environmental_diagnostics()` - Correlation heatmap with significance

### Utility Functions
* Added `validate_simulation_parameters()` - Parameter validation and warnings
* Added `calculate_data_completeness()` - Overall data completeness score
* Added `assess_sample_size()` - Context-specific sample size assessment
* Added `calculate_evenness()` - Pielou's evenness for communities
* Added `calculate_autocorrelation()` - Lag-1 autocorrelation
* Added `estimate_spatial_autocorrelation()` - Moran's I approximation
* Added `reshape_cor_matrix()` - Format correlation matrices for ggplot

### S3 Methods
* Added `print.holobiont_data()` - Print method for holobiont objects
* Added `print.holobiont_diagnostics()` - Detailed diagnostic output
* Added `print.gea_power_analysis()` - Power analysis summary
* Added `print.gea_assumptions()` - GEA assumption check results
* Added `summary.gea_power_analysis()` - Summary with sample size recommendation
* Added `plot.gea_power_analysis()` - Power visualization
* Added `plot.holobiont_diagnostics()` - Diagnostic visualization
* Added `plot.gea_assumptions()` - Assumption check plots

### Dependencies
* Added: coala (coalescent simulation)
* Added: vegan (RDA analysis)
* Added: MASS (multivariate distributions)
* Added: parallel (parallel processing)
* Suggested: LEA (LFMM analysis)
* Suggested: irlba (fast SVD for large matrices)
* Core: ggplot2, dplyr, tidyr
* Removed: corrplot, readr, writexl (from previous version)

### Key Features
* **Realistic simulations**: Coalescent-based genetics with population structure
* **Flexible design**: 50+ customizable parameters for any GEA scenario
* **Holobiont focus**: Explicitly models host-microbiome-environment interactions
* **Parallel processing**: Multi-core support for fast power calculations
* **Comprehensive diagnostics**: Unique checks not available in other packages
  - G-E independence testing
  - Signal-to-noise ratio estimation
  - Integrated data quality across all components
  - Method-specific assumption validation
* **Publication-ready**: Built-in visualizations with ggplot2

### Use Cases
* **Experimental planning**: Determine required sample size for GEA studies
* **Grant proposals**: Justify sampling design with quantitative power analysis
* **Method comparison**: Compare statistical power across GEA methods
* **Data quality**: Assess real datasets before expensive analyses
* **Education**: Understand factors affecting GEA power (effect size, structure, gradient)

### Documentation
* Complete function documentation with roxygen2
* Detailed parameter descriptions
* Extensive examples for common scenarios
* Updated README with new workflow
* Package focuses on: genotype-environment associations, holobiont biology, power analysis

### Removed Functions (from 0.0.0.9000)
* `envirodata()`, `phenotydata()`, `genotydata()` - Replaced by holobiont framework
* `showtidydata()`, `createtidydata()`, `savetidydata()` - Out of scope
* `mergeomics()`, `calculatestress()`, `compareconditions()` - Out of scope
* `aggregatetemporal()`, `detectoutliers()`, `correlatephenoenv()` - Out of scope
* `runpca()`, `summarizebygroup()` - Available in other packages

---

# phytometr 0.0.0.9000

## Initial Release (Deprecated)
* Initial development version - basic data generation and exploration
* See version 0.1.0.9000 for current functionality

## Note
Version 0.0.0.9000 contained mock data generation functions that have been completely
replaced with a comprehensive GEA power analysis framework in version 0.1.0.9000.
The package focus has shifted from general data exploration to specialized power
analysis for genotype-environment association studies in holobiont systems.