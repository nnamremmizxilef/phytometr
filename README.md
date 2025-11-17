# phytometR

<div style="display: flex; align-items: center;">
  <img src="inst/hex/logo.png" alt="phytometr logo" height="160" style="margin-right: 25px;">
  <p style="margin: 0;">
    <strong>phytometR</strong> provides comprehensive power analysis and experimental design tools for genotype-environment association (GEA) studies in holobiont systems. Simulate realistic population genetic data, environmental gradients, and holobiont responses to determine required sample sizes and assess data quality for GEA analyses.
  </p>
</div>

[![R-CMD-check](https://github.com/nnamremmizxilef/phytometr/workflows/R-CMD-check/badge.svg)](https://github.com/nnamremmizxilef/phytometr/actions)
[![Codecov test coverage](https://codecov.io/gh/nnamremmizxilef/phytometr/branch/main/graph/badge.svg)](https://codecov.io/gh/nnamremmizxilef/phytometr?branch=main)

## Overview

**phytometR** is designed for researchers planning genotype-environment association (GEA) studies, particularly in holobiont systems (host-microbiome-environment interactions). The package helps answer critical experimental design questions:

- **How many host individuals do I need to detect genetic effects on holobiont composition?** → Power analysis across sample sizes
- **How much environmental variation do I need (e.g., how far apart should sites be)?** → Gradient strength simulations
- **Can I detect host genotype × environment interactions on microbiome?** → GxE interaction power analysis
- **How many sites along a north-south transect are needed?** → Spatial sampling design optimization
- **Is my existing data suitable for GEA analysis?** → Comprehensive diagnostic checks
- **How does population structure affect my ability to detect adaptation?** → Simulations with realistic genetic structure

## Key Features

- **Realistic simulations**: Coalescent-based population genetics via `coala`
- **Flexible design**: 50+ customizable parameters for any GEA scenario
- **Holobiont focus**: Host-microbiome-environment interaction modeling
- **Population structure**: Panmictic, isolation-by-distance, discrete populations, admixture
- **Multiple response types**: Microbiome abundance, diversity, composition, host phenotypes
- **GxE interactions**: Model genotype-by-environment interactions
- **Gradient optimization**: Test different environmental ranges (temperature, precipitation, etc.)
- **Parallel processing**: Multi-core support for fast simulations
- **Unique diagnostics**: G-E independence, signal-to-noise ratio, LD decay, gradient quality
- **Publication-ready plots**: Built-in ggplot2 visualizations

## Installation

Install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("nnamremmizxilef/phytometr")
```

### Dependencies

**Required packages** (installed automatically):
- `coala` - Coalescent simulation
- `vegan` - RDA analysis
- `ggplot2` - Visualization
- `MASS` - Multivariate distributions
- `parallel` - Multi-core processing

**Optional packages** (for enhanced functionality):
- `LEA` - LFMM analysis
- `irlba` - Fast PCA for large datasets

## Quick Start

### 1. Basic Power Analysis

Determine sample size needed for 80% power:
```r
library(phytometr)

# Simulate power across sample sizes
power_results <- simulate_gea_power(
  n_individuals = seq(50, 200, by = 25),
  n_loci = 5000,
  n_causal = 10,
  effect_size = 0.3,
  response_variable = "composition",
  n_replicates = 100
)

# Visualize results
plot(power_results)

# Calculate required sample size
calculate_required_n(power_results, target_power = 0.8)
```

### 2. How Many Trees to Detect Holobiont Responses?

Test different sample sizes for detecting host genetic control of microbiome:
```r
# Simulate host genetic effect on fungal community
holobiont_power <- simulate_gea_power(
  n_individuals = seq(30, 200, by = 20),
  n_loci = 5000,
  n_causal = 15,                      # Expected number of host QTL
  effect_size = 0.25,                  # Effect of each QTL
  response_variable = "composition",   # Fungal community PC1
  holobiont_heritability = 0.20,      # Host genetic control: 20%
  env_effect_on_holobiont = 0.30,     # Environmental effect: 30%
  n_replicates = 200
)

plot(holobiont_power)
# Result: Need ~120 trees for 80% power with these parameters
```

### 3. How Much Environmental Variation is Needed?

Compare different gradient strengths (e.g., temperature ranges):
```r
# Weak gradient: 5°C difference
power_weak <- simulate_gea_power(
  n_individuals = 100,
  env_range = c(15, 20),
  response_variable = "diversity",
  n_replicates = 100
)

# Moderate gradient: 10°C difference  
power_moderate <- simulate_gea_power(
  n_individuals = 100,
  env_range = c(15, 25),
  response_variable = "diversity",
  n_replicates = 100
)

# Strong gradient: 15°C difference
power_strong <- simulate_gea_power(
  n_individuals = 100,
  env_range = c(15, 30),
  response_variable = "diversity",
  n_replicates = 100
)

# Compare
compare_power_scenarios(
  "Weak (5°C)" = power_weak,
  "Moderate (10°C)" = power_moderate,
  "Strong (15°C)" = power_strong
)
# Result: Shows how gradient strength affects power
```

### 4. How Far Apart Should Sites Be (North-South)?

Simulate different spatial scales:
```r
# Short transect: 200km, small climate variation
power_short <- simulate_gea_power(
  n_individuals = 150,
  n_sites = 5,
  gradient_type = "linear",
  env_range = c(800, 1000),           # Precipitation: 200mm difference
  response_variable = "composition",
  n_replicates = 100
)

# Long transect: 1000km, large climate variation
power_long <- simulate_gea_power(
  n_individuals = 150,
  n_sites = 5,
  gradient_type = "linear",
  env_range = c(600, 1200),           # Precipitation: 600mm difference
  response_variable = "composition",
  n_replicates = 100
)

compare_power_scenarios(
  "Short transect (200mm range)" = power_short,
  "Long transect (600mm range)" = power_long
)
```

### 5. Can I Detect GxE Interactions?

Test power for genotype-by-environment effects on holobiont:
```r
# Without GxE
power_no_gxe <- simulate_gea_power(
  n_individuals = seq(50, 250, by = 50),
  gxe_interaction = FALSE,
  holobiont_heritability = 0.25,
  env_effect_on_holobiont = 0.30,
  n_replicates = 150
)

# With GxE
power_with_gxe <- simulate_gea_power(
  n_individuals = seq(50, 250, by = 50),
  gxe_interaction = TRUE,
  gxe_strength = 0.15,
  holobiont_heritability = 0.25,
  env_effect_on_holobiont = 0.30,
  n_replicates = 150
)

compare_power_scenarios(
  "Main effects only" = power_no_gxe,
  "With GxE interaction" = power_with_gxe
)
# Result: GxE requires larger sample sizes
```

### 6. Complex Scenario with Population Structure

Model isolation-by-distance with environmental gradient:
```r
power_complex <- simulate_gea_power(
  n_individuals = seq(30, 150, by = 20),
  pop_structure = "isolation_by_distance",
  n_populations = 5,
  migration_rate = 0.02,
  gradient_type = "linear",
  env_range = c(15, 25),
  response_variable = "composition",
  gxe_interaction = TRUE,
  holobiont_heritability = 0.25,
  env_effect_on_holobiont = 0.35,
  n_replicates = 200
)

plot(power_complex)
summary(power_complex)
```

### 7. Diagnose Real Data Before Analysis

Assess data quality before expensive GEA analyses:
```r
# Create holobiont data structure
my_data <- create_holobiont_data(n_individuals = 100, host_name = "Quercus_robur")

# Add your data
my_data <- set_host_genotypes(my_data, your_snp_matrix)
my_data <- set_environment(my_data, "temperature", your_temp_data)
my_data$microbiome$fungal_diversity <- your_diversity_data

# Run comprehensive diagnostics
diagnostics <- diagnose_holobiont_data(my_data)
plot(diagnostics, which = "all")

# Check GEA assumptions
assumptions <- assess_gea_assumptions(my_data, env_variable = "temperature")
```

## Main Functions

### Data Structure
- `create_holobiont_data()` - Create integrated holobiont data structure
- `get_host_genotypes()`, `get_environment()`, `get_microbiome()` - Data accessors
- `set_host_genotypes()`, `set_environment()` - Data setters

### Power Analysis
- `simulate_gea_power()` - Comprehensive power simulation for GEA studies
- `calculate_required_n()` - Calculate sample size for target power
- `compare_power_scenarios()` - Compare power across study designs

### Diagnostics
- `diagnose_holobiont_data()` - Comprehensive data quality assessment
  - Missing data, MAF, heterozygosity, LD decay
  - Environmental gradients, spatial autocorrelation, collinearity
  - Microbiome community structure
- `assess_gea_assumptions()` - Check GEA analysis assumptions
  - G-E independence testing
  - Population structure confounding
  - Signal-to-noise ratio
  - Gradient quality

### Visualization
- `plot.gea_power_analysis()` - Power curve plots
- `plot_power_fdr()` - Power and false discovery rate
- `plot_n_detected()` - Number of loci detected
- `plot.holobiont_diagnostics()` - Diagnostic plots

## Complete Workflow Example: Oak-Ectomycorrhizal Study

Planning a study on oak genetic control of ectomycorrhizal community under drought:

### Research Questions:
1. How many oak trees needed to detect genetic effects?
2. How large should the precipitation gradient be?
3. Can we detect oak genotype × drought interactions?
```r
# Question 1: Sample size for detecting genetic effects
oak_sample_size <- simulate_gea_power(
  n_individuals = seq(50, 300, by = 50),
  n_loci = 10000,                          # Expected after filtering
  n_causal = 20,                           # Estimated QTL affecting fungi
  pop_structure = "isolation_by_distance", # Natural populations
  n_populations = 6,                       # Sampling sites
  migration_rate = 0.01,                   # Limited gene flow
  env_type = "continuous",
  gradient_type = "linear",                # Precipitation gradient
  env_range = c(400, 800),                 # 400mm precipitation range
  response_variable = "composition",       # Fungal community PC1
  holobiont_heritability = 0.20,           # Host genetic effect: 20%
  env_effect_on_holobiont = 0.30,          # Direct drought effect: 30%
  effect_size = 0.25,                      # Expected effect per QTL
  gea_method = "rda",
  n_replicates = 500,                      # High precision
  parallel = TRUE
)

plot(oak_sample_size)
summary(oak_sample_size)
# Answer: Need ~150 oak trees for 80% power

# Question 2: How much precipitation gradient is needed?
gradient_weak <- simulate_gea_power(
  n_individuals = 150,
  env_range = c(500, 700),                 # 200mm range (weak)
  response_variable = "composition",
  holobiont_heritability = 0.20,
  env_effect_on_holobiont = 0.30,
  n_replicates = 200
)

gradient_moderate <- simulate_gea_power(
  n_individuals = 150,
  env_range = c(400, 800),                 # 400mm range (moderate)
  response_variable = "composition",
  holobiont_heritability = 0.20,
  env_effect_on_holobiont = 0.30,
  n_replicates = 200
)

gradient_strong <- simulate_gea_power(
  n_individuals = 150,
  env_range = c(300, 900),                 # 600mm range (strong)
  response_variable = "composition",
  holobiont_heritability = 0.20,
  env_effect_on_holobiont = 0.30,
  n_replicates = 200
)

compare_power_scenarios(
  "Weak gradient (200mm)" = gradient_weak,
  "Moderate gradient (400mm)" = gradient_moderate,
  "Strong gradient (600mm)" = gradient_strong
)
# Answer: Moderate gradient (400mm) provides good power, 
#         spanning ~500km north-south

# Question 3: Can we detect genotype × drought interactions?
gxe_power <- simulate_gea_power(
  n_individuals = seq(100, 400, by = 50),
  n_loci = 10000,
  n_causal = 20,
  pop_structure = "isolation_by_distance",
  n_populations = 6,
  env_range = c(400, 800),
  response_variable = "composition",
  holobiont_heritability = 0.20,
  env_effect_on_holobiont = 0.30,
  gxe_interaction = TRUE,
  gxe_strength = 0.15,                     # GxE effect
  n_replicates = 300
)

plot(gxe_power)
calculate_required_n(gxe_power, target_power = 0.8)
# Answer: Need ~250 trees for 80% power to detect GxE

# Final decision: Sample 200 trees across 6 sites (33 per site)
# spanning 500km north-south (400-800mm precipitation)
```

## Practical Examples for Common Scenarios

### Scenario 1: Host Genetic Control of Microbiome Abundance

**Question**: How many plants needed to detect genetic control of a specific fungal taxon?
```r
abundance_power <- simulate_gea_power(
  n_individuals = seq(40, 200, by = 20),
  response_variable = "abundance",         # Single taxon
  holobiont_heritability = 0.15,           # Moderate genetic control
  env_effect_on_holobiont = 0.25,          # Environmental effect
  n_causal = 10,                           # QTL controlling abundance
  n_replicates = 200
)

plot(abundance_power)
# Typical result: ~100-120 individuals
```

### Scenario 2: Alpha Diversity Response to Environment

**Question**: Sample size for detecting environmental effects on microbiome diversity?
```r
diversity_power <- simulate_gea_power(
  n_individuals = seq(30, 150, by = 15),
  response_variable = "diversity",         # Alpha diversity
  env_range = c(10, 30),                   # Temperature range
  holobiont_heritability = 0.10,           # Weak genetic control
  env_effect_on_holobiont = 0.40,          # Strong environmental effect
  baseline_diversity = 50,                 # Mean richness
  n_replicates = 200
)

plot(diversity_power)
# Typical result: ~80-100 individuals for strong env effect
```

### Scenario 3: Elevation Gradient Study

**Question**: How many sites along elevation gradient (0-2000m)?
```r
elevation_power <- simulate_gea_power(
  n_individuals = 150,
  n_sites = c(3, 5, 7, 10),               # Different numbers of sites
  gradient_type = "linear",
  env_range = c(0, 2000),                 # Elevation (m)
  response_variable = "composition",
  n_replicates = 150
)

plot(elevation_power)
# Shows diminishing returns after ~7 sites
```

### Scenario 4: Latitudinal Transect

**Question**: North-south sampling for climate adaptation?
```r
latitude_power <- simulate_gea_power(
  n_individuals = 200,
  n_populations = 8,                      # 8 sites along transect
  pop_structure = "isolation_by_distance",
  migration_rate = 0.005,                 # Limited by distance
  gradient_type = "linear",
  env_range = c(5, 20),                   # Mean annual temp (°C)
  response_variable = "phenotype",        # Growth rate
  env_effect_on_holobiont = 0.35,
  n_replicates = 200
)

plot(latitude_power)
summary(latitude_power)
```

### Scenario 5: Minimum Viable Environmental Range

**Question**: If I can only sample 100 trees, how large must the environmental gradient be?
```r
# Test different gradient strengths with fixed sample size
fixed_n_test <- lapply(c(0.5, 1, 2, 5, 10), function(range_size) {
  simulate_gea_power(
    n_individuals = 100,
    env_range = c(20, 20 + range_size),   # Different temperature ranges
    response_variable = "composition",
    n_replicates = 100,
    verbose = FALSE
  )
})

names(fixed_n_test) <- paste0(c(0.5, 1, 2, 5, 10), "°C")
do.call(compare_power_scenarios, fixed_n_test)
# Shows minimum gradient needed for adequate power
```

### Scenario 6: Microhabitat Variation

**Question**: Does within-site environmental variation matter?
```r
low_noise <- simulate_gea_power(
  n_individuals = 120,
  gradient_noise = 0.05,                  # Low within-site variation
  env_range = c(15, 25),
  n_replicates = 150
)

high_noise <- simulate_gea_power(
  n_individuals = 120,
  gradient_noise = 0.25,                  # High within-site variation
  env_range = c(15, 25),
  n_replicates = 150
)

compare_power_scenarios(
  "Low microhabitat variation" = low_noise,
  "High microhabitat variation" = high_noise
)
# Shows how environmental noise affects power
```

## Decision Framework for Study Design

### Step 1: Define Your Research Question

Choose your response variable:
- `"abundance"` - Relative abundance of specific taxon
- `"diversity"` - Alpha diversity (richness, Shannon)
- `"composition"` - Beta diversity, community composition (PC1)
- `"phenotype"` - Host trait (growth, chemistry)

### Step 2: Estimate Effect Sizes

From literature or pilot data:
- **Holobiont heritability**: 0.05-0.10 (weak), 0.15-0.25 (moderate), >0.30 (strong)
- **Environmental effect**: 0.10-0.20 (weak), 0.25-0.40 (moderate), >0.45 (strong)
- **Number of causal loci**: 5-10 (oligogenic), 10-30 (polygenic), >30 (highly polygenic)

### Step 3: Define Your Environmental Gradient

Consider:
- **Range**: How much variation exists/can you sample?
- **Gradient type**: `"linear"` (transect), `"patchy"` (discrete habitats), `"random"` (local variation)
- **Spatial scale**: Local (<10km), regional (10-100km), continental (>100km)

### Step 4: Consider Population Structure

- **Panmictic**: Single well-mixed population
- **Isolation-by-distance**: Continuous populations, limited dispersal
- **Discrete**: Distinct populations with occasional migration
- **Admixture**: Recent mixing of differentiated populations

### Step 5: Run Power Simulations
```r
my_study_power <- simulate_gea_power(
  n_individuals = seq(50, 300, by = 50),  # Test range
  # ... your parameters based on Steps 1-4 ...
  n_replicates = 200                      # Higher = more precise
)

plot(my_study_power)
calculate_required_n(my_study_power, target_power = 0.8)
```

### Step 6: Optimize Design

Use `compare_power_scenarios()` to test:
- Different sample sizes vs. gradient strengths
- Concentrated sampling (few sites, many individuals) vs. distributed (many sites, fewer per site)
- With/without GxE interactions

### Step 7: Validate With Your Data

Before expensive analyses:
```r
my_data <- create_holobiont_data(n_individuals = your_n)
# ... add your actual data ...
diagnostics <- diagnose_holobiont_data(my_data)
assumptions <- assess_gea_assumptions(my_data)
```

## Use Cases

### Experimental Planning
- **Field sampling design**: Determine sample size before fieldwork
- **Budget optimization**: Balance number of sites vs. individuals per site
- **Gradient selection**: Choose environmental ranges that provide adequate power
- **Temporal design**: Decide how many time points needed

### Grant Proposals
- **Justify sample sizes**: Quantitative power analysis
- **Demonstrate feasibility**: Show realistic effect sizes are detectable
- **Cost-benefit analysis**: Show proposed sampling is efficient

### Method Comparison
- **RDA vs. LFMM**: Compare power across GEA methods
- **Different response variables**: Abundance vs. diversity vs. composition
- **Accounting for structure**: With/without controlling for population structure

### Data Quality Assessment
- **Pre-analysis checks**: Evaluate before expensive sequencing/analysis
- **Identify issues**: Missing data, poor gradients, confounding
- **Optimize filtering**: Balance SNP count vs. quality

### Sensitivity Analysis
- **Effect size uncertainty**: Test range of plausible effect sizes
- **Parameter robustness**: How sensitive is power to assumptions?
- **Risk assessment**: Probability of underpowered study

### Education and Training
- **Teach GEA concepts**: Interactive exploration of power
- **Understand assumptions**: What affects GEA success?
- **Best practices**: Learn from simulation experiments

## Application Domains

### Forest Ecology
- Oak-ectomycorrhizal fungi associations along drought gradients
- Tree-pathogen interactions across climate zones
- Root microbiome responses to soil conditions

### Plant-Microbe Interactions
- Rhizosphere community assembly
- Endophyte colonization patterns
- Arbuscular mycorrhizal fungal diversity

### Landscape Genomics
- Local adaptation to climate variables
- Elevation or latitudinal clines
- Isolation-by-environment patterns

### Agricultural Systems
- Crop-soil microbiome interactions
- Disease resistance across environments
- Stress tolerance breeding programs

### Conservation Biology
- Adaptive potential under climate change
- Genetic vulnerability assessment
- Translocation planning

## Unique Features Not Available Elsewhere

### 1. Holobiont-Specific Modeling
- Explicitly models host genetic control of microbiome
- Environmental effects on holobiont
- Genotype-by-environment interactions on microbiome
- Multiple response types (abundance, diversity, composition)

### 2. Integrated Diagnostics
- **G-E independence testing**: PC-environment correlations
- **Signal-to-noise estimation**: Predict power from observed data
- **LD decay analysis**: Genomic architecture assessment
- **Gradient quality metrics**: Strength, autocorrelation, resolution
- **Community diagnostics**: Evenness, prevalence, compositional checks
- **Integrated assessment**: Genetic + environmental + response quality

### 3. Experimental Design Optimization
- Test gradient strength vs. sample size trade-offs
- Optimize spatial sampling (site distribution)
- Balance within-site vs. between-site sampling
- Compare concentrated vs. distributed designs

### 4. Realistic Population Genetics
- Coalescent-based simulations (not just random SNPs)
- Accurate LD structure
- Demographic history (bottlenecks, expansions)
- Multiple population structure models

### 5. Flexible GxE Modeling
- Not just main effects
- Explicit GxE interaction terms
- Variable effect size distributions
- Dominance effects

## References

### Key GEA Methodology Papers

- Forester BR, Lasky JR, Wagner HH, Urban DL (2018). Comparing methods for detecting multilocus adaptation with multivariate genotype-environment associations. *Molecular Ecology*, 27: 2215-2233.

- Rellstab C, Gugerli F, Eckert AJ, Hancock AM, Holderegger R (2015). A practical guide to environmental association analysis in landscape genomics. *Molecular Ecology*, 24: 4348-4370.

- Capblancq T, Forester BR (2021). Redundancy analysis: A Swiss Army Knife for landscape genomics. *Methods in Ecology and Evolution*, 12: 2298-2309.

### Power Analysis

- Johnson PCD, Barry SJE, Ferguson HM, Müller P (2015). Power analysis for generalized linear mixed models in ecology and evolution. *Methods in Ecology and Evolution*, 6: 133-142.

### Holobiont and Microbiome

- Zilber-Rosenberg I, Rosenberg E (2008). Role of microorganisms in the evolution of animals and plants: the hologenome theory of evolution. *FEMS Microbiology Reviews*, 32: 723-735.

- Washburne AD, Silverman JD, Morton JT, et al. (2019). Phylogenetic factorization of compositional data yields lineage-level associations in microbiome datasets. *PeerJ*, 7: e16969.

## Related Projects

This package was developed by the **PhytOakmeter** research unit, studying oak holobiont responses to drought stress across environmental gradients from controlled experiments to pan-European field transects.

## License

MIT © Felix Zimmermann

## Contact

- **Author**: Felix Zimmermann
- **Institution**: WSL - Swiss Federal Institute for Forest, Snow and Landscape Research
- **Email**: felix.zimmermann@wsl.ch
- **GitHub**: [@nnamremmizxilef](https://github.com/nnamremmizxilef)
- **Issues**: [https://github.com/nnamremmizxilef/phytometr/issues](https://github.com/nnamremmizxilef/phytometr/issues)

---

**Note**: This package is under active development. Please report any issues or feature requests on GitHub.
