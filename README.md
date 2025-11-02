# phytometR

<div style="display: flex; align-items: center;">

  <img src="inst/hex/logo.png" alt="phytometR logo" height="160" style="margin-right: 25px;">

  <p style="margin: 0;">
    <strong>phytometR</strong> provides tools for analyzing multi-omics data integrated with environmental measurements in ecological and evolutionary studies.  
    The package facilitates investigation of organism responses to environmental stressors across spatial and temporal gradients.
  </p>

</div>





## Features

- **Data generation**: Create customizable example datasets for testing and development
- **Multi-omics integration**: Merge genomic, phenotypic, and environmental data
- **Data transformation**: Convert between wide and tidy formats with flexible options
- **Environmental analysis**: Calculate stress indices, aggregate temporal data, detect outliers
- **Statistical analysis**: Compare conditions, correlate traits with environment, run PCA
- **Visualization**: Distribution plots, correlation heatmaps, PCA biplots, comparison plots

## Installation

Install the development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("nnamremmizxilef/phytometr")
```

## Quick Start

### Generate Example Data
```r
library(phytometr)

# Create datasets with default parameters
env <- envirodata()
pheno <- phenotydata()
geno <- genotydata()

# Customize data generation
env_large <- envirodata(n_sites = 20, n_months = 24, temp_range = c(0, 30))
pheno_small <- phenotydata(n_trees = 100, height_range = c(80, 120))
geno_dense <- genotydata(n_trees = 50, n_snps = 20)
```

### Explore and Transform Data
```r
# Visualize distributions and correlations
showtidydata(env)

# Convert to tidy format
tidy_env <- createtidydata(env, keep_cols = "site_id")

# Save tidy data
savetidydata(env, "environment_tidy.csv", keep_cols = c("site_id", "month"))
```

### Merge Multi-Omics Data
```r
# Merge all three data types
merged <- mergeomics(enviro = env, pheno = pheno, geno = geno)
head(merged)
```

### Analyze Stress Responses
```r
# Calculate drought stress index
env_stress <- calculatestress(env, stress.type = "drought")

# Compare sites statistically
comparison <- compareconditions(env, 
                                group.by = "site_id",
                                variables = c("temp_mean", "precipitation"))
print(comparison$results)
print(comparison$plots$temp_mean)
```

### Correlate Phenotype with Environment
```r
# Link traits to environmental conditions
correlations <- correlatephenoenv(pheno, env)
print(correlations$correlation_matrix)
print(correlations$plot)
```

### Dimensionality Reduction
```r
# Run PCA on genotypic data
pca_results <- runpca(geno, n.components = 5, color.by = "site_id")
print(pca_results$plots$scree)
print(pca_results$plots$biplot)
```

### Temporal Aggregation
```r
# Aggregate environmental data monthly
monthly <- aggregatetemporal(env, 
                             time.col = "month",
                             variables = c("temp_mean", "precipitation"),
                             group.by = "site_id",
                             agg.fun = "mean")
```

### Summary Statistics
```r
# Calculate grouped summaries
stats <- summarizebygroup(env,
                         group.vars = "site_id",
                         summary.vars = c("temp_mean", "precipitation"),
                         funs = c("mean", "sd", "min", "max"))
```

## Main Functions

### Data Generation
- `envirodata()` - Generate environmental monitoring data
- `phenotydata()` - Generate phenotypic measurements
- `genotydata()` - Generate genotypic/omics data

### Data Transformation
- `showtidydata()` - Visualize distributions and correlations
- `createtidydata()` - Convert to tidy long format
- `savetidydata()` - Save tidy data to file

### Data Integration
- `mergeomics()` - Merge multiple data types

### Analysis
- `calculatestress()` - Calculate stress indices
- `compareconditions()` - Statistical comparisons across groups
- `aggregatetemporal()` - Temporal aggregation
- `detectoutliers()` - Outlier detection
- `correlatephenoenv()` - Phenotype-environment correlations
- `runpca()` - Principal Component Analysis
- `summarizebygroup()` - Grouped summary statistics

## Use Cases

- Analyzing acclimation and adaptation mechanisms
- Integrating multi-omics datasets with environmental data
- Investigating stress responses (drought, herbivory, temperature)
- Studying organism-microbiome interactions
- Linking field monitoring with controlled experiments
- Quality control and data exploration

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT © Felix Zimmermann

## Citation

If you use phytometR in your research, please cite:
```
Zimmermann, F. (2025). phytometR: Analysis Tools for Tree Holobiont 
Acclimation and Adaptation Research. R package version 0.0.0.9000.
https://github.com/nnamremmizxilef/phytometr
```