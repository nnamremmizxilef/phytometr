# phytometr 0.0.0.9000

## Initial Release

* Initial development version of phytometR package
* Package setup with MIT license and GitHub repository

## Data Generation Functions

* Added `envirodata()` - Generate customizable environmental monitoring data
  - Options for number of sites, months, temperature range, precipitation
* Added `phenotydata()` - Generate phenotypic measurements from oak saplings
  - Customizable tree numbers, sites, height ranges
* Added `genotydata()` - Generate multi-omics data
  - SNP markers, epigenetic marks, gene expression, microbiome diversity

## Data Exploration Functions

* Added `showtidydata()` - Exploratory data analysis
  - Distribution plots and correlation heatmaps
  - Pearson and Spearman correlation options
* Added `createtidydata()` - Convert data to tidy long format
  - Flexible variable selection and column naming
  - Optional timestamp and sorting
* Added `savetidydata()` - Save tidy data to file
  - Multiple formats: csv, tsv, rds, xlsx
  - Preserves identifier columns

## Analysis Functions

* Added `mergeomics()` - Merge environmental, phenotypic, and genotypic datasets
* Added `calculatestress()` - Calculate composite stress indices
  - Drought, heat, cold, and herbivory stress types
* Added `compareconditions()` - Statistical comparison across groups
  - ANOVA, Kruskal-Wallis, and t-test options
  - Automated visualization
* Added `aggregatetemporal()` - Aggregate data over time periods
  - Weekly, monthly, seasonal, yearly aggregation
* Added `detectoutliers()` - Identify outliers using multiple methods
  - IQR, z-score, and MAD methods
* Added `correlatephenoenv()` - Correlate phenotype with environment
  - Correlation matrix and heatmap visualization
* Added `runpca()` - Principal Component Analysis
  - Scree plots and biplots
  - Customizable number of components
* Added `summarizebygroup()` - Calculate grouped summary statistics
  - Multiple summary functions (mean, sd, median, min, max, n)

## Dependencies

* Core: dplyr, tidyr, ggplot2, corrplot, rlang, readr
* Suggested: writexl (for Excel export)

## Documentation

* Created package hex sticker logo
* Added comprehensive README with installation and examples
* Complete function documentation with examples
* All functions tested with example datasets