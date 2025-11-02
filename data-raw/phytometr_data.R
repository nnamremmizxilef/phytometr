# Create example datasets for phytometr package
set.seed(123)

# 1. Environmental data
envirodata <- data.frame(
  site_id = rep(paste0("Site_", 1:10), each = 12),
  month = rep(1:12, times = 10),
  year = 2024,
  temp_mean = rnorm(120, 15, 8) + rep(sin(seq(0, 2*pi, length.out = 12)) * 10, 10),
  temp_max = rnorm(120, 22, 10) + rep(sin(seq(0, 2*pi, length.out = 12)) * 12, 10),
  temp_min = rnorm(120, 8, 6) + rep(sin(seq(0, 2*pi, length.out = 12)) * 8, 10),
  precipitation = abs(rnorm(120, 60, 30)),
  soil_moisture = runif(120, 10, 40),
  vpd = abs(rnorm(120, 1.2, 0.6)),
  par = abs(rnorm(120, 400, 150)),
  drought_index = runif(120, 0, 1),
  herbivory_aboveground = rpois(120, 3),
  herbivory_belowground = rpois(120, 2)
)

# 2. Phenotypic data
phenotydata <- data.frame(
  tree_id = paste0("Tree_", 1:50),
  site_id = rep(paste0("Site_", 1:10), each = 5),
  clone = "DF159",
  height_cm = rnorm(50, 180, 30),
  dbh_mm = rnorm(50, 45, 10),
  leaf_area_cm2 = rnorm(50, 120, 25),
  root_biomass_g = abs(rnorm(50, 85, 20)),
  shoot_biomass_g = abs(rnorm(50, 150, 35)),
  stomatal_conductance = abs(rnorm(50, 0.3, 0.1)),
  photosynthesis_rate = abs(rnorm(50, 12, 3)),
  water_use_efficiency = abs(rnorm(50, 4.5, 1.2)),
  stress_tolerance_index = runif(50, 0, 1)
)

# 3. Genotypic data (simplified SNP-like data)
genotydata <- data.frame(
  tree_id = paste0("Tree_", 1:50),
  site_id = rep(paste0("Site_", 1:10), each = 5),
  snp_001 = sample(c(0, 1, 2), 50, replace = TRUE),
  snp_002 = sample(c(0, 1, 2), 50, replace = TRUE),
  snp_003 = sample(c(0, 1, 2), 50, replace = TRUE),
  snp_004 = sample(c(0, 1, 2), 50, replace = TRUE),
  snp_005 = sample(c(0, 1, 2), 50, replace = TRUE),
  epigenetic_mark_1 = rnorm(50, 0, 1),
  epigenetic_mark_2 = rnorm(50, 0, 1),
  gene_expression_drought = abs(rnorm(50, 100, 30)),
  gene_expression_herbivory = abs(rnorm(50, 80, 25)),
  microbiome_diversity = abs(rnorm(50, 3.5, 0.8))
)

# Save all datasets
usethis::use_data(envirodata, overwrite = TRUE)
usethis::use_data(phenotydata, overwrite = TRUE)
usethis::use_data(genotydata, overwrite = TRUE)
