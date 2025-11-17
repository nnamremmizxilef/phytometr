# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  # ggplot2 aesthetics
  "n_individuals",
  "power",
  "power_se",
  "fdr",
  "fdr_se", 
  "n_detected",
  "scenario",
  "value",
  "metric",
  "Var1",
  "Var2",
  "individual",
  "heterozygosity",
  "missing_rate",
  "maf",
  "se",
  
  # Internal functions
  "run_gea_analysis"
))