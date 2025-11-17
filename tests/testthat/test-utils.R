test_that("validate_simulation_parameters catches errors", {
  expect_error(validate_simulation_parameters(
    n_individuals = 10,
    n_loci = 100,
    n_causal = 150  # More causal than total loci
  ))
  
  expect_error(validate_simulation_parameters(
    n_individuals = 10,
    n_loci = 100,
    n_causal = -5  # Negative causal loci
  ))
})

test_that("calculate_required_n works", {
  # Create mock power results
  power_df <- data.frame(
    n_individuals = c(50, 100, 150),
    power = c(0.5, 0.75, 0.85),
    power_se = c(0.05, 0.04, 0.03)
  )
  class(power_df) <- c("gea_power_analysis", "data.frame")
  
  required_n <- calculate_required_n(power_df, target_power = 0.8)
  expect_equal(required_n, 150)
})