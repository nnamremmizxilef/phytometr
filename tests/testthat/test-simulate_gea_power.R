test_that("simulate_gea_power runs with minimal parameters", {
  skip_on_cran()  # Slow test
  
  power <- simulate_gea_power(
    n_individuals = c(20, 30),
    n_loci = 100,
    n_causal = 5,
    n_replicates = 2,
    parallel = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(power, "gea_power_analysis")
  expect_true("power" %in% names(power))
  expect_equal(nrow(power), 2)
})

test_that("simulate_gea_power validates parameters", {
  expect_error(simulate_gea_power(n_individuals = -10, parallel = FALSE))
  expect_error(simulate_gea_power(n_causal = 100, n_loci = 50, parallel = FALSE))
  expect_error(simulate_gea_power(effect_size = -0.5, parallel = FALSE))
})

test_that("simulate_gea_power handles different population structures", {
  skip_on_cran()
  
  power_panmictic <- simulate_gea_power(
    n_individuals = 30,
    n_loci = 50,
    pop_structure = "panmictic",
    n_replicates = 2,
    parallel = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(power_panmictic, "gea_power_analysis")
})