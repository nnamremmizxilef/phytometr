test_that("create_holobiont_data works", {
  holo <- create_holobiont_data(n_individuals = 50)
  
  expect_s3_class(holo, "holobiont_data")
  expect_equal(length(holo$host_id), 50)
  expect_null(holo$host_genotypes)
  expect_true(is.list(holo$environment))
})

test_that("create_holobiont_data validates input", {
  expect_error(create_holobiont_data(n_individuals = -5))
  expect_error(create_holobiont_data(n_individuals = "fifty"))
})

test_that("setters work correctly", {
  holo <- create_holobiont_data(n_individuals = 10)
  snps <- matrix(sample(0:2, 100, replace = TRUE), nrow = 10, ncol = 10)
  
  holo <- set_host_genotypes(holo, snps)
  expect_equal(dim(holo$host_genotypes), c(10, 10))
  
  temp <- rnorm(10, mean = 20, sd = 3)
  holo <- set_environment(holo, "temperature", temp)
  expect_equal(length(holo$environment$temperature), 10)
})

test_that("setters validate dimensions", {
  holo <- create_holobiont_data(n_individuals = 10)
  wrong_snps <- matrix(0, nrow = 5, ncol = 10)  # Wrong number of individuals
  
  expect_error(set_host_genotypes(holo, wrong_snps))
  expect_error(set_environment(holo, "temp", 1:5))  # Wrong length
})