test_that("diagnose_holobiont_data works", {
  holo <- create_holobiont_data(n_individuals = 50)
  snps <- matrix(sample(0:2, 500, replace = TRUE), nrow = 50, ncol = 10)
  holo <- set_host_genotypes(holo, snps)
  holo <- set_environment(holo, "temperature", rnorm(50, 20, 3))
  
  diag <- diagnose_holobiont_data(holo, verbose = FALSE)
  
  expect_s3_class(diag, "holobiont_diagnostics")
  expect_true(!is.null(diag$genetic))
  expect_true(!is.null(diag$environment))
})

test_that("diagnose_holobiont_data detects issues", {
  holo <- create_holobiont_data(n_individuals = 50)
  
  # High missing data
  snps <- matrix(sample(c(0:2, NA), 500, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.4)), 
                nrow = 50, ncol = 10)
  holo <- set_host_genotypes(holo, snps)
  
  diag <- diagnose_holobiont_data(holo, verbose = FALSE)
  expect_gt(diag$genetic$missing_rate, 0.2)
})