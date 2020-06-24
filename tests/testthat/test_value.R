context("Test numerical correctness of functions.")
set.seed(1234)

test_that("Estimates are reasonable when data is Poisson", {
  counts <- matrix(rpois(10000, lambda=50), nrow=100, ncol=100)
  m1 <- newFit(counts)
  expect_true(all(getPhi(m1) < 1e-4))
  
  m2 <- newFit(counts, commondispersion = FALSE)
  expect_true(all(getPhi(m2) < 1e-4))
  

  expect_true(abs(mean(getMu(m1)) - 50) < 1)
  expect_true(abs(mean(getMu(m2)) - 50) < 1)
})

test_that("Estimates are reasonable when data is Negative Binomial", {
  counts <- matrix(rnbinom(10000, mu=50, size = 10), nrow=100, ncol=100)
  
  m1 <- newFit(counts, commondispersion = TRUE)
  
  expect_true(abs(mean(getMu(m1)) - 50) < 1)
  expect_true(abs(mean(getTheta(m1)) - 10) < 1)
})