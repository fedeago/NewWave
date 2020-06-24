context("Test newModel class accessors.")
set.seed(123)

test_that("dimensions", {
  x <- newmodel()
  expect_equal(nSamples(x), NROW(getX(x)))
  expect_equal(nFeatures(x), NROW(getV(x)))
  expect_equal(nFactors(x), NCOL(getW(x)))
})

test_that("getters return the promised values", {
  x <- newmodel()
  expect_equal(getX(x), x@X)
  expect_equal(getV(x), x@V)
  expect_equal(getW(x), x@W)
  
  expect_equal(getBeta(x), x@beta)
  expect_equal(getGamma(x), x@gamma)
  expect_equal(getAlpha(x), x@alpha)
  
  expect_equal(getZeta(x), x@zeta)
  expect_equal(getTheta(x), exp(x@zeta))
  expect_equal(getPhi(x), exp(-x@zeta))
})