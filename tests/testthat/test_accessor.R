context("Test newModel class accessors.")
set.seed(123)

test_that("dimensions", {
    x <- newmodel()
    expect_equal(numberSamples(x), NROW(newX(x)))
    expect_equal(numberFeatures(x), NROW(newV(x)))
    expect_equal(numberFactors(x), NCOL(newW(x)))
})

test_that("getters return the promised values", {
    x <- newmodel()
    expect_equal(newX(x), x@X)
    expect_equal(newV(x), x@V)
    expect_equal(newW(x), x@W)
    
    expect_equal(newBeta(x), x@beta)
    expect_equal(newGamma(x), x@gamma)
    expect_equal(newAlpha(x), x@alpha)
    
    expect_equal(newZeta(x), x@zeta)
    expect_equal(newTheta(x), exp(x@zeta))
    expect_equal(newPhi(x), exp(-x@zeta))
})