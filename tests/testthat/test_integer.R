context("Test newmodel class accessors.")
set.seed(1234)


test_that("Non-integer numbers are caught", {
  se <- SummarizedExperiment(matrix(rnorm(60, mean=5), nrow=10, ncol=6),
                             colData = data.frame(bio = gl(2, 3)))
  
  expect_error(m1 <- newWave(se), "The input matrix should contain only whole numbers.")
  expect_error(m2 <- newWave(se), "The input matrix should contain only whole numbers.")
})