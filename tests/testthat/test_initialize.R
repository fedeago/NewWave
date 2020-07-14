context("Test simulation, initialization and epsilon.")
set.seed(1234)


test_that("W can be estimated from random matrix (no signal)", {

    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(rnorm(nS*K),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(rnorm(nG), ncol = 1)
    mm = newmodel(X = X,
                    V = V,
                    W = W,
                    gamma = matrix(1, ncol = nS),
                    beta = matrix(5, ncol=nG))

    my_data <- newSim(mm)

    expect_equal(W, newW(mm))
    expect_equal(X, newX(mm))
    expect_equal(V, newV(mm))

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts),
                                rowDat = data.frame(V),
                                colDat = data.frame(X))
    sf = newFit(SE,
                V = "~V - 1",
                K = 2)

    round(cor(cbind(W, sf@W)), 2)
})

test_that("W can be estimated from two-dimensional signal (no signal)", {
    Sys.setenv("R_TESTS" = "")
    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(c(rnorm(50, mean=5, sd=.1), rnorm(50, mean=1, sd=.1),
                rnorm(10, mean=5, sd=.1), rnorm(90, mean=1, sd=.1)),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(rnorm(nG), ncol = 1)
    mm = newmodel(X = X,
                V = V,
                W = W,
                gamma = matrix(1, ncol = nS),
                beta = matrix(5, ncol=nG))

    my_data <- newSim(mm)

    expect_equal(W, newW(mm))
    expect_equal(X, newX(mm))
    expect_equal(V, newV(mm))

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts),
                            rowDat = data.frame(V),
                            colDat = data.frame(X))
    sf = newFit(SE,
                V = "~V - 1",
                K = 2)

    round(cor(cbind(W, sf@W)), 2)
})

test_that("Initialization works with large epsilon", {

    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(c(rnorm(50, mean=5, sd=.1), rnorm(50, mean=1, sd=.1),
                rnorm(10, mean=5, sd=.1), rnorm(90, mean=1, sd=.1)),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(1, ncol = 1, nrow = nG)
    mm = newmodel(X = X,
                V = V,
                W = W,
                gamma = matrix(1, ncol = nS),
                beta = matrix(5, ncol=nG))

    my_data <- newSim(mm)

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts))

    expect_silent(sf <- newFit(SE, K = 2, epsilon = 1e4))
    expect_silent(sf <- newFit(SE, K = 2, epsilon = 1e5))
    expect_silent(sf <- newFit(SE, K = 2, epsilon = 1e12))


})
