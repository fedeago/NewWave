#' Initialize an object of class newmodel
#' 
#' @param X matrix. The design matrix containing sample-level covariates, one
#'   sample per row.
#' @param V matrix. The design matrix containing gene-level covariates, one gene
#'   per row.
#' @param W matrix. The factors of sample-level latent factors.
#' @param beta matrix or NULL. The coefficients of X in the regression of mu.
#' @param gamma matrix or NULL. The coefficients of V in the regression of
#'   mu.
#' @param alpha matrix or NULL. The coefficients of W in the regression of
#'   mu.
#' @param zeta numeric. A vector of log of inverse dispersion parameters.
#' @param epsilon nonnegative scalar. Regularization parameter.
#' @param epsilon_beta nonnegative scalar. Regularization parameter for
#'   beta.
#' @param epsilon_gamma nonnegative scalar. Regularization parameter for
#'   gamma.
#' @param epsilon_W nonnegative scalar. Regularization parameter for W.
#' @param epsilon_alpha nonnegative scalar. Regularization parameter for alpha
#' @param epsilon_zeta nonnegative scalar. Regularization parameter for zeta.
#' @param n integer. Number of samples.
#' @param J integer. Number of genes.
#' @param K integer. Number of latent factors.
#' 
#' @export
#' 
#' 
#' @details This is a wrapper around the new() function to create an
#'   instance of class \code{newmodel}. Rarely, the user will need to create a
#'   \code{newmodel} object from scratch, as tipically this is the result of
#'   \code{\link{newFit}}.
#' 
#' @details If any of \code{X}, \code{V}, \code{W} matrices are passed,
#'   \code{n}, \code{J}, and \code{K} are inferred. Alternatively, the user can
#'   specify one or more of \code{n}, \code{J}, and \code{K}.
#' 
#' @details The regularization parameters can be set by a unique parameter
#'   \code{epsilon} or specific values for the different regularization
#'   parameters can also be provided.
#'   If only \code{epsilon} is specified, the other parameters take the
#'   following values:
#'   \itemize{
#'   \item epsilon_beta = epsilon/J
#'   \item epsilon_gamma = epsilon/n
#'   \item epsilon_W = epsilon/n
#'   \item epsilon_alpha = epsilon/J
#'   \item epsilon_zeta = epsilon
#'   }
#'   We empirically found that large values of \code{epsilon} provide a more
#'   stable estimation of \code{W}.
#' 
#' @details A call with no argument has the following default values: \code{n =
#'   50}, \code{J = 100}, \code{K = 0}, \code{epsilon=J}.
#' 
#' @details Although it is possible to create new instances of the class by
#'   calling this function, this is not the most common way of creating
#'   \code{newmodel} objects. The main use of the class is within the
#'   \code{\link{newFit}} function.
#' 
#' @return an object of class \code{\linkS4class{newmodel}}.
#' 
#' @examples
#' a <- newmodel()
#' numberSamples(a)
#' numberFeatures(a)
#' numberFactors(a)

newmodel <- function(X, V, W, beta,
                       gamma,alpha, zeta, epsilon,
                       epsilon_beta, epsilon_gamma, epsilon_W, epsilon_alpha,
                       epsilon_zeta, n, J, K) {

   # Find n (default 50), J (default 100), K (default 0)
   if (missing(n)) {
           if (!missing(X)) {
               n <- NROW(X)
           } else if (!missing(gamma)) {
               n <- NCOL(gamma)
           }  else if (!missing(W)) {
               n <- NROW(W)
           } else {
               n <- 50
           }
       }
           if (missing(J)) {
               if (!missing(V)) {
               J <- NROW(V)
           } else if (!missing(beta)) {
               J <- NCOL(beta)
           }  else if (!missing(alpha)) {
               J <- NCOL(alpha)
           }  else if (!missing(zeta)) {
               J <- length(zeta)
           } else {
               J <- 100
           }
       }
       if (missing(K)) {
           if (!missing(W)) {
               K <- NCOL(W)
           } else {
               K <- 0
           }
       }

    # Set the different slots for the matrices
       if(missing(X)) {
           X <- matrix(1, nrow=n, ncol=1)
       }

       if (missing(V)) {
           V <- matrix(1, nrow=J, ncol=1)
       }


       X_intercept <- FALSE
       if(ncol(X) > 0) {
           X_intercept <- TRUE
       }


       V_intercept <- FALSE
       if(ncol(V) > 0) {
           V_intercept <- TRUE

       }

       if (missing(W)) {
           W <- matrix(0, nrow=n , ncol=K)
       }

       if (missing(beta)) {
           beta <- matrix(0, nrow=ncol(X), ncol=J)
       }
       if (missing(gamma)) {
           gamma <- matrix(0, nrow=ncol(V), ncol=n)
       }
       if (missing(alpha)) {
           alpha <- matrix(0, nrow=K , ncol=J)
       }
       if (missing(zeta)) {
           zeta <- numeric(J)
       }

       # Regularization parameters
       if (missing(epsilon)) {
           epsilon <- J
       }
       if (missing(epsilon_beta)) {
           epsilon_beta <- epsilon/J
       }
       if (missing(epsilon_gamma)) {
           epsilon_gamma <- epsilon/n
       }
       if (missing(epsilon_W)) {
           epsilon_W <- epsilon/n
       }
       if (missing(epsilon_alpha)) {
           epsilon_alpha <- epsilon/J
       }
       if (missing(epsilon_zeta)) {
           epsilon_zeta <- epsilon
       }

       obj <- new(Class="newmodel",
               X = X, V = V, X_intercept = X_intercept,
               V_intercept = V_intercept, W = W, beta = beta,
               gamma = gamma, alpha = alpha, zeta = zeta,
               epsilon_beta = epsilon_beta,
               epsilon_gamma = epsilon_gamma,
               epsilon_W = epsilon_W, epsilon_alpha = epsilon_alpha,
               epsilon_zeta = epsilon_zeta)

       validObject(obj) # call of the inspector
       return(obj)
}

#' @export
#' @describeIn newmodel show useful info on the object.
#'
#' @param object an object of class \code{newmodel}.
setMethod("show", "newmodel",
          function(object) {
              cat(paste0("Object of class newmodel.\n",
                numberSamples(object), " samples; ", numberFeatures(object),
                         " genes.\n",
                         NCOL(newX(object)),
                         " sample-level covariate(s) (mu); ",
                         NCOL(newV(object)),
                         " gene-level covariate(s) (mu); ",
                         numberFactors(object), " latent factor(s).\n"))
          }
)


################################################################
# Extract various informations and variables from a newmodel #
################################################################

#' @export
#' @describeIn newmodel returns the number of samples.
#' @param x an object of class \code{newmodel}.
setMethod("numberSamples", "newmodel",
          function(x) {
              return(NROW(x@X))
          }
)

#' @export
#' @describeIn newmodel returns the number of features.
setMethod("numberFeatures", "newmodel",
          function(x) {
              return(NROW(x@V))
          }
)

#' @export
#' @describeIn newmodel returns the number of latent factors.
setMethod("numberFactors", "newmodel",
          function(x) {
              return(NCOL(x@W))
          }
)


#' @export
#' @describeIn newmodel returns the sample-level design matrix for mu.
setMethod("newX", "newmodel",
          function(object) {
              return(object@X)
          }
)

#' @export
#' @describeIn newmodel returns the gene-level design matrix for mu.
setMethod("newV", "newmodel",
          function(object) {
              return(object@V)
          }
)

#' @export
#' @describeIn newmodel returns the logarithm of the mean of the non-zero
#'   component.
setMethod("newLogMu", "newmodel",
          function(object) {
              return(newX(object) %*% object@beta +
                         t(newV(object) %*% object@gamma) +
                         object@W %*% object@alpha)
          }
)

#' @export
#' @describeIn newmodel returns the mean of the non-zero component.
setMethod("newMu", "newmodel",
    function(object) {
        return(exp(newLogMu(object)))
    }
)

#' @export
#' @describeIn newmodel returns the log of the inverse of the dispersion
#'   parameter.
setMethod("newZeta", "newmodel",
          function(object) {
              return(object@zeta)
          }
)

#' @export
#' @describeIn newmodel returns the dispersion parameter.
setMethod("newPhi", "newmodel",
          function(object) {
              return(exp(-object@zeta))
          }
)

#' @export
#' @describeIn newmodel returns the inverse of the dispersion parameter.
setMethod("newTheta", "newmodel",
          function(object) {
              return(exp(object@zeta))
          }
)

#' @export
#' @describeIn newmodel returns the regularization parameters for
#'   \code{beta}.
setMethod("newEpsilon_beta", "newmodel",
          function(object) {
              e <- rep(object@epsilon_beta, ncol(object@X))
              if (object@X_intercept) {
                  e[1] <- 0
              }
              e
          }
)

#' @export
#' @describeIn newmodel returns the regularization parameters for
#'   \code{gamma}.
setMethod("newEpsilon_gamma", "newmodel",
          function(object) {
              e <- rep(object@epsilon_gamma, ncol(object@V))
              if (object@V_intercept) {
                  e[1] <- 0
              }
              e
          }
)



#' @export
#' @describeIn newmodel returns the regularization parameters for
#'   \code{W}.
setMethod("newEpsilon_W", "newmodel",
          function(object) {
              rep(object@epsilon_W, numberFactors(object))
          }
)

#' @export
#' @describeIn newmodel returns the regularization parameters for
#'   \code{alpha}.
setMethod("newEpsilon_alpha", "newmodel",
          function(object) {
              rep(object@epsilon_alpha, numberFactors(object))
          }
)

#' @export
#' @describeIn newmodel returns the regularization parameters for
#'   \code{zeta}.
setMethod("newEpsilon_zeta", "newmodel",
          function(object) {
              object@epsilon_zeta
          }
)

#' @export
#' @describeIn newmodel returns the matrix W of inferred sample-level
#'   covariates.
setMethod("newW", "newmodel",
          function(object) {
              object@W
          }
)

#' @export
#' @describeIn newmodel returns the matrix beta of inferred parameters.
setMethod("newBeta", "newmodel",
          function(object) {
              object@beta
          }
)

#' @export
#' @describeIn newmodel returns the matrix gamma of inferred parameters.
setMethod("newGamma", "newmodel",
          function(object) {
              object@gamma
          }
)


#' @export
#' @describeIn newmodel returns the matrix alpha of inferred parameters.
setMethod("newAlpha", "newmodel",
          function(object) {
              object@alpha
          }
)


#' @export
#' @describeIn numberParams returns the total number of parameters in the model.
setMethod("numberParams", "newmodel",
          function(model) {

              X <- newX(model)

              V <- newV(model)

              n <- numberSamples(model)
              J <- numberFeatures(model)
              K <- numberFactors(model)

              M <- NCOL(X)

              L <- NCOL(V)
              ndisp <- length(unique(newZeta(model)))

              J * (M) + n * (L) + 2 * K * J + n * K + ndisp
          }
)

########################
# Other useful methods #
########################

#' @export
#' @describeIn newSim simulate from a nb distribution.
#'
setMethod(
    f="newSim",
    signature="newmodel",
    definition=function(object, seed) {

        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            runif(1)
        }

        if (missing(seed)) {
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        } else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }

        mu <- newMu(object)
        theta <- newTheta(object)
        n <- numberSamples(object)
        J <- numberFeatures(object)

        # Simulate negative binomial with the mean matrix and dispersion
        # parameters
        i <- seq(n*J)
        datanb <- rnbinom(length(i), mu = mu[i], size = theta[ceiling(i/n)])
        data.nb <- matrix(datanb, nrow = n)


        # Matrix of counts
        counts <- data.nb

        # Fraction of zeros in the matrix
        zero.fraction <- sum(counts == 0) / (n*J)

        ret <- list(counts = t(counts),
                    zeroFraction = zero.fraction)
        attr(ret, "seed") <- RNGstate
        ret
    }
)

#' @export
#' @describeIn newloglik return the log-likelihood of the nb model.
setMethod(
    f="newloglik",
    signature=c("newmodel","matrix"),
    definition=function(model, x) {
        nb.loglik(t(x), newMu(model),
                  rep(newTheta(model), rep(ncol(x),nrow(x))))
    }
)

#' @export
#' @describeIn newAIC returns the AIC of the NB model.
setMethod(
    f="newAIC",
    signature=c("newmodel","matrix"),
    definition=function(model, x) {
    if ((numberSamples(model) != ncol(x))|(numberFeatures(model) != nrow(x))) {
            stop("x and model should have the same dimensions!")
        }
        k <- numberParams(model)
        ll <- newloglik(model, x)
        return(2*k - 2*ll)
    }
)

#' @export
#' @describeIn newBIC returns the BIC of the NB model.
setMethod(
    f="newBIC",
    signature=c("newmodel","matrix"),
    definition=function(model, x) {
        n <- numberSamples(model)
        if ((n != ncol(x))|(numberFeatures(model) != nrow(x))) {
            stop("x and model should have the same dimensions!")
        }
        k <- numberParams(model)
        ll <- newloglik(model, x)
        return(log(n)*k - 2*ll)
    }
)

#' @export
#' @describeIn newpenalty return the penalization.
setMethod(
    f="newpenalty",
    signature="newmodel",
    definition=function(model) {
        sum(newEpsilon_alpha(model)*(model@alpha)^2)/2 +
        sum(newEpsilon_beta(model)*(model@beta)^2)/2 +
        sum(newEpsilon_gamma(model)*(model@gamma)^2)/2 +
        sum(newEpsilon_W(model)*t(model@W)^2)/2 +
        newEpsilon_zeta(model)*var(newZeta(model))/2
    }
)

# Copied on 5/14/2019 from log1pexp of the copula package (v. 0.999.19) by
# Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan,
# Johanna G. Neslehova
# Copied here to avoid dependence on gsl which causes troubles.
log1pexp <- function (x, c0 = -37, c1 = 18, c2 = 33.3)
{
    if (has.na <- any(ina <- is.na(x))) {
        y <- x
        x <- x[ok <- !ina]
    }
    r <- exp(x)
    if (any(i <- c0 < x & (i1 <- x <= c1)))
        r[i] <- log1p(r[i])
    if (any(i <- !i1 & (i2 <- x <= c2)))
        r[i] <- x[i] + 1/r[i]
    if (any(i3 <- !i2))
        r[i3] <- x[i3]
    if (has.na) {
        y[ok] <- r
        y
    }
    else r
}