#' Generic function that returns the number of latent factors
#'
#' Given an object that describes a dataset or a model involving latent factors,
#' this function returns the number of latent factors.
#' @param x an object that describes a dataset or a model involving latent
#'  factors
#' @return the number of latent factors
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberFactors(a)
#' @export
setGeneric("numberFactors", function(x) standardGeneric("numberFactors"))


#' Generic function that returns the number of samples
#'
#' Given an object that describes a dataset or a model involving samples,
#' this function returns the number of samples.
#' @param x an object that describes a dataset or a model involving samples
#' @return the number of samples
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberSamples(a)
#' @export
setGeneric("numberSamples", function(x) standardGeneric("numberSamples"))

#' Generic function that returns the number of features
#'
#' Given an object that describes a dataset or a model involving features,
#' this function returns the number of features
#' @param x an object that describes a dataset or a model involving features
#' @return the number of features
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberFeatures(a)
#' @export
setGeneric("numberFeatures", function(x) standardGeneric("numberFeatures"))

#' Generic function that returns the number of params
#'
#' Given an object that describes a dataset or a model involving params,
#' this function returns the number of params
#' @param x an object that describes a dataset or a model involving params
#' @return the number of params
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberParams(a)
#' @export
setGeneric("numberParams", function(x) standardGeneric("numberParams"))

#' Returns the matrix of mean parameters
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the matrix of mean parameters.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the matrix of mean parameters
#' @details Note that although the user interface of \code{\link{newFit}}
#' requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#' samples in row and genes in column). Hence the parameter matrix returned by
#' this function is of n x J dimensions.
#' @examples
#' a <- newmodel(n=5, J=10)
#' newMu(a)
#' @export
setGeneric("newMu", function(object) standardGeneric("newMu"))

#' Returns the matrix of logarithm of mean parameters
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the matrix of logarithm of mean parameters.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the matrix of logarithms of mean parameters
#' @details Note that although the user interface of \code{\link{newFit}}
#'  requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#'  samples in row and genes in column). Hence the parameter matrix returned by
#'  this function is of n x J dimensions.
#' @examples
#' a <- newmodel(n=5, J=10)
#' newLogMu(a)
#' @export
setGeneric("newLogMu", function(object) standardGeneric("newLogMu"))

#' Returns the vector of dispersion parameters
#'
#' Given an object that describes a matrix of negative binomial negative binomial
#' distributions, returns the vector of dispersion parameters \code{phi}.
#' @param object an object that describes a matrix of negative binomial.
#'  distributions.
#' @return the vector of dispersion parameters
#' @examples
#' a <- newmodel(n=5, J=10)
#' newPhi(a)
#' @export
setGeneric("newPhi", function(object) standardGeneric("newPhi"))

#' Returns the vector of inverse dispersion parameters
#'
#' Given an object that describes a matrix of negative binomial negative binomial
#' distributions, returns the vector of inverse dispersion parameters
#' \code{theta}.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the vector of inverse dispersion parameters theta
#' @examples
#' a <- newmodel(n=5, J=10)
#' newTheta(a)
#' @export
setGeneric("newTheta", function(object) standardGeneric("newTheta"))

#' Returns the vector of log of inverse dispersion parameters
#'
#' Given an object that describes a matrix of negative binomial negative binomial
#' distributions, returns the vector \code{zeta} of log of inverse dispersion
#' parameters
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the vector \code{zeta} of log of inverse dispersion parameters
#' @examples
#' a <- newmodel(n=5, J=10)
#' newZeta(a)
#' @export
setGeneric("newZeta", function(object) standardGeneric("newZeta"))

#' Returns the low-dimensional matrix of inferred sample-level covariates W
#'
#' Given an object that contains the fit of a nb-WaVE model, returns the
#' matrix \code{W} of low-dimensional matrix of inferred sample-level
#' covariates.
#'
#' @param object a \code{\linkS4class{newmodel}} object, typically the result
#'  of \code{\link{newFit}}.
#' @return the matrix \code{W} of inferred sample-level covariates.
#' @examples
#' a <- newmodel(n=5, J=10)
#' newW(a)
#' @export
setGeneric("newW", function(object) standardGeneric("newW"))

#' Simulate counts from a negative binomial model
#'
#' Given an object that describes negative binomial distribution,
#' simulate counts from the distribution.
#' @param object an object that describes a matrix of  negative
#'  binomial.
#' @param seed an optional integer to specify how the random number generator
#'  should be initialized with a call to \code{set.seed}. If missing, the
#'  random generator state is not changed.
#' @param ... additional arguments.
#' @return A list with the following elements.
#'  \itemize{
#'  \item{counts}{the matrix with the simulated counts.}
#'  \item{dataNB}{the data simulated from the negative binomial.}
#'  \item{dataDropouts}{the data simulated from the binomial process.}
#'  \item{zeroFraction}{the fraction of zeros.}
#'  }
#' @examples
#' a <- newmodel(n=5, J=10)
#' newSim(a)
#' @export
setGeneric("newSim",function(object, seed, ...) standardGeneric("newSim"))



#' Compute the penalty of a model
#'
#' Given a statistical model with regularization parameters, compute the
#' penalty.
#' @param model an object that describes a statistical model with regularization
#'  parameters.
#' @return The penalty of the model.
#' @examples
#' m <- newmodel(K=2)
#' newpenalty(m)
#' @export
setGeneric("newpenalty", function(model) standardGeneric("newpenalty"))

#' Fit a nb regression model
#'
#' Given an object with the data, it fits a nb model.
#'
#' @param Y The data (genes in rows, samples in columns).
#' @param ... Additional parameters to describe the model, see
#'  \code{\link{newmodel}}.
#' @return An object of class \code{newmodel} that has been fitted by penalized
#'  maximum likelihood on the data.
#' @export
setGeneric("newFit", function(Y, ...) standardGeneric("newFit"))

#' Returns the sample-level design matrix for mu
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the sample-level design matrix for mu
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @param ... Additional parameters.
#' @return the sample-level design matrix for mu
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newX(a)
setGeneric("newX", function(object, ...) standardGeneric("newX"))


#' Returns the gene-level design matrix for mu
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the gene-level design matrix for mu
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @param ... Additional parameters.
#' @return the gene-level design matrix for mu
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newV(a)
setGeneric("newV", function(object, ...) standardGeneric("newV"))


#' Returns the matrix of paramters beta
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the matrix of parameters associated with X
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @param ... Additional parameters.
#' @return the matrix of beta parameters
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newBeta(a)
setGeneric("newBeta", function(object, ...) standardGeneric("newBeta"))


#' Returns the matrix of paramters gamma
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the matrix of parameters associated with V
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @param ... Additional parameters.
#' @return the matrix of gamma parameters
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newGamma(a)
setGeneric("newGamma", function(object, ...) standardGeneric("newGamma"))


#' Returns the matrix of paramters alpha
#'
#' Given an object that describes a matrix of negative binomial distributions,
#' returns the matrix of parameters associated with W for the mean part (mu)
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @param ... Additional parameters.
#' @return the matrix of alpha parameters
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newAlpha(a)
setGeneric("newAlpha", function(object, ...) standardGeneric("newAlpha"))


#' Returns the vector of regularization parameter for beta
#'
#' Given an object describing a nb model, returns a vector of size the number
#' of rows in the parameter \code{beta} with the regularization parameters
#' associated to each row.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the regularization parameters for \code{beta}.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newEpsilon_beta(a)
setGeneric("newEpsilon_beta",
           function(object) standardGeneric("newEpsilon_beta"))

#' Returns the vector of regularization parameter for gamma
#'
#' Given an object describing a nb model, returns a vector of size the number
#' of columns in the parameter \code{gamma} with the regularization
#' parameters associated to each row.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the regularization parameters for \code{gamma}.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newEpsilon_gamma(a)
setGeneric("newEpsilon_gamma",
           function(object) standardGeneric("newEpsilon_gamma"))



#' Returns the vector of regularization parameter for W
#'
#' Given an object describing a nb model, returns a vector of size the number
#' of columns in the parameter \code{W} with the regularization
#' parameters associated to each column.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the regularization parameters for \code{W}.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newEpsilon_W(a)
setGeneric("newEpsilon_W", function(object) standardGeneric("newEpsilon_W"))

#' Returns the vector of regularization parameter for alpha
#'
#' Given an object describing a nb model, returns a vector of size the number
#' of rows in the parameter \code{alpha} with the regularization parameters
#' associated to each row.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the regularization parameters for \code{alpha}.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newEpsilon_alpha(a)
setGeneric("newEpsilon_alpha",
           function(object) standardGeneric("newEpsilon_alpha"))

#' Returns the regularization parameter for the dispersion parameter
#'
#' The regularization parameter penalizes the variance of zeta, the log of
#' the dispersion parameters across samples.
#' @param object an object that describes a matrix of negative binomial
#'  distributions.
#' @return the regularization parameters for \code{zeta}.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' newEpsilon_zeta(a)
setGeneric("newEpsilon_zeta",
           function(object) standardGeneric("newEpsilon_zeta"))

#' Generic function that returns the number of features
#'
#' Given an object that describes a dataset or a model, it returns the number of
#' features.
#' @param x an object that describes a dataset or a model.
#' @return the number of features.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberFeatures(a)
setGeneric(
    name = "numberFeatures",
    def = function(x) {
        standardGeneric("numberFeatures")
    }
)

#' Generic function that returns the number of samples
#'
#' Given an object that describes a model or a dataset, it returns the number of
#' samples.
#' @param x an object that describes a dataset or a model.
#' @return the number of samples.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberSamples(a)
setGeneric(
    name = "numberSamples",
    def = function(x) {
        standardGeneric("numberSamples")
    }
)

#' Generic function that returns the total number of parameters of the model
#'
#' Given an object that describes a model or a dataset, it returns total number
#'  of parameters of the model.
#' @param model an object that describes a dataset or a model.
#' @return the total number of parameters of the model.
#' @export
#' @examples
#' a <- newmodel(n=5, J=10)
#' numberParams(a)
setGeneric(
    name = "numberParams",
    def = function(model) {
        standardGeneric("numberParams")
    }
)

#' Perform dimensionality reduction using a nb regression model with
#' gene and cell-level covariates.
#'
#' Given an object with the data, it performs dimensionality reduction using
#' a nb regression model with gene and cell-level covariates.
#'
#' @param Y The data (genes in rows, samples in columns). Currently implemented
#'  only for \code{SummarizedExperiment}.
#' @param ... Additional parameters to describe the model, see
#'  \code{\link{newmodel}}.
#' @return An object of class \code{SingleCellExperiment}; the dimensionality
#'  reduced matrix is stored in the \code{reducedDims} slot and optionally
#'  normalized values and residuals are added in the list of assays.
#' @export
setGeneric("newWave", function(Y, ...) standardGeneric("newWave"))

#' Compute the log-likelihood of a model given some data
#'
#' Given a statistical model and some data, this function computes the
#' log-likelihood of the model given the data, i.e., the log-probability of the
#' data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @param ... additional arguments.
#' @return The log-likelihood of the model given the data.
#' @examples
#' m <- newmodel(n=5, J=10)
#' x <- newSim(m)
#' newloglik(m, x$counts)
#' @export
setGeneric("newloglik", function(model, x, ...) standardGeneric("newloglik"))

#' Compute the AIC of a model given some data
#'
#' Given a statistical model and some data, this function computes the AIC
#' of the model given the data, i.e., the AIC of the data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @return the AIC of the model.
#' @examples
#' m <- newmodel(n=5, J=10)
#' x <- newSim(m)
#' newAIC(m, x$counts)
#' @export
setGeneric(
    name = "newAIC",
    def = function(model, x) {
        standardGeneric("newAIC")
    }
)

#' Compute the BIC of a model given some data
#'
#' Given a statistical model and some data, this function computes the BIC
#' of the model given the data, i.e., the BIC of the data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @return the BIC of the model.
#' @examples
#' m <- newmodel(n=5, J=10)
#' x <- newSim(m)
#' newBIC(m, x$counts)
#' @export
setGeneric(
    name = "newBIC",
    def = function(model, x) {
        standardGeneric("newBIC")
    }
)


