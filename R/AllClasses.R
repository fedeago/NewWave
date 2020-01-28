#' Class ZinbModel
#'
#' Objects of this class store all the values needed to work with a
#' negative binomial model, as described in the vignette.
#' They contain all information to fit a model by penalized maximum likelihood
#' or simulate data from a model.
#'
#' @slot X matrix. The design matrix containing sample-level covariates, one
#'   sample per row.
#' @slot V matrix. The design matrix containing gene-level covariates, one gene
#'   per row.
#' @slot O matrix. The offset matrix.
#' @slot X_intercept logical. TRUE if X contains an intercept.
#' @slot V_intercept logical. TRUE if V contains an intercept.
#' @slot W matrix. The factors of sample-level latent factors.
#' @slot beta matrix or NULL. The coefficients of X in the regression.
#' @slot gamma matrix or NULL. The coefficients of V in the regression.
#' @slot alpha  matrix. The weight of sample-level latent factors.
#' @slot zeta numeric. A vector of log of inverse dispersion parameters.
#' @slot epsilon_beta nonnegative scalar. Regularization parameter for
#'   beta
#' @slot epsilon_gamma nonnegative scalar. Regularization parameter for
#'   gamma
#' @slot epsilon_W nonnegative scalar. Regularization parameter for W
#' @slot epsilon_alpha nonnegative scalar. Regularization parameter for alpha
#' @slot epsilon_zeta nonnegative scalar. Regularization parameter for zeta
#'
#'
#' @details For the full description of the model see the model vignette.
#'   Internally, the slots are checked so that the matrices are of the
#'   appropriate dimensions: in particular, \code{X}, \code{O}
#'   and \code{W} need to have \code{n} rows, \code{V} needs to have \code{J}
#'   rows, \code{zeta} must be of length \code{J}.
#' @name ZinbModel-class
#' @import methods
#' @exportClass ZinbModel
#' @aliases ZinbModel
#'
#' @return \code{nSamples} returns the number of samples; \code{nFeatures}
#' returns the number of features; \code{nFactors} returns the number of latent
#' factors.
#'
setClass(
    Class = "nbModel",
    slots = list(X = "matrix",
                 V = "matrix",
                 X_intercept = "logical",
                 V_intercept = "logical",
                 W = "matrix",
                 beta = "matrix",
                 gamma = "matrix",
                 alpha = "matrix",
                 zeta = "numeric",
                 epsilon_beta = "numeric",
                 epsilon_gamma = "numeric",
                 epsilon_W = "numeric",
                 epsilon_alpha = "numeric",
                 epsilon_zeta = "numeric")
)

setValidity("nbModel", function(object){
    n <- NROW(getX(object)) # number of samples
    J <- NROW(getV(object)) # number of genes
    K <- NCOL(getW(object)) # number of latent factors

    if(K > n) {
        return("Cannot have more latent factors than samples.")
    }
    if(K > J) {
        return("Cannot have more latent factors than genes.")
    }
    if(NROW(getW(object)) != n) {
        return("W must have n rows!")
    }
    if(NCOL(getX(object)) != NROW(getBeta(object))){
        return("beta must have the same number of rows as there are columns in X!")
    }
    if(NCOL(getV(object)) != NROW(getGamma(object))){
        return("gamma must have the same number of rows as there are columns in V!")
    }
    if(NCOL(getBeta(object)) != J) {
        return("beta must have J columns!")
    }
    if(NCOL(getGamma(object)) != n) {
        return("gamma must have n columns!")
    }
    if(NCOL(getAlpha(object)) != J) {
        return("alpha must have J columns!")
    }
    if(NROW(getAlpha(object)) != K) {
        return("alpha must have K rows!")
    }
    if(length(getZeta(object)) != J) {
        return("zeta must have length J!")
    }
    if((length(object@epsilon_beta) != 1) || (object@epsilon_beta < 0)) {
        return("epsilon_beta must be a nonnegative scalar !")
    }
    if((length(object@epsilon_gamma) != 1) || (object@epsilon_gamma < 0)) {
        return("epsilon_gamma must be a nonnegative scalar !")
    }
    if((length(object@epsilon_W) != 1) || (object@epsilon_W < 0)) {
        return("epsilon_W must be a nonnegative scalar !")
    }
    if((length(object@epsilon_alpha) != 1) || (object@epsilon_alpha < 0)) {
        return("epsilon_alpha must be a nonnegative scalar !")
    }
    if((length(object@epsilon_zeta) != 1) || (object@epsilon_zeta < 0)) {
        return("epsilon_zeta must be a nonnegative scalar !")
    }

    return(TRUE)
}
)
