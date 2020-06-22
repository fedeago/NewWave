#' Log-likelihood of the zero-inflated negative binomial model for each entry
#' in the matrix of counts
#'
#' Given a matrix of counts, this function computes the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (NB) model. For each count, the NB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param model the nb model
#' @param x the matrix of counts
#' @return the matrix of log-likelihood of the model.
#' @importFrom stats dnbinom
nb.loglik.matrix <- function(model, x) {
  mu <- getMu(model)
  theta <- getTheta(model)
  theta_mat <- matrix(rep(theta, each = nrow(x)), ncol = ncol(x))
  pi <- getPi(model)
  lik <- pi * (x == 0) + (1 - pi) * dnbinom(x, size = theta_mat, mu = mu)
  lik[lik == 0] <- min(lik[lik != 0]) #to avoid log lik to be infinite
  log(lik)
}

#' Observational weights of the zero-inflated negative binomial model for each entry
#' in the matrix of counts
#'
#' Given a matrix of counts, this function computes the
#' observational weights of the counts under a zero-inflated negative binomial
#' (NB) model. For each count, the NB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param model the nb model
#' @param x the matrix of counts
#' @return the matrix of observational weights computed from the model.
#' @importFrom stats dnbinom
#' @export
computeObservationalWeights <- function(model, x){
  mu <- getMu(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  nbwg <- t(nbwg)
  nbwg[x > 0] <- 1
  nbwg[nbwg < 1e-15] <- 1e-15
  nbwg
}



#' Deviance residuals of negative binomial model
#'
#' Given a matrix of counts, this function computes the
#' deviance residuals under a negative binomial
#' (NB) model.
#'
#' @param model the nb model
#' @param x the matrix of counts n cells by J genes
#' @param ignoreW logical, if true matrix \code{W} is ignored. Default is TRUE.
#' @export
#' @return the matrix of deviance residuals of the model.
computeDevianceResiduals <- function(model, x, ignoreW = TRUE) {
  
  # this makes a copy of "model" -- is there a more efficient way?
  if (ignoreW) {
    model@W <- matrix(0, ncol = nFactors(model), nrow = nSamples(model))
  }
  
  mu_hat <- getMu(model)
  ll <- nb.loglik.matrix(model, x)
  sign <- 1*(x - x_hat > 0)
  sign[sign == 0] <- -1
  sign * sqrt(-2 * ll)
}


CleanEnvir <- function() {
  objs <- ls(pos = ".GlobalEnv")
  rm(list =c(Y_sh,X_sh,V_sh,mu,beta_sh, alpha_sh, gamma_sh, W_sh, zeta_sh), pos = ".GlobalEnv")
}

#' @describeIn nbwave Y is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @export
#'
#' @param X The design matrix containing sample-level covariates, one sample per
#'   row. If missing, X will contain only an intercept. If Y is a
#'   SummarizedExperiment object, X can be a formula using the variables in the
#'   colData slot of Y.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept. If Y is a
#'   SummarizedExperiment object, V can be a formula using the variables in the
#'   rowData slot of Y.
#' @param K integer. Number of latent factors. Specify \code{K = 0} if only
#'   computing observational weights.
#' @param fitted_model a \code{\link{nbModel}} object.
#' @param which_assay numeric or character. Which assay of Y to use. If missing,
#'   if `assayNames(Y)` contains "counts" then that is used. Otherwise, the
#'   first assay is used.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose Print helpful messages.
#' @param nb.repeat.initialize Number of iterations for the initialization of
#'   beta_mu and gamma_mu.
#' @param maxiter.optimize maximum number of iterations for the optimization
#'   step (default 25).
#' @param stop.epsilon.optimize stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#' @param normalizedValues indicates wether or not you want to compute
#' normalized values for the counts after adjusting for gene and cell-level
#' covariates.
#' @param residuals indicates wether or not you want to compute the residuals
#' of the NB model. Deviance residuals are computed.
#' @param observationalWeights indicates whether to compute the observational
#'   weights for differential expression (see vignette).
#'
#' @details For visualization (heatmaps, ...), please use the normalized values.
#' It corresponds to the deviance residuals when the \code{W} is not included
#' in the model but the gene and cell-level covariates are. As a results, when
#' \code{W} is not included in the model, the deviance residuals should capture
#' the biology. Note that we do not recommend to use the normalized values for
#' any downstream analysis (such as clustering, or differential expression), but
#' only for visualization.
#'
#' @details If one has already fitted a model using \code{\link{nbModel}},
#' the object containing such model can be used as input of \code{nbwave} to
#' save the resulting W into a \code{SummarizedExperiment} and optionally
#' compute residuals and normalized values, without the need for re-fitting the
#' model.
#'
#' @details By default \code{nbwave} uses all genes to estimate \code{W}.
#'   However, we recommend to use the top 1,000 most variable genes for this
#'   step. In general, a user can specify any custom set of genes to be used to
#'   estimate \code{W}, by specifying either a vector of gene names, or a single
#'   character string corresponding to a column of the \code{rowData}.
#'
#' @details Note that if both \code{which_genes} is specified and at least one
#'   among \code{observationalWeights}, \code{imputedValues}, \code{residuals},
#'   and \code{normalizedValues} is \code{TRUE}, the model needs to be fit
#'   twice.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- nbwave(se, X="~bio")
setMethod("nbwave", "SummarizedExperiment",
          function(Y, X, V, K=2, which_assay,
                   commondispersion = TRUE, verbose=FALSE,
                   maxiter_optimize=100,
                   stop_epsilon=.0001,  children = 1,
                   random_init = FALSE, random_start = FALSE,
                   n_gene_disp = NULL,
                   n_cell_par = NULL, n_gene_par = NULL,
                   normalizedValues = FALSE, residuals = FALSE,
                   observationalWeights = FALSE, cross_batch = FALSE,...) {
              
            
            
              fitted_model <- nbFit(Y, X, V, K,
                                      which_assay, commondispersion,
                                      verbose,
                                      maxiter_optimize, stop_epsilon,
                                      children, random_init, 
                                      random_start, n_gene_disp,
                                      n_cell_par, n_gene_par, cross_batch, ...)
              
              
              
              out <- as(Y, "SingleCellExperiment")
              
              if (nFactors(fitted_model) > 0){
                W <- getW(fitted_model)
                colnames(W) <- paste0('W', seq_len(nFactors(fitted_model)))
                reducedDim(out, "nbwave") <- W
              }
              
              if(missing(which_assay)) {
                if("counts" %in% assayNames(Y)) {
                  dataY <- assay(Y, "counts")
                } else {
                  warning("No assay named `counts`, using first assay.",
                          "Use `assay` to specify a different assay.")
                  dataY <- assay(Y)
                }
              } else {
                if(!(is.character(which_assay) | is.numeric(which_assay))) {
                  stop("assay needs to be a numeric or character specifying which assay to use")
                } else {
                  dataY <- assay(Y, which_assay)
                }
              }
              
              refit <- any(c(observationalWeights, normalizedValues, residuals)) &
                (nrow(Y) != nFeatures(fitted_model)) & K > 0
              
              if(refit) {
                fitted_model <- nbFit(Y, X, V, K, which_assay,
                                      commondispersion, verbose,
                                      maxiter_optimize,
                                      stop_epsilon, children,
                                      random_init, random_start,
                                      n_gene_disp,
                                      n_cell_par, n_gene_par, cross_batch,...)
              }
              
              if (normalizedValues){
                norm <- computeDevianceResiduals(fitted_model, t(dataY),
                                                 ignoreW = TRUE)
                assay(out, "normalizedValues") <- t(norm)
              }
              
              if (residuals){
                devres <- computeDevianceResiduals(fitted_model, t(dataY),
                                                   ignoreW = FALSE)
                assay(out, "residuals") <- t(devres)
              }
              
              
              if (observationalWeights) {
                weights <- computeObservationalWeights(fitted_model, dataY)
                dimnames(weights) <- dimnames(out)
                assay(out, "weights") <- weights
              }
              
              #CleanEnvir()
              
              return(out)
          }
)
