#' Log-likelihood of the zero-inflated negative binomial model for each entry
#' in the matrix of counts
#'
#' Given a matrix of counts, this function computes the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (NB) model. For each count, the NB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param model the newmodel
#' @param x the matrix of counts
#' @return the matrix of log-likelihood of the model.
#' @importFrom stats dnbinom
nb.loglik.matrix <- function(model, x) {
    mu <- newMu(model)
    theta <- newTheta(model)
    theta_mat <- matrix(rep(theta, each = nrow(x)), ncol = ncol(x))
    lik <- dnbinom(x, size = theta_mat, mu = mu)
    lik[lik == 0] <- min(lik[lik != 0]) #to avoid log lik to be infinite
    log(lik)
}



#' @describeIn newWave Y is a
#'   \code{SummarizedExperiment}.
#' @export
#' 
#' @param Y The SummarizedExperiment with the data
#' @param X The design matrix containing sample-level covariates, one sample 
#'   per row. If missing, X will contain only an intercept. If Y is a
#'   SummarizedExperiment object, X can be a formula using the variables in the
#'   colData slot of Y.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept. If Y is a
#'   SummarizedExperiment object, V can be a formula using the variables in the
#'   rowData slot of Y.
#' @param K integer. Number of latent factors(default 2).
#' @param which_assay numeric or character. Which assay of Y to use. 
#'   If missing, if `assayNames(Y)` contains "counts" then that is used.
#'   Otherwise, the first assay is used.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose Print helpful messages(default FALSE).
#' @param maxiter_optimize maximum number of iterations for the optimization
#'   step (default 100).
#' @param stop_epsilon stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#' @param children number of cores of the used cluster(default 1)
#' @param random_init if TRUE no initializations is done(default FALSE)
#' @param random_start if TRUE the setup of parameters is a random samplig
#' (default FALSE)
#' @param n_gene_disp number of genes used in mini-batch dispersion 
#' estimation approach(default NULL > all genes are used)
#' @param n_cell_par number of cells used in mini-batch cells related 
#' parameters estimation approach(default NULL > all cells are used)
#' @param n_gene_par number of genes used in mini-batch genes related 
#' parameters estimation approach(default NULL > all genes are used)
#'
#' @details For visualization (heatmaps, ...), please use the normalized values.
#' It corresponds to the deviance residuals when the \code{W} is not included
#' in the model but the gene and cell-level covariates are. As a results, when
#' \code{W} is not included in the model, the deviance residuals should capture
#' the biology. Note that we do not recommend to use the normalized values for
#' any downstream analysis (such as clustering, or differential expression), 
#' but only for visualization.
#'
#' @details If one has already fitted a model using \code{\link{newmodel}},
#' the object containing such model can be used as input of \code{newWave} to
#' save the resulting W into a \code{SummarizedExperiment} and optionally
#' compute residuals and normalized values, without the need for re-fitting the
#' model.
#'
#' @details By default \code{newWave} uses all genes to estimate \code{W}.
#'   However, we recommend to use the top 1,000 most variable genes for this
#'   step. In general, a user can specify any custom set of genes to be used to
#'   estimate \code{W}, by specifying either a vector of gene names, or a
#'   single character string corresponding to a column of the \code{rowData}.
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
#' m <- newWave(se, X="~bio")
setMethod("newWave", "SummarizedExperiment",
       function(Y, X, V, K=2, which_assay,
                    commondispersion = TRUE, verbose=FALSE,
                    maxiter_optimize=100,
                    stop_epsilon=.0001,  children = 1,
                    random_init = FALSE, random_start = FALSE,
                    n_gene_disp = NULL,
                    n_cell_par = NULL, n_gene_par = NULL, ...) {
              
            
            
        fitted_model <- newFit(Y, X, V, K,
                            which_assay, commondispersion,
                            verbose,
                            maxiter_optimize, stop_epsilon,
                            children, random_init, 
                            random_start, n_gene_disp,
                            n_cell_par, n_gene_par, ...)
              
              
              
        out <- as(Y, "SingleCellExperiment")
              
            if (numberFactors(fitted_model) > 0){
                W <- newW(fitted_model)
                colnames(W) <- paste0('W', seq_len(numberFactors(fitted_model)))
                reducedDim(out, "newWave") <- W
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
                    stop("assay needs to be a numeric or
                         character specifying which assay to use")
                } else {
                    dataY <- assay(Y, which_assay)
                }
            }
            
                return(out)
           }
)
