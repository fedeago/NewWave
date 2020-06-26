#' @describeIn newFit Y is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @export
#'
#' @import SharedObject
#' @import SummarizedExperiment
#' @importFrom stats model.matrix as.formula
#' 
#' @param which_assay numeric or character. Which assay of Y to use (only if Y
#'   is a SummarizedExperiment).
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- newFit(se, X=model.matrix(~bio, data=colData(se)))
setMethod("newFit", "SummarizedExperiment",
          function(Y, X, V, K = 2, which_assay,
                   commondispersion = T, verbose=FALSE,
                   maxiter_optimize = 100,
                   stop_epsilon=.0001, children = 1,
                   random_init = FALSE, random_start = FALSE,
                   n_gene_disp = NULL,
                   n_cell_par = NULL, n_gene_par = NULL,
                   cross_batch = F, ... ) {

              if(missing(which_assay)) {
                  if("counts" %in% assayNames(Y)) {
                      dataY <- assay(Y, "counts")
                  } else {
                      dataY <- assay(Y)
                  }
              } else {
                  if(!(is.character(which_assay) | is.numeric(which_assay))) {
                      stop("Assay needs to be a numeric or character specifying which assay to use")
                  } else {
                      dataY <- assay(Y, which_assay)
                  }
              }

              if(!missing(X)) {
                  if(!is.matrix(X)) {
                      tryCatch({
                          f <- as.formula(X)
                          X <- model.matrix(f, data=colData(Y))
                      },
                      error = function(e) {
                          stop("X must be a matrix or a formula with variables in colData(Y)")
                      })
                  }
              }

              if(!missing(V)) {
                  if(!is.matrix(V)) {
                      tryCatch({
                          f <- as.formula(V)
                          V <- model.matrix(f, data=rowData(Y))
                      },
                      error = function(e) {
                          stop("V must be a matrix or a formula with variables in rowData(Y)")
                      })
                  }
              }

              # Apply newFit on the assay of SummarizedExperiment
              res <- newFit(Y = dataY, X = X, V = V, K = K,
                           commondispersion = commondispersion, 
                           verbose = verbose,
                           maxiter_optimize = maxiter_optimize,
                           stop_epsilon = stop_epsilon, children = children,
                           random_init = random_init, random_start = random_start,
                           n_gene_disp = n_gene_disp , n_cell_par = n_cell_par,
                           n_gene_par = n_gene_par,cross_batch=cross_batch)
              
              
               
              return(res)
})


#' @describeIn newFit Y is a matrix of counts (genes in rows).
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
#' @param K integer. Number of latent factors..
#' @param verbose Print helpful messages.
#' @param nb.repeat.initialize Number of iterations for the initialization of
#'   beta and gamma.
#' @param maxiter.optimize maximum number of iterations for the optimization
#'   step (default 25).
#' @param stop.epsilon.optimize stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#'
#' @details By default, i.e., if no arguments other than \code{Y} are passed,
#'   the model is fitted with an intercept for the regression across-samples and
#'   one intercept for the regression across genes.
#'
#' @details If Y is a Summarized experiment, the function uses the assay named
#'   "counts", if any, or the first assay.
#'
#' @seealso \code{\link[stats]{model.matrix}}.
#'
#' @import BiocParallel
#' @examples
#' bio <- gl(2, 3)
#' m <- newFit(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'              X=model.matrix(~bio))


setMethod("newFit", "matrix",
          function(Y, X, V, K = 2,
                   commondispersion = T, verbose=FALSE,
                   maxiter_optimize = 100,
                   stop_epsilon=.0001, children = 1,
                   random_init = FALSE, random_start = FALSE,
                   n_gene_disp = NULL,
                   n_cell_par = NULL, n_gene_par = NULL,
                   cross_batch = F, ... ) {

    # Check that Y contains whole numbers
    if(!all(.is_wholenumber(Y))) {
        stop("The input matrix should contain only whole numbers.")
    }

    # Transpose Y: UI wants genes in rows, internals genes in columns!
    Y <- t(Y)
    
    if(any(rowSums(Y) == 0)) {
      stop("Sample ", which(rowSums(Y) == 0)[1], " has only 0 counts!")
    }
    
    if(any(colSums(Y) == 0)) {
      stop("Gene ", which(colSums(Y) == 0)[1], " has only 0 counts!")
    }
    
    # Create a newmodel object
    m <- newmodel(n=NROW(Y), J=NCOL(Y), K=K, X=X, V=V)

    cl <- makePSOCKcluster(children)
    on.exit(stopCluster(cl), add = TRUE)
    on.exit(CleanEnvir(), add=TRUE)
    # Exporting values to the main and the child process
    
    # If the set the value of parameters is zero we must do the initialization
    if(random_init){ random_start = T}
    
    setup(cluster = cl, model = m, random_start = random_start, children = children,
          random_init = random_init, verbose = verbose, Y = Y,mode = "matrix")
    
    # Initializize value

    if (!random_init){

      initialization(cluster = cl, children = children, model = m,
                     verbose = verbose)

    }
    
    orthog <- (nFactors(m)>0)

    # Optimize value

    info <- optimization(cluster = cl, children = children, model = m, max_iter = maxiter_optimize,
                         orthog = orthog, stop_epsilon = stop_epsilon,
                         commondispersion = commondispersion,  n_gene_disp = n_gene_disp,
                         n_cell_par = n_cell_par, n_gene_par = n_gene_par, verbose =  verbose,
                         mode = "matrix", cross_batch = cross_batch)
    
    return(info)
})



setMethod("newFit", "DelayedMatrix",
          function(Y, X, V, K,
                   commondispersion, verbose,
                   nb_repeat, maxiter_optimize,
                   stop_epsilon, children,
                   random_start, n_gene_disp,
                   n_cell_par, n_gene_par,
                   ... ) {
            
            # if(!all(.is_wholenumber(Y))) {
            #   stop("The input matrix should contain only whole numbers.")
            # }
            
            # Transpose Y: UI wants genes in rows, internals genes in columns!
            Y <- t(Y)
            
            # Create a newmodel object
            m <- newmodel(n=NROW(Y), J=NCOL(Y), K=K, X=X)
            
            cl <- makePSOCKcluster(children)
            on.exit(stopCluster(cl), add = TRUE)
            clusterEvalQ(cl, library(DelayedArray))
            
            # If the set the value of parameters is zero we must do the initialization
            random_init = T
            random_start = T
            
            # Exporting values to the main and the child process
            setup(cluster = cl, model = m, random_start = random_start, children = children,
                  random_init = random_init, verbose = verbose, Y = Y, mode = "Deleyed")
            
            
            orthog <- (nFactors(m)>0)
            
            # Optimize value
            
            info <- optimization(cluster = cl, children = children, model = m, max_iter = maxiter_optimize,
                                 orthog = orthog, stop_epsilon = stop_epsilon,
                                 commondispersion = commondispersion,  n_gene_disp = n_gene_disp,
                                 n_cell_par = n_cell_par, n_gene_par = n_gene_par, verbose =  verbose,
                                 mode = "Delayed")
            
            return(info)
          })
            
#' @describeIn newFit Y is a sparse matrix of counts (genes in rows).
#' @export
#'
#' @details Currently, if Y is a sparseMatrix, this calls the newFit method on
#'   as.matrix(Y)
#'
#' @importClassesFrom Matrix dgCMatrix
setMethod("newFit", "dgCMatrix",
          function(Y, ...) {
            newFit(as.matrix(Y), ...)
          })


#' Setup the parameters of a Negative Binomial regression model
#'
#' There are different type of starting values:
#' 1.You can set all values to 0.
#' 2.You can sample values from a gaussian distribution or a
#' chisq distibution for the dispersion parameters,
#'
#' It creates different shared object and exort them to the father
#' and the child process.
#' @param m The model of class newmodel
#' @param Y The matrix of counts.
#' @param random_start The type of start.
#' @param children Number of child process.
#' @examples
#' Y <- matrix(rpois(60, lambda=2), 6, 10)
#' bio <- gl(2, 3)
#' time <- rnorm(6)
#' gc <- rnorm(10)
#' m <- newmodel(Y, X=model.matrix(~bio + time), V=model.matrix(~gc), K=1)
#' m <- nbInitialize(m, Y)
#'

setup <- function(cluster, model, random_start, children,
                  random_init, verbose, Y, mode) {

  ptm <- proc.time()
  
  if (mode == "matrix"){
    Y_sh <<- share(Y)
    if(!random_init){
    L_sh <<- share(log1p(Y_sh))
    clusterExport(cluster,"L_sh",envir = environment())
    }
  } else {
    Y_sh <<- Y
    L_sh <<- log1p(Y_sh)
  }
  X_sh <<- share(model@X)
  V_sh <<- share(model@V)
  if (!random_start){
    beta_sh <<- share(model@beta, copyOnWrite=FALSE)
    alpha_sh <<- share(model@alpha, copyOnWrite=FALSE)
    W_sh <<- share(model@W, copyOnWrite=FALSE)
    gamma_sh <<- share(model@gamma, copyOnWrite=FALSE)
    zeta_sh <<- share(model@zeta, copyOnWrite=FALSE)
  } else {
    beta_sh <<- share(matrix(rnorm(ncol(X_sh)*ncol(Y_sh)), nrow = ncol(X_sh)), copyOnWrite=FALSE)
    alpha_sh <<- share(matrix(rnorm(nFactors(model)*ncol(Y_sh)), nrow = nFactors(model)), copyOnWrite=FALSE)
    W_sh <<- share(matrix(rnorm(nrow(Y_sh)*nFactors(model)), nrow = nrow(Y_sh)), copyOnWrite=FALSE)
    gamma_sh <<- share(matrix(rnorm(nrow(Y_sh)*ncol(V_sh)), nrow = ncol(V_sh)), copyOnWrite=FALSE)
    zeta_sh <<- share(rep(rnorm(1), length = ncol(Y_sh)), copyOnWrite=FALSE)
  }
  epsilon_gamma <- getEpsilon_gamma(model)
  epsilon_beta <- getEpsilon_beta(model)
  epsilonright <- c(getEpsilon_beta(model), getEpsilon_alpha(model))
  epsilonleft <- c(getEpsilon_gamma(model), getEpsilon_W(model))

  n = nSamples(model)
  J = nFeatures(model)
  
  
  #Cluster
  fun_path <- system.file("function.R", package = "NewWave")
  clusterExport(cluster, c("beta_sh" ,"alpha_sh","Y_sh","X_sh","W_sh","V_sh","gamma_sh","zeta_sh", "epsilonright",
                           "epsilonleft","children", "epsilon_gamma",
                           "epsilon_beta", "fun_path"),
                envir = environment())
  clusterEvalQ(cluster, source(fun_path))

  
  if(verbose){
  cat("Time of setup\n")
  init.t <- proc.time()
  print(init.t-ptm)
  }
}


initialization <- function(cluster, children, model, verbose){

  ptm <- proc.time()

    clusterApply(cluster, seq.int(children), "gamma_init")

    clusterApply(cluster, seq.int(children), "beta_init")

    
  D <- L_sh - (X_sh %*% beta_sh) - t(V_sh %*% gamma_sh)
  
  
  set.seed(1234)
  R <- irlba::irlba(D, nu=nFactors(model), nv=nFactors(model))


  # Orthogonalize to get W and alpha
  W_sh[] <- (getEpsilon_alpha(model) / getEpsilon_W(model))[1]^(1/4) *
    R$u %*% diag(sqrt(R$d[1:nFactors(model)]), nrow = length(R$d[1:nFactors(model)]))
  alpha_sh[] <- (getEpsilon_W(model)/getEpsilon_alpha(model))[1]^(1/4) *
    diag(sqrt(R$d[1:nFactors(model)]),nrow = length(R$d[1:nFactors(model)])) %*% t(R$v)

  if(verbose){
  cat("Time of initialization\n")
  init.t <- proc.time()
  print(init.t-ptm)
  }
}


#' Optimize the parameters of a Negative Binomial regression model
#'
#' The parameters of the model given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument.
#' @param m The model of class newmodel
#' @param Y The matrix of counts.
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @return An object of class newmodel similar to the one given as argument
#'   with modified parameters alpha, beta, gamma, W.
#' @examples
#' Y = matrix(10, 3, 5)
#' m = newmodel(n=NROW(Y), J=NCOL(Y))
#' m = nbInitialize(m, Y)
#' m = nbOptimize(m, Y)
#' @export

optimization <- function(cluster, children, model ,
                         max_iter, stop_epsilon,
                         n_gene_disp,
                         n_cell_par, n_gene_par, orthog,
                         commondispersion, verbose, mode, cross_batch){

  iter = 0
  
  total.lik=rep(NA,max_iter)
  
  
  mu_sh <<- share(exp(getX(model) %*% beta_sh + t(getV(model) %*% gamma_sh) +
                     W_sh %*% alpha_sh))
  clusterExport(cl = cluster, "mu_sh",
                envir = environment())
  total.lik[1] <- ll_calc(mu = mu_sh, model  = model, Y_sh = as.matrix(Y_sh), z = zeta_sh,
                          alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
  
  for (iter in seq.int(max_iter)){

    ptm <- proc.time()

    
    if(iter > 1){
      
      if(abs((total.lik[iter]-total.lik[iter-1]) / total.lik[iter-1])<stop_epsilon) break
      
      mu_sh[] <- exp(getX(model) %*% beta_sh + t(getV(model) %*% gamma_sh) + W_sh %*% alpha_sh)
    }
    
    if(verbose){
    message("Iteration ",iter)
    message("penalized log-likelihood = ",  total.lik[iter])
    }

    optimd(ncol(Y_sh), mu = mu_sh, cluster = cluster,
           children = children, commondispersion = commondispersion,
           num_gene = n_gene_disp, iter = iter)

    if(verbose){
    cat("Time of dispersion optimization\n")
    print(proc.time()-ptm)

    l_pen <- ll_calc(mu = mu_sh, model  = model, Y_sh = Y_sh, z = zeta_sh,
                     alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
    message("after optimize dispersion = ",  l_pen)
    }

    ptm <- proc.time()

    if(mode == "matrix"){
      clusterApply(cluster, seq.int(children), "optimr",  num_gene = n_gene_par, cross_batch = cross_batch, iter = iter)
    } else clusterApply(cluster, seq.int(children), "optimr_delayed",  num_gene = n_gene_par)


    if(verbose){
    cat("Time of right optimization\n")
    print(proc.time()-ptm)
    itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                    W_sh %*% alpha_sh)
    l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh, z = zeta_sh,
                     alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
    message("after right optimization= ",  l_pen)
    }

    if (orthog) {
      o <- orthogonalizeTraceNorm(W_sh, alpha_sh, model@epsilon_W, model@epsilon_alpha)
      W_sh[] <- o$U
      alpha_sh[] <- o$V
    }

    if(verbose){
    itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                      W_sh %*% alpha_sh)
    l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh, z = zeta_sh,
                     alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
    message("after orthogonalization = ",  l_pen)
    }

    
    ptm <- proc.time()

    if(mode == "matrix"){
    llikelihood <- clusterApply(cluster, seq.int(children), "optiml" , num_cell = n_cell_par, cross_batch = cross_batch, iter = iter)
    } else clusterApply(cluster, seq.int(children), "optiml_delayed" , num_cell = n_cell_par)
    
    total.lik[iter+1] <-sum(unlist(llikelihood))
    
    if(verbose){
    cat("Time of left optimization\n")
    print(proc.time()-ptm)
    itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                    W_sh %*% alpha_sh)
    l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh, z = zeta_sh,
                     alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
    message("after left optimization= ",  l_pen)
    }

    if (orthog) {
      o <- orthogonalizeTraceNorm(W_sh, alpha_sh, model@epsilon_W, model@epsilon_alpha)
      W_sh[] <- o$U
      alpha_sh[] <- o$V
    }

    if(verbose){
    itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                      W_sh %*% alpha_sh)
    l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh, z = zeta_sh,
                     alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)
    message("after orthogonalization = ",  l_pen)
    }

    m <- newmodel(X = model@X, V = model@V,
                   W = W_sh, beta = beta_sh,
                   gamma = gamma_sh,
                   alpha = alpha_sh, zeta = zeta_sh,
                   epsilon_beta = model@epsilon_beta,
                   epsilon_gamma = model@epsilon_gamma,
                   epsilon_W = model@epsilon_W, epsilon_alpha = model@epsilon_alpha,
                   epsilon_zeta = model@epsilon_zeta)
  }

  return(m)
}

optimd <- function(J, mu, cluster, children, num_gene = NULL, commondispersion, iter){
  
  if (commondispersion || iter == 1){
    
    genes = seq.int(J)
    
    if(!is.null(num_gene) && iter != 1){
      
      genes <- sample(x = J, size = num_gene)
      mu <- mu[,genes]
      
    }
    
    g=optimize(f=nb.loglik.dispersion, Y=as.matrix(Y_sh[,genes]), mu=mu,
               maximum=TRUE,
               interval=c(-50,50))
    
    zeta_sh[] <- rep(g$maximum,J)
    
  } else {
    
    clusterApply(cluster, seq.int(children), "optim_genwise_dispersion",num_gene = num_gene, iter = iter)
  }
  
}


orthogonalizeTraceNorm <- function(U, V, a=1, b=1) {

  # do QR of U
  U.qr <- qr (U)
  U.Q <- qr.Q (U.qr)
  U.R <- qr.R (U.qr)

  # do QR of t(V)
  V.qr <- qr (t(V))
  V.Q <- qr.Q (V.qr)
  V.R <- qr.R (V.qr)

  # do SVD of the U.R %*% t(V.R) matrix to have orthog %*% diag %*% orthog
  A <- svd( U.R %*% t(V.R) )

  # Scaling factor
  s <- (a/b)^{1/4}

  # orthogonalized U
  U2 <- 1/s * U.Q %*% A$u %*% sqrt(diag(A$d,length(A$d)))

  # orthogonalized V and W
  V2 <- s * t( V.Q %*% A$v %*% sqrt(diag(A$d,length(A$d))) )

  list(U=U2, V=V2)
}


nb.loglik <- function(Y, mu, theta) {
  
  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  sum(logPnb)
  
}

nb.loglik.dispersion <- function(zeta, Y, mu){

  nb.loglik(Y, mu, exp(zeta))

}


ll_calc <- function(mu, model, Y_sh, z, alpha , beta, gamma, W, commondispersion){
  
  theta <- exp(z)
  
  loglik <- nb.loglik(Y_sh, mu, rep(theta, rep(nrow(Y_sh),ncol(Y_sh))))
  
  zeta_pen <- ifelse(commondispersion == T,  getEpsilon_zeta(model)*var(z)/2, 0)
  
  penalty <- sum(getEpsilon_alpha(model) * (alpha)^2)/2 +
    sum(getEpsilon_beta(model) * (beta)^2)/2 +
    sum(getEpsilon_gamma(model)*(gamma)^2)/2 +
    sum(getEpsilon_W(model)*t(W)^2)/2 +
    zeta_pen
  
  loglik - penalty
}


CleanEnvir <- function() {
  objs <- ls(pos = ".GlobalEnv")
  rm(list = objs[grep("_sh$", objs)], pos = ".GlobalEnv")
}


.is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}