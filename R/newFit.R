#' @describeIn newFit Y is a
#'   \code{SummarizedExperiment}.
#' @export
#'
#' @import SharedObject
#' @import parallel
#' @importFrom stats model.matrix as.formula optimize rnbinom rnorm runif var
#' @param Y The SummarizedExperiment with the data
#' @param X The design matrix containing sample-level covariates, one sample per
#'   row. If missing, X will contain only an intercept. If Y is a
#'   SummarizedExperiment object, X can be a formula using the variables in the
#'   colData slot of Y.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept. If Y is a
#'   SummarizedExperiment object, V can be a formula using the variables in the
#'   rowData slot of Y.
#' @param K integer. Number of latent factors(default 2).
#' @param which_assay numeric or character. Which assay of Y to use. If missing,
#'   if `assayNames(Y)` contains "counts" then that is used. Otherwise, the
#'   first assay is used.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose Print helpful messages(default FALSE).
#' @param maxiter_optimize maximum number of iterations for the optimization
#'   step (default 100).
#' @param stop_epsilon stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#' @param children number of cores of the used cluster(default 1)
#' @param random_init if TRUE no initializations is done(default FALSE)
#' @param random_start if TRUE the setup of parameters is
#'  a random samplig(default FALSE)
#' @param n_gene_disp number of genes used in mini-batch dispersion estimation
#'  approach(default NULL > all genes are used)
#' @param n_cell_par number of cells used in mini-batch
#'  cell's related parameters estimation approach
#'  (default NULL > all cells are used)
#' @param n_gene_par number of genes used in mini-batch
#'  gene's related parameters estimation approach
#'  (default NULL > all genes are used)
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- newFit(se, X=model.matrix(~bio, data=colData(se)))
setMethod("newFit", "SummarizedExperiment",
          function(Y, X, V, K = 2, which_assay,
                    commondispersion = TRUE, verbose=FALSE,
                    maxiter_optimize = 100,
                    stop_epsilon=.0001, children = 1,
                    random_init = FALSE, random_start = FALSE,
                    n_gene_disp = NULL,
                    n_cell_par = NULL, n_gene_par = NULL, ... ) {

            if(missing(which_assay)) {
                if("counts" %in% assayNames(Y)) {
                    dataY <- assay(Y, "counts")
                } else {
                    dataY <- assay(Y)
                }
            } else {
                if(!(is.character(which_assay) | is.numeric(which_assay))) {
                    stop("Assay needs to be a numeric or character specifying
                         which assay to use")
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
                        stop("X must be a matrix or a formula with variables in
                             colData(Y)")
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
                        stop("V must be a matrix or a formula with variables in
                             rowData(Y)")
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
                        n_gene_par = n_gene_par)
              
            return(res)
})


#' @describeIn newFit Y is a matrix of counts (genes in rows).
#' @export
#'
#' @param Y The matrix with the data
#' @param X The design matrix containing sample-level covariates, one sample per
#'   row. If missing, X will contain only an intercept.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept.
#' @param K integer. Number of latent factors(default 2).
#' @param which_assay numeric or character. Which assay of Y to use. If missing,
#'   if `assayNames(Y)` contains "counts" then that is used. Otherwise, the
#'   first assay is used.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose Print helpful messages(default FALSE).
#' @param maxiter_optimize maximum number of iterations for the optimization
#'   step (default 100).
#' @param stop_epsilon stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#' @param children number of cores of the used cluster(default 1)
#' @param random_init if TRUE no initializations is done(default FALSE)
#' @param random_start if TRUE the setup of parameters is
#'  a random samplig(default FALSE)
#' @param n_gene_disp number of genes used in mini-batch dispersion estimation
#'  approach(default NULL > all genes are used)
#' @param n_cell_par number of cells used in mini-batch cell's related
#'  parameters estimation approach(default NULL > all cells are used)
#' @param n_gene_par number of genes used in mini-batch gene's related
#'  parameters estimation approach(default NULL > all genes are used)
#'
#' @details By default, i.e., if no arguments other than \code{Y} are passed,
#'   the model is fitted with an intercept for the regression across-samples and
#'   one intercept for the regression across genes.
#'
#' @details If Y is a Summarized experiment, the function uses the assay named
#'   "counts", if any, or the first assay.
#'
#' @seealso \code{\link[stats]{model.matrix}}.
#' @import irlba
#' @examples
#' bio <- gl(2, 3)
#' m <- newFit(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'              X=model.matrix(~bio))


setMethod("newFit", "matrix",
        function(Y, X, V, K = 2,
                commondispersion = TRUE, verbose=FALSE,
                maxiter_optimize = 100,
                stop_epsilon=.0001, children = 1,
                random_init = FALSE, random_start = FALSE,
                n_gene_disp = NULL,
                n_cell_par = NULL, n_gene_par = NULL, ... ) {
      
    # Check that Y contains whole numbers
    if(!all(.is_wholenumber(Y))) {
        stop("The input matrix should contain only whole numbers.")
    }

    # Transpose Y: UI wants genes in rows, internals genes in columns!
    Y_sh <- share(t(Y))
    
    if(any(rowSums(Y_sh) == 0)) {
        stop("Sample ", which(rowSums(Y_sh) == 0)[1], " has only 0 counts!")
    }
    
    if(any(colSums(Y_sh) == 0)) {
        stop("Gene ", which(colSums(Y_sh) == 0)[1], " has only 0 counts!")
    }
    
    # Create a newmodel object
    m <- newmodel(n=NROW(Y_sh), J=NCOL(Y_sh), K=K, X=X, V=V)

    cl <- makePSOCKcluster(children)
    on.exit(stopCluster(cl), add = TRUE)
    # Exporting values to the main and the child process
    
    # If the set the value of parameters is zero we must do the initialization
    if(random_init){ random_start = TRUE}
    
    m <- setup(cluster = cl, model = m, random_start = random_start,
        children = children, random_init = random_init, verbose = verbose,
        Y_sh = Y_sh)
    
    # Initializize value

    if (!random_init){

        initialization(cluster = cl, children = children, model = m,
                    verbose = verbose, Y = Y_sh)

    }
    
    
    # Optimize value

    info <- optimization(Y = Y_sh, cluster = cl, children = children, model = m,
                        max_iter = maxiter_optimize,
                        stop_epsilon = stop_epsilon,
                        commondispersion = commondispersion,
                        n_gene_disp = n_gene_disp, n_cell_par = n_cell_par,
                        n_gene_par = n_gene_par, verbose =  verbose)
    
    return(info)
})

#' @describeIn newFit Y is a DeleyedMatrix of counts (genes in rows).
#' @export
#'
#'
#' @import DelayedArray
#' @import BiocSingular

setMethod("newFit", "DelayedMatrix",
          function(Y, X, V, K = 2,
                   commondispersion = TRUE, verbose=FALSE,
                   maxiter_optimize = 100,
                   stop_epsilon=.0001, children = 1,
                   random_init = FALSE, random_start = FALSE,
                   n_gene_disp = NULL,
                   n_cell_par = NULL, n_gene_par = NULL, ... ) {
          
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
    # random_init = TRUE
    # random_start = TRUE

    # Exporting values to the main and the child process
    m <- setup(cluster = cl, model = m, random_start = random_start,
        children = children, random_init = random_init, verbose = verbose,
        Y_sh = Y)
    
    delayed_initialization(cluster = cl, children = children, model = m,
                   verbose = verbose, Y = Y)

    # Optimize value

    info <- delayed_optimization(Y=Y, cluster = cl, children = children,
        model = m, max_iter = maxiter_optimize,
        stop_epsilon = stop_epsilon,
        commondispersion = commondispersion,  n_gene_disp = n_gene_disp,
        n_cell_par = n_cell_par, n_gene_par = n_gene_par, verbose =  verbose)

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
#
#' There are different type of starting values:
#' 1.You can set all values to 0.
#' 2.You can sample values from a gaussian distribution or a
#' chisq distibution for the dispersion parameters,
#'
#' It creates different shared object and epxort them to the father
#' and the child process.
#' @param cluster the PSOCK cluster object
#' @param model The model of class newmodel
#' @param random_start if TRUE the setup of parameters is a 
#'   random samplig(default FALSE)
#' @param children Number of child process.
#' @param random_init if TRUE no initializations is done(default FALSE)
#' @param verbose Print helpful messages(default FALSE).
#' @param Y matrix of counts
#' @return A object of class newModel
#' @keywords internal


setup <- function(cluster, model, random_start, children,
                  random_init, verbose, Y_sh) {

    ptm <- proc.time()
    
    X_sh <- SharedObject::share(newX(model))
    V_sh <- SharedObject::share(newV(model))

    if (!random_start){
        beta_sh <- SharedObject::share(newBeta(model), copyOnWrite=FALSE)
        alpha_sh <- SharedObject::share(newAlpha(model), copyOnWrite=FALSE)
        W_sh <- SharedObject::share(newW(model), copyOnWrite=FALSE)
        gamma_sh <- SharedObject::share(newGamma(model), copyOnWrite=FALSE)
        zeta_sh <- SharedObject::share(newZeta(model), copyOnWrite=FALSE)
    } else {
        beta_sh <- SharedObject::share(matrix(rnorm(ncol(X_sh)*ncol(Y_sh)), 
            nrow = ncol(X_sh)), copyOnWrite=FALSE)
        alpha_sh <- SharedObject::share(matrix(rnorm(
            numberFactors(model)*ncol(Y_sh)), nrow = numberFactors(model)),
            copyOnWrite=FALSE)
        W_sh <- SharedObject::share(matrix(rnorm(nrow(Y_sh)*numberFactors(model)), 
            nrow = nrow(Y_sh)), copyOnWrite=FALSE)
        gamma_sh <- SharedObject::share(matrix(rnorm(nrow(Y_sh)*ncol(V_sh)), 
            nrow = ncol(V_sh)), copyOnWrite=FALSE)
        zeta_sh <- SharedObject::share(rep(rnorm(1), 
            length = ncol(Y_sh)), copyOnWrite=FALSE)
    }
    
    epsilon_gamma <- newEpsilon_gamma(model)
    epsilon_beta <- newEpsilon_beta(model)
    epsilon_alpha <- newEpsilon_alpha(model)
    epsilon_W <- newEpsilon_W(model)
    epsilonright <- c(newEpsilon_beta(model), newEpsilon_alpha(model))
    epsilonleft <- c(newEpsilon_gamma(model), newEpsilon_W(model))
    
    fun_path <- system.file("function.R", package = "NewWave")
    clusterExport(cluster, c("beta_sh" ,"alpha_sh","Y_sh","X_sh","W_sh","V_sh",
        "gamma_sh","zeta_sh", "epsilonright",
        "epsilonleft","children", "epsilon_gamma",
        "epsilon_beta", "epsilon_alpha", "epsilon_W", "fun_path"),
        envir = environment())
    
    clusterEvalQ(cluster, source(fun_path))
    
    m <- newmodel(X = X_sh, V = V_sh,
                  W = W_sh, beta = beta_sh,
                  gamma = gamma_sh,
                  alpha = alpha_sh, zeta = zeta_sh,
                  epsilon_beta = model@epsilon_beta,
                  epsilon_gamma = model@epsilon_gamma,
                  epsilon_W = model@epsilon_W,
                  epsilon_alpha = model@epsilon_alpha,
                  epsilon_zeta = model@epsilon_zeta)

    if(verbose){
        cat("Time of setup\n")
        init.t <- proc.time()
        print(init.t-ptm)
    }
    
    return(m)
}

#' Initialize the parameters of a Negative Binomial regression model
#'
#' It initialize gamma and beta using a Ridge Regression and W and alpha using PCA
#' @param cluster The PSOCK cluster
#' @param children Number of child process.
#' @param model The newmodel object
#' @param verbose Print proc time
#' @return It does not return anything, the parameters
#'  are update internally as these are shared object
#' @keywords internal
initialization <- function(cluster, children, model, verbose, Y){
    
  
    ptm <- proc.time()
    
    L_sh <- log1p(Y)
    
    clusterExport(cluster,"L_sh",envir = environment())
    
    clusterApply(cluster, seq.int(children), "gamma_init")

    clusterApply(cluster, seq.int(children), "beta_init")

  
    D <- L_sh - (newX(model) %*% newBeta(model)) - t(newV(model) %*% newGamma(model))


    
    R <- irlba::irlba(D, nu=numberFactors(model), nv=numberFactors(model))
    
    

    # Orthogonalize to get W and alpha
    model@W[] <- (newEpsilon_alpha(model) / newEpsilon_W(model))[1]^(1/4) *
        R$u %*% diag(sqrt(R$d[seq.int(numberFactors(model))]),
        nrow = length(R$d[seq.int(numberFactors(model))]))
    
    model@alpha[] <- (newEpsilon_W(model)/newEpsilon_alpha(model))[1]^(1/4) *
        diag(sqrt(R$d[seq.int(numberFactors(model))]),
        nrow = length(R$d[seq.int(numberFactors(model))])) %*% t(R$v)
    
    
    if(verbose){
        cat("Time of initialization\n")
        init.t <- proc.time()
        print(init.t-ptm)
    }
}


#' Initialize the parameters of a Negative Binomial regression model with a 
#' DelayedArray object
#'
#' It initialize gamma and beta using a Ridge Regression and W and alpha using PCA
#' @param cluster The PSOCK cluster
#' @param children Number of child process.
#' @param model The newmodel object
#' @param verbose Print proc time
#' @param Y Is the data matrix
#' @return It does not return anything, the parameters
#'  are update internally as these are shared object
#' @keywords internal

delayed_initialization <- function(cluster, children, model, verbose, Y){
  
  
  ptm <- proc.time()
  
  L_sh <- log1p(Y)
  
  clusterExport(cluster,"L_sh",envir = environment())
  
  clusterApply(cluster, seq.int(children), "delayed_gamma_init")
  
  clusterApply(cluster, seq.int(children), "delayed_beta_init")

  cl <- as(cluster,"SnowParam")
  R <- BiocSingular::runExactSVD(L_sh, k=numberFactors(model), BPPARAM = cl)
  
  
  # Orthogonalize to get W and alpha
  model@W[] <- (newEpsilon_alpha(model) / newEpsilon_W(model))[1]^(1/4) *
    R$u %*% diag(sqrt(R$d[seq.int(numberFactors(model))]),
                 nrow = length(R$d[seq.int(numberFactors(model))]))
  
  model@alpha[] <- (newEpsilon_W(model)/newEpsilon_alpha(model))[1]^(1/4) *
    diag(sqrt(R$d[seq.int(numberFactors(model))]),
         nrow = length(R$d[seq.int(numberFactors(model))])) %*% t(R$v)
  
  
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
#' @param cluster The PSOCK cluster
#' @param children Number of child process
#' @param model newmodel item
#' @param max_iter maximum number of iterations
#' @param stop_epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon
#' @param n_gene_disp number of genes used in mini-batch dispersion estimation
#'   approach(default NULL > all genes are used)
#' @param n_cell_par number of cells used in mini-batch cell's related
#'   parameters estimation approach(default NULL > all cells are used)
#' @param n_gene_par number of genes used in mini-batch gene's related
#'   parameters estimation approach(default NULL > all genes are used)
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose print information (default FALSE)
#' @return An object of class newmodel similar to the one given as argument
#'   with modified parameters alpha, beta, gamma, W.
#' @keywords internal
#' 
optimization <- function(Y, cluster, children, model ,
                        max_iter, stop_epsilon,
                        n_gene_disp,
                        n_cell_par, n_gene_par,
                        commondispersion, verbose){

    iter = 0

    total.lik=rep(NA,max_iter)
    
    Y_sh <- Y
    alpha_sh <- newAlpha(model)
    beta_sh <-newBeta(model)
    gamma_sh <- newGamma(model)
    W_sh <- newW(model)
    X_sh <- newX(model)
    V_sh <- newV(model)
    zeta_sh <-newZeta(model)
    
    mu_sh <- SharedObject::share(exp(X_sh %*% beta_sh +
                t(V_sh %*% gamma_sh) +  W_sh %*% alpha_sh))
    clusterExport(cl = cluster, "mu_sh",
                envir = environment())
    total.lik[1] <- ll_calc(mu = mu_sh, model  = model, Y_sh = Y_sh,
                z = zeta_sh,
                alpha_sh, beta_sh, gamma_sh, W_sh, commondispersion)

    for (iter in seq.int(max_iter)){

        ptm <- proc.time()


        if(iter > 1){
  
            mu_sh[] <- exp(newX(model) %*% beta_sh + t(newV(model) %*% gamma_sh) + 
                W_sh %*% alpha_sh)
  
            total.lik[iter] <- ll_calc(mu = mu_sh, model  = model,
                Y_sh = Y_sh, z = zeta_sh, alpha_sh, beta_sh,
                gamma_sh, W_sh, commondispersion)
  
            if(abs((total.lik[iter]-total.lik[iter-1]) / total.lik[iter-1])<stop_epsilon){
              m <- newmodel(X = unshare(newX(model)), V = unshare(newV(model)),
                            W = unshare(W_sh), beta = unshare(beta_sh),
                            gamma = unshare(gamma_sh),
                            alpha = unshare(alpha_sh), zeta = unshare(zeta_sh),
                            epsilon_beta = model@epsilon_beta,
                            epsilon_gamma = model@epsilon_gamma,
                            epsilon_W = model@epsilon_W,
                            epsilon_alpha = model@epsilon_alpha,
                            epsilon_zeta = model@epsilon_zeta)
              break
            }
  
        }

        if(verbose){
          
            message("Iteration ",iter)
            message("penalized log-likelihood = ",  total.lik[iter])
            
         }

        optimd(Y_sh, mu = mu_sh, cluster = cluster,
               children = children, commondispersion = commondispersion,
               num_gene = n_gene_disp, iter = iter, alpha_sh = alpha_sh,
               beta_sh = beta_sh, gamma_sh = gamma_sh, 
               W_sh = W_sh, zeta_sh = zeta_sh, model, total.lik)

        if(verbose){
            cat("Time of dispersion optimization\n")
            print(proc.time()-ptm)

            l_pen <- ll_calc(mu = mu_sh, model  = model, Y_sh = Y_sh,
                z = zeta_sh, alpha_sh, beta_sh, gamma_sh, W_sh,
                commondispersion)
            
            message("after optimize dispersion = ", l_pen )
          
          }

        ptm <- proc.time()

        
        clusterApply(cluster, seq.int(children), "optimr",  
                num_gene = n_gene_par, iter = iter)


        if(verbose){
            cat("Time of right optimization\n")
            print(proc.time()-ptm)
            itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                W_sh %*% alpha_sh)
            l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh,
                z = zeta_sh, alpha_sh, beta_sh, gamma_sh, W_sh,
                commondispersion)
            message("after right optimization= ",  l_pen)
        }

        o <- orthogonalizeTraceNorm(W_sh, alpha_sh, newEpsilon_W(model),
                                    newEpsilon_alpha(model))
        
        W_sh[] <- o$U
        alpha_sh[] <- o$V


        if(verbose){
            itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                W_sh %*% alpha_sh)
            l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh,
                z = zeta_sh, alpha_sh, beta_sh, gamma_sh, W_sh,
                commondispersion)
            message("after orthogonalization = ",  l_pen)
        }


        ptm <- proc.time()

        
        clusterApply(cluster, seq.int(children), "optiml" , 
            num_cell = n_cell_par, iter = iter) 


        if(verbose){
            cat("Time of left optimization\n")
            print(proc.time()-ptm)
            itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                W_sh %*% alpha_sh)
            l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh,
                z = zeta_sh, alpha_sh, beta_sh, gamma_sh, W_sh,
                commondispersion)
            message("after left optimization= ",  l_pen)
        }

        o <- orthogonalizeTraceNorm(W_sh, alpha_sh, newEpsilon_W(model),
                                    newEpsilon_alpha(model))
        W_sh[] <- o$U
        alpha_sh[] <- o$V


        if(verbose){
            itermu <- exp(X_sh %*% beta_sh + t(V_sh %*% gamma_sh) +
                W_sh %*% alpha_sh)
            l_pen <- ll_calc(mu = itermu, model  = model, Y_sh = Y_sh,
                z = zeta_sh, alpha_sh, beta_sh, gamma_sh, W_sh,
                commondispersion)
            message("after orthogonalization = ",  l_pen)
        }

        m <- newmodel(X = newX(model), V = newV(model),
            W = W_sh, beta = beta_sh,
            gamma = gamma_sh,
            alpha = alpha_sh, zeta = zeta_sh,
            epsilon_beta = model@epsilon_beta,
            epsilon_gamma = model@epsilon_gamma,
            epsilon_W = model@epsilon_W,
            epsilon_alpha = model@epsilon_alpha,
            epsilon_zeta = model@epsilon_zeta)
    }

    return(m)
}

#' Optimize the parameters of a Negative Binomial regression model with Delayed Array object
#'
#' The parameters of the model given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument.
#' @param cluster The PSOCK cluster
#' @param children Number of child process
#' @param model newmodel item
#' @param max_iter maximum number of iterations
#' @param stop_epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon
#' @param n_gene_disp number of genes used in mini-batch dispersion estimation
#'   approach(default NULL > all genes are used)
#' @param n_cell_par number of cells used in mini-batch cell's related
#'   parameters estimation approach(default NULL > all cells are used)
#' @param n_gene_par number of genes used in mini-batch gene's related
#'   parameters estimation approach(default NULL > all genes are used)
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param verbose print information (default FALSE)
#' @return An object of class newmodel similar to the one given as argument
#'   with modified parameters alpha, beta, gamma, W.
#' @keywords internal
#' 
delayed_optimization <- function(Y, cluster, children, model ,
                         max_iter, stop_epsilon,
                         n_gene_disp,
                         n_cell_par, n_gene_par,
                         commondispersion, verbose){
  
  iter = 0
  
  total.lik=rep(NA,max_iter)
  
  Y_sh <- Y
  alpha_sh <- newAlpha(model)
  beta_sh <-newBeta(model)
  gamma_sh <- newGamma(model)
  W_sh <- newW(model)
  X_sh <- newX(model)
  V_sh <- newV(model)
  zeta_sh <-newZeta(model)
  
  
  total.lik[1] <- delayed_calc(cluster, children, model)
  
  for (iter in seq.int(max_iter)){
    
    ptm <- proc.time()
    
    
    if(iter > 1){
      
      total.lik[iter] <- delayed_calc(cluster, children, model)
      
      if(abs((total.lik[iter]-total.lik[iter-1]) / total.lik[iter-1])<stop_epsilon) break
      
      
    }
    
    if(verbose){
      
      message("Iteration ",iter)
      message("penalized log-likelihood = ",  total.lik[iter])
    }
    
    clusterApply(cluster, seq.int(children), "optim_genwise_dispersion_delayed",
                 Y,num_gene = n_gene_disp, iter = iter)
    
    if(verbose){
      cat("Time of dispersion optimization\n")
      print(proc.time()-ptm)
      message("after optimize dispersion = ",  delayed_calc(cluster, children, model))
    }
    
    ptm <- proc.time()
    
    clusterApply(cluster, seq.int(children), "optimr_delayed",  
                      num_gene = n_gene_par, iter = iter)
    
    
    if(verbose){
      cat("Time of right optimization\n")
      print(proc.time()-ptm)
      message("after right optimization= ",  delayed_calc(cluster, children, model))
    }
    
    o <- orthogonalizeTraceNorm(W_sh, alpha_sh, newEpsilon_W(model),
                                newEpsilon_alpha(model))
    
    W_sh[] <- o$U
    alpha_sh[] <- o$V
    
    
    if(verbose){
      message("after orthogonalization = ",  delayed_calc(cluster, children, model))
    }
    
    
    ptm <- proc.time()
    
    clusterApply(cluster, seq.int(children), "optiml_delayed" , 
                      num_cell = n_cell_par, iter = iter)
    
    
    if(verbose){
      cat("Time of left optimization\n")
      print(proc.time()-ptm)
      message("after left optimization= ",  delayed_calc(cluster, children, model))
      
    }
    
    o <- orthogonalizeTraceNorm(W_sh, alpha_sh, newEpsilon_W(model),
                                newEpsilon_alpha(model))
    W_sh[] <- o$U
    alpha_sh[] <- o$V
    
    
    if(verbose){
    
      message("after orthogonalization = ",  delayed_calc(cluster, children, model))
    }
    
    m <- newmodel(X = newX(model), V = newV(model),
                  W = W_sh, beta = beta_sh,
                  gamma = gamma_sh,
                  alpha = alpha_sh, zeta = zeta_sh,
                  epsilon_beta = model@epsilon_beta,
                  epsilon_gamma = model@epsilon_gamma,
                  epsilon_W = model@epsilon_W,
                  epsilon_alpha = model@epsilon_alpha,
                  epsilon_zeta = model@epsilon_zeta)
  }
  
  return(m)
}

optimd <- function(Y_sh, mu, cluster, children, num_gene = NULL, commondispersion,
                   iter, alpha_sh,beta_sh,gamma_sh, 
                   W_sh, zeta_sh, model, total.lik){
    
    J <- ncol(Y_sh)
    
    if (commondispersion || iter == 1){
  
        genes = seq.int(J)
  
    if(!is.null(num_gene) && iter != 1){
    
        genes <- sample(x = J, size = num_gene)
        
        g=optimize(f=nb.loglik.dispersion, Y=Y_sh[,genes], mu=mu[,genes],
                   maximum=TRUE,
                   interval=c(-50,50))
        
        l_pen <- ll_calc(mu = mu, model  = model, Y_sh = Y_sh,
                         z = rep(g$maximum,J), alpha_sh, beta_sh, gamma_sh, W_sh,
                         commondispersion)
        
        if(l_pen > total.lik[iter]){
          zeta_sh[] <- rep(g$maximum,J)
        }
    
    } else{
    
    
        g=optimize(f=nb.loglik.dispersion, Y=Y_sh[,genes], mu=mu[,genes],
                   maximum=TRUE,
                   interval=c(-50,50))
        
        
          zeta_sh[] <- rep(g$maximum,J)
        
    }
  
    } else {
  
        clusterApply(cluster, seq.int(children), "optim_genwise_dispersion",
                     num_gene = num_gene, iter = iter)
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



ll_calc <- function(mu, model, Y_sh, z, alpha , beta, gamma, W,
                    commondispersion){
  
    theta <- exp(z)
  
    loglik <- nb.loglik(Y_sh, mu, rep(theta, rep(nrow(Y_sh),ncol(Y_sh))))
  
    zeta_pen <- ifelse(commondispersion == TRUE,
                       newEpsilon_zeta(model)*var(z)/2, 0)
  
    penalty <- sum(newEpsilon_alpha(model) * (alpha)^2)/2 +
       sum(newEpsilon_beta(model) * (beta)^2)/2 +
       sum(newEpsilon_gamma(model)*(gamma)^2)/2 +
       sum(newEpsilon_W(model)*t(W)^2)/2 +
       zeta_pen
  
    loglik - penalty
}

delayed_calc <- function(cluster,children, m){
  
 ll <-sum(unlist(clusterApply(cluster, seq.int(children), "delayed_ll")))
 penalty <- sum(newEpsilon_alpha(m) * (newAlpha(m))^2)/2 +
   sum(newEpsilon_beta(m) * (newBeta(m))^2)/2 +
   sum(newEpsilon_gamma(m)*(newGamma(m))^2)/2 +
   sum(newEpsilon_W(m)*t(newW(m))^2)/2
 ll - penalty
}


.is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}
