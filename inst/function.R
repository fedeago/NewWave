solveRidgeRegression <- function(mat, y, beta=rep(0,NCOL(x)), epsilon=1e-6) {

    # loglik
    f <-function(b) {
        eta <- mat %*% b
        l <- sum((y-eta)^2)/2
        l + sum(epsilon*b^2)/2
    }
  
    # gradient of loglik
    g <-function(b) {
        eta <- mat %*% b
        l <- t(mat) %*% (eta - y)
        l + epsilon*b
    }
  
    # optimize
    m <- optim( fn=f , gr=g , par=beta, control=list(trace=0) , method="BFGS")
    m$par

}

gamma_init <- function(k) {

    Xbeta <- X_sh %*% beta_sh
    step <- ceiling(nrow(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, nrow(Y_sh))
    intervall <- seq.int(from = j1, to = j2)
    
    lapply(intervall,function(x){
      gamma_sh[,x] <- solveRidgeRegression( mat = V_sh, y = L_sh[x,] - Xbeta[x,],
                 beta = gamma_sh[,x],
                 epsilon = epsilon_gamma)
      return()})


}

delayed_gamma_init <- function(k) {
  
  step <- ceiling(nrow(L_sh) / children)
  j1 <- (k-1) * step + 1
  j2 <- min(k * step, nrow(L_sh))
  intervall <- seq.int(from = j1, to = j2)
  
  DelayedArray::blockApply(
    x = L_sh[intervall,],
    FUN = over_gamma,
    grid = DelayedArray::RegularArrayGrid(
      refdim = dim(L_sh[intervall,]),
      spacings = c( 1L, ncol(L_sh[intervall,]))),
    BPPARAM = BiocParallel::SerialParam(),
    V_sh = V_sh, X_sh =X_sh[intervall,,drop=FALSE], beta_sh = beta_sh,
    gamma_sh = gamma_sh, epsilon_gamma = epsilon_gamma,
    intervall = intervall)

  
}


over_gamma <- function(x, V_sh, X_sh, beta_sh, gamma_sh, epsilon_gamma,
                       intervall){
  j = attr(x,"block_id")
  Xbeta <- X_sh[j,] %*% beta_sh
  
  out <- solveRidgeRegression(mat = V_sh, y = as.vector(x - Xbeta),
                              beta = gamma_sh[,intervall[j]],
                              epsilon = epsilon_gamma)
  gamma_sh[,intervall[j]] <- out
}


beta_init <- function(k) {

    Vgamma <- t(V_sh %*% gamma_sh)
    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))
    intervall <- seq.int(from = j1, to = j2)
    
    
    lapply(intervall, function(x){
      beta_sh[,x]<-solveRidgeRegression(mat=X_sh, y=L_sh[,x] - Vgamma[,x],
      beta = beta_sh[,x],
      epsilon = epsilon_beta)
      return()})
}


delayed_beta_init <- function(k) {
  
  step <- ceiling(ncol(L_sh) / children)
  j1 <- (k-1) * step + 1
  j2 <- min(k * step, ncol(L_sh))
  intervall <- seq.int(from = j1, to = j2)
  
  DelayedArray::blockApply(
    x = L_sh[,intervall],
    FUN = over_beta,
    grid = DelayedArray::RegularArrayGrid(
      refdim = dim(L_sh[,intervall]),
      spacings = c(nrow(L_sh[,intervall]), 1L)),
    BPPARAM = BiocParallel::SerialParam(),
    X_sh = X_sh, V_sh = V_sh[intervall,,drop=FALSE], gamma_sh = gamma_sh,
    beta_sh = beta_sh, epsilon_beta = epsilon_beta,
    intervall = intervall)
  
  
}


over_beta <- function(x, X_sh, V_sh,gamma_sh, beta_sh, epsilon_beta,
                      intervall){
  j = attr(x,"block_id")

  Vgamma <- t(V_sh[j,] %*% gamma_sh)
  out <- solveRidgeRegression(mat=X_sh, y=as.vector(x - Vgamma),
                              beta = beta_sh[,intervall[j]],
                              epsilon = epsilon_beta)
  
  beta_sh[,intervall[j]] <- out
}


nb.loglik <- function(Y, mu, theta) {

    # log-probabilities of counts under the NB model
    logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    sum(logPnb)

}

delayed_ll <- function(k){
  step <- ceiling(ncol(Y_sh) / children)
  j1 <- (k-1) * step + 1
  j2 <- min(k * step, ncol(Y_sh))
  intervall <- seq.int(from = j1, to = j2)
  
  val <- sum(unlist(DelayedArray::blockApply(
    x = Y_sh[,intervall],
    FUN = over_ll,
    grid = DelayedArray::RegularArrayGrid(
      refdim = dim(Y_sh[,intervall]),
      spacings = c(nrow(Y_sh[,intervall]), 1L)),
    BPPARAM = BiocParallel::SerialParam(),
    X = X_sh, V = V_sh[intervall,,drop=FALSE], gamma = gamma_sh,
    beta = beta_sh[,intervall,drop=FALSE], W = W_sh, alpha = alpha_sh[,intervall],
    zeta = zeta_sh[intervall], intervall = intervall)))
  
  val
}

over_ll <- function(x, X, V, gamma, beta, W, alpha, zeta, intervall){
  
  j = attr(x,"block_id")
  
  theta <- exp(zeta[j])

  mu <- exp(X %*% beta[,j] + t(V[j,] %*% gamma) +
              W %*% alpha[,j])
  
  loglik <- nb.loglik(x, mu=mu, theta)
  
  loglik 
  
}
nb.regression.parseModel <- function(par, A.mu, B.mu, C.mu) {

    n <- nrow(A.mu)
    logMu <- C.mu
    dim.par <- rep(0,2)
    start.par <- rep(NA,2)
    i <- 0
  
    j <- ncol(A.mu)
    if (j>0) {
        logMu <- logMu + A.mu %*% par[(i+1):(i+j)]
        dim.par[1] <- j
        start.par[1] <- i+1
        i <- i+j
    }
  
  
    j <- ncol(B.mu)
    if (j>0) {
        logMu <- logMu + B.mu %*% par[(i+1):(i+j)]
        dim.par[2] <- j
        start.par[2] <- i+1
    }
  
    return(list(logMu=logMu, dim.par=dim.par,
                start.par=start.par))
}



nb.loglik.regression <- function(par, Y,
                                A.mu = matrix(nrow=length(Y), ncol=0),
                                B.mu = matrix(nrow=length(Y), ncol=0),
                                C.mu = matrix(0, nrow=length(Y), ncol=1),
                                C.theta = matrix(0, nrow=length(Y), ncol=1),
                                epsilon=0) {

    # Parse the model
    r <- nb.regression.parseModel(par=par,
                                A.mu = A.mu,
                                B.mu = B.mu,
                                C.mu = C.mu)
  
    # Call the log likelihood function
    z <- nb.loglik(Y, exp(r$logMu), exp(C.theta))
  
    # Penalty
    z <- z - sum(epsilon*par^2)/2
    z
}

nb.loglik.regression.gradient <- function(par, Y,
                                        A.mu = matrix(nrow=length(Y), ncol=0),
                                        B.mu = matrix(nrow=length(Y), ncol=0),
                                        C.mu = matrix(0, nrow=length(Y),
                                                        ncol=1),
                                        C.theta = matrix(0, nrow=length(Y),
                                                        ncol=1),
                                        epsilon=0) {

    # Parse the model
    r <- nb.regression.parseModel(par=par,
                                  A.mu = A.mu,
                                  B.mu = B.mu,
                                  C.mu = C.mu)
   
    
    Y=as.vector(Y)
    theta <- exp(C.theta)
    mu <- exp(r$logMu)
    n <- length(Y)
  
    # Check what we need to compute,
    # depending on the variables over which we optimize
    need.wres.mu <- r$dim.par[1] >0 || r$dim.par[2] >0
  
    # Compute the partial derivatives we need
    ## w.r.t. mu
    if (need.wres.mu) {
        wres_mu <- numeric(length = n)
        wres_mu <- Y - mu *
            (Y + theta)/(mu + theta)
        wres_mu <- as.vector(wres_mu)
    }
  
  
    # Make gradient
    grad <- numeric(0)
  
    ## w.r.t. a_mu
    if (r$dim.par[1] >0) {
        istart <- r$start.par[1]
        iend <- r$start.par[1]+r$dim.par[1]-1
        grad <- c(grad , colSums(wres_mu * A.mu) -
                    epsilon[istart:iend]*par[istart:iend])
    }
  
  
    ## w.r.t. b
    if (r$dim.par[2] >0) {
        istart <- r$start.par[2]
        iend <- r$start.par[2]+r$dim.par[2]-1
        grad <- c(grad , colSums(wres_mu * B.mu) -
                    epsilon[istart:iend]*par[istart:iend])
    }
  
    grad
}


nb.loglik.dispersion <- function(zeta, Y, mu){
  
    nb.loglik(Y, mu, exp(zeta))
}

nb.loglik.dispersion.gradient <- function(zeta, Y, mu) {
  
    theta <- exp(zeta)
    grad <- theta * (digamma(Y + theta) - digamma(theta) +
                        zeta - log(mu + theta) + 1 -
                        (Y + theta)/(mu + theta) ) 
    grad <-sum(grad)
}



optim_genwise_dispersion <- function(k, num_gene, iter) {
  
    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))
  
    if(is.null(num_gene) || iter == 1){
  
        intervall <- seq.int(from = j1, to = j2)
  
    } else {
  
        intervall <- sample(x = seq.int(from = j1, to = j2),
                        size = num_gene/children)
  
    }
    
    lapply(intervall, f_temp_d)
    
}



locfun <- function(par, Y, mu){
  nb.loglik.dispersion(zeta = par, Y, mu)
}

locgrad <- function(par, Y, mu){
  nb.loglik.dispersion.gradient(zeta = par, Y, mu)
}

f_temp_d <- function(x){
  
  
    zeta_sh[x] <- optim(fn=locfun,
                        gr=locgrad,
                        par=zeta_sh[x],
                        Y = Y_sh[,x],
                        mu = mu_sh[,x],
                        control=list(fnscale=-1,trace=0),
                        method="BFGS")$par
    return()
  }

optimr <- function(k, num_gene, iter) {

    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))

    if(is.null(num_gene) || iter == 1){

        intervall <- seq.int(from = j1, to = j2)

    } else {

        intervall <- sample(x = seq.int(from = j1, to = j2),
                          size = num_gene/children)
    }

    lapply(intervall, f_temp_r )
    
}

f_temp_r <- function(x){
  out <- optimright_fun_nb(
    beta_sh[,x,drop=FALSE], alpha_sh[,x,drop=FALSE],
    Y_sh[,x,drop=FALSE], X_sh[,,drop=FALSE],
    W_sh[,,drop=FALSE], V_sh[x,,drop=FALSE],
    gamma_sh[,, drop=FALSE], zeta_sh[x],
    nrow(Y_sh), epsilonright)
  
  params <- split_params(out, "right")
  beta_sh[,x] <- params$beta
  alpha_sh[,x] <- params$alpha
  return()
}
optimr_delayed <- function(k, num_gene) {
    
    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))
    
    if(is.null(num_gene) || iter == 1){
      
      intervall <- seq.int(from = j1, to = j2)
      
    } else {
      
      intervall <- sample(x = seq.int(from = j1, to = j2),
                          size = num_gene/children)
    }
    
    DelayedArray::blockApply(
        x = Y_sh[,intervall],
        FUN = over_optr,
        grid = DelayedArray::RegularArrayGrid(
            refdim = dim(Y_sh[,intervall]),
            spacings = c(nrow(Y_sh[,intervall]), 1L)),
        BPPARAM = BiocParallel::SerialParam(),
        beta = beta_sh[,intervall, drop = FALSE],
        alpha = alpha_sh[,intervall, drop = FALSE],X = X_sh,
        W = W_sh, V = V_sh[intervall,,drop = FALSE], gamma = gamma_sh, 
        zeta = zeta_sh[intervall], n = nrow(Y_sh),epsilonright = epsilonright,
        intervall = intervall
    )
    
  }


over_optr <- function(x, beta, alpha, X , W,
                      V , gamma , zeta,
                      n , epsilonright, intervall ){
    j = attr(x,"block_id")
    out <- optimright_fun_nb(
      beta[,j], alpha[,j], x, X,
      W, V[j,], gamma, zeta[j],
      n, epsilonright)
    
    params <- split_params(out, "right")
    beta_sh[,intervall[j]] <- params$beta
    alpha_sh[,intervall[j]] <- params$alpha
}


optimright_fun_nb <- function(beta, alpha, Y, X, W,
                              V, gamma, zeta, n, epsilonright) {
    optim( fn=nb.loglik.regression,
        gr=nb.loglik.regression.gradient,
        par=c(beta, alpha),
        Y=Y,
        A.mu=cbind(X, W),
        C.mu=t(V %*% gamma),
        C.theta=matrix(zeta, nrow =n, ncol=1),
        epsilon=epsilonright,
        control=list(fnscale = -1,trace=0),
        method="BFGS")$par
}

optiml <- function(k, num_cell, iter){

    step <- ceiling( nrow(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, nrow(Y_sh))

    if(is.null(num_cell) || iter == 1 ){

        intervall <- seq.int(from = j1, to = j2)

    } else {

        intervall <- sample(x = seq.int(from = j1, to = j2),
                            size = num_cell/children)

    }

    lapply(intervall, f_temp_l )
    
}

f_temp_l <- function(x){
  out <- optimleft_fun_nb(gamma_sh[,x],
                          W_sh[x,], Y_sh[x,] , V_sh, alpha_sh,
                          X_sh[x,], beta_sh, zeta_sh, epsilonleft)
  
  par <- split_params(out, eq = "left")
  gamma_sh[,x] <- par$gamma
  W_sh[x,] <- par$W
  return()
}
optiml_delayed <- function(k, num_cell){
    
    step <- ceiling( nrow(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, nrow(Y_sh))
    
    if(is.null(num_cell)){
      
        intervall <- seq.int(from = j1, to = j2)
      
    } else {
      
        intervall <- sample(x = seq.int(from = j1, to = j2),
                        size = num_cell/children)
      
    }
    
    DelayedArray::blockApply(
        x = Y_sh[intervall,],
        FUN = over_optl,
        grid = DelayedArray::RegularArrayGrid(
          refdim = dim(Y_sh[intervall,]),
          spacings = c( 1L, ncol(Y_sh[intervall,]))),
        BPPARAM = BiocParallel::SerialParam(),
        gamma = gamma_sh[,intervall, drop = FALSE],
        W = W_sh[intervall,, drop = FALSE], V = V_sh,
        alpha = alpha_sh, X = X_sh[intervall,, drop = FALSE], beta = beta_sh,  
        zeta = zeta_sh, epsilonleft = epsilonleft,
        intervall = intervall)
    
  }

over_optl <- function(x, gamma,  W, V , alpha,
                      X ,beta, zeta ,
                      epsilonleft, intervall){
    i = attr(x,"block_id")
    out <- optimleft_fun_nb(gamma[,i],
                            W[i,], x , V, alpha,
                            X[i,], beta, zeta, epsilonleft)
    
    params <- split_params(out, eq = "left")
    gamma_sh[,intervall[i]] <- params$gamma
    W_sh[intervall[i],] <- params$W
}


optimleft_fun_nb <- function(gamma, W, Y, V, alpha,
                             X, beta, zeta, epsilonleft) {
    optim( fn=nb.loglik.regression,
        gr=nb.loglik.regression.gradient,
        par=c(gamma, t(W)),
        Y=t(Y),
        A.mu=V,
        B.mu=t(alpha),
        C.mu=t(X%*%beta),
        C.theta=zeta,
        epsilon=epsilonleft,
        control=list(fnscale = -1,trace=0),
        method="BFGS")$par
}


split_params <- function(merged, eq = NA) {

    if(eq == "right"){
        beta_num <- nrow(beta_sh)
        alpha_num <- nrow(alpha_sh)
    
        list(
            beta = merged[seq.int(from = 1, length.out = beta_num)],
            alpha = merged[seq.int(from = beta_num+1, length.out = alpha_num)]
        )
    } else {
        gamma_num <- nrow(gamma_sh)
        W_num <- ncol(W_sh)
    
        list(
            gamma = merged[seq.int(from = 1, length.out = gamma_num)],
            W = merged[seq.int(from = gamma_num+1, length.out = W_num)]
        )
    }
}

