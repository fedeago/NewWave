solveRidgeRegression <- function(x, y, beta=rep(0,NCOL(x)), epsilon=1e-6) {

  # loglik
  f <-function(b) {
    eta <- x %*% b
    l <- sum((y-eta)^2)/2
    l + sum(epsilon*b^2)/2
  }

  # gradient of loglik
  g <-function(b) {
    eta <- x %*% b
    l <- t(x) %*% (eta - y)
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
  
  for (j in intervall) {
    out <- solveRidgeRegression(x = V_sh, y = L_sh[j,] - Xbeta[j,],
                                 beta = gamma_sh[,j],
                                epsilon = epsilon_gamma)

    gamma_sh[,j] <- out
  }

}

beta_init <- function(k) {

  Vgamma <- t(V_sh %*% gamma_sh)
  step <- ceiling(ncol(Y_sh) / children)
  j1 <- (k-1) * step + 1
  j2 <- min(k * step, ncol(Y_sh))
  intervall <- seq.int(from = j1, to = j2)
  
  for (j in intervall) {
    out <- solveRidgeRegression(x=X_sh, y=L_sh[,j] - Vgamma[, j],
                                 beta = beta_sh[,j],
                                epsilon = epsilon_beta)


    beta_sh[,j] <- out
  }
}

nb.loglik <- function(Y, mu, theta) {
  
  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  sum(logPnb)
  
}






nb.regression.parseModel <- function(alpha, A.mu, B.mu, C.mu) {

  n <- nrow(A.mu)
  logMu <- C.mu
  dim.alpha <- rep(0,2)
  start.alpha <- rep(NA,2)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
    start.alpha[1] <- i+1
    i <- i+j
  }


  j <- ncol(B.mu)
  if (j>0) {
    logMu <- logMu + B.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[2] <- j
    start.alpha[2] <- i+1
  }

  return(list(logMu=logMu, dim.alpha=dim.alpha,
              start.alpha=start.alpha))
}


nb.loglik.regression <- function(alpha, Y,
                                 A.mu = matrix(nrow=length(Y), ncol=0),
                                 B.mu = matrix(nrow=length(Y), ncol=0),
                                 C.mu = matrix(0, nrow=length(Y), ncol=1),
                                 C.theta = matrix(0, nrow=length(Y), ncol=1),
                                 epsilon=0) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                A.mu = A.mu,
                                B.mu = B.mu,
                                C.mu = C.mu)

  # Call the log likelihood function
  z <- nb.loglik(Y, exp(r$logMu), exp(C.theta))

  # Penalty
  z <- z - sum(epsilon*alpha^2)/2
  z
}

nb.loglik.regression.gradient <- function(alpha, Y,
                                          A.mu = matrix(nrow=length(Y), ncol=0),
                                          B.mu = matrix(nrow=length(Y), ncol=0),
                                          C.mu = matrix(0, nrow=length(Y),
                                                        ncol=1),
                                          C.theta = matrix(0, nrow=length(Y),
                                                           ncol=1),
                                          epsilon=0) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                A.mu = A.mu,
                                B.mu = B.mu,
                                C.mu = C.mu)
  Y=as.vector(Y)
  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)

  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0 || r$dim.alpha[2] >0

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
  if (r$dim.alpha[1] >0) {
    istart <- r$start.alpha[1]
    iend <- r$start.alpha[1]+r$dim.alpha[1]-1
    grad <- c(grad , colSums(wres_mu * A.mu) -
                epsilon[istart:iend]*alpha[istart:iend])
  }


  ## w.r.t. b
  if (r$dim.alpha[2] >0) {
    istart <- r$start.alpha[2]
    iend <- r$start.alpha[2]+r$dim.alpha[2]-1
    grad <- c(grad , colSums(wres_mu * B.mu) -
                epsilon[istart:iend]*alpha[istart:iend])
  }

  grad
}


nb.loglik.dispersion <- function(zeta, Y, mu){
  
  # zeta <- rbind(zeta)[rep(1,nrow(Y)),]
  nb.loglik(Y, mu, exp(zeta))
  
}

nb.loglik.dispersion.gradient <- function(zeta, Y, mu) {
  
  theta <- exp(zeta)
  
  grad <- theta * (digamma(Y + theta) - digamma(theta) +
                         log(theta) - log(mu + theta) + 1 -
                         (Y + theta)/(mu + theta) ) 
  grad <-sum(grad)
}



optim_genwise_dispersion <- function(k, num_gene, iter) {
  
  locfun <- function(zeta, Y, mu){
    nb.loglik.dispersion(zeta, Y, mu)
  }
  
  locgrad <- function(zeta, Y, mu){
    nb.loglik.dispersion.gradient(zeta, Y, mu)
  }
  
  step <- ceiling(ncol(Y_sh) / children)
  j1 <- (k-1) * step + 1
  j2 <- min(k * step, ncol(Y_sh))

  if(is.null(num_gene)){

    intervall <- seq.int(from = j1, to = j2)

  } else {

    intervall <- sample(x = seq.int(from = j1, to = j2), size = num_gene/children)

  }
  
  out <- list()
  
  for (j in intervall){
    
    out[[j]] <- optim(fn=locfun , gr=locgrad,
                      par=zeta_sh[j],Y = Y_sh[,j], mu = mu[,j],
                      control=list(fnscale=-1,trace=0), method="BFGS")
    zeta_sh[j] <- out[[j]]$par
  }
  
  # likelihood <- t(matrix(unlist(out), nrow = length(out[[length(out)]])))
  # code <- Sys.getpid()
  # filepath <- paste0("~/Scrivania/prove_zinb_par/like_genewise_dispersion/likelihood_",iter,"_",code,".RDS")
  # saveRDS(likelihood[,2], file = filepath )
  # 
  
}

optimr <- function(k, num_gene) {

    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))

    if(is.null(num_gene)){

      intervall <- seq.int(from = j1, to = j2)

    } else {

      intervall <- sample(x = seq.int(from = j1, to = j2), size = num_gene/children)
    }

    
    for ( j in intervall){

      out <- optimright_fun_nb(
        beta_sh[,j], alpha_sh[,j], Y_sh[,j], X_sh,
        W_sh, V_sh[j,], gamma_sh, zeta_sh[j],
        nrow(Y_sh), epsilonright)

      params <- split_params(out, "right")
      beta_sh[,j] <- params$beta
      alpha_sh[,j] <- params$alpha


  }


}

optimr_delayed <- function(k, num_gene) {
    
    step <- ceiling(ncol(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, ncol(Y_sh))
    
    if(is.null(num_gene)){
      
      intervall <- seq.int(from = j1, to = j2)
      
    } else {
      
      intervall <- sample(x = seq.int(from = j1, to = j2), size = num_gene/children)
    }
    
    blockApply(
      x = Y_sh[,intervall],
      FUN = over_optr,
      grid = RegularArrayGrid(
        refdim = dim(Y_sh[,intervall]),
        spacings = c(nrow(Y_sh[,intervall]), 1L)),
      BPPARAM = SerialParam(),
      beta = beta_sh[,intervall, drop = F], alpha = alpha_sh[,intervall, drop = F], X = X_sh,
      W = W_sh, V = V_sh[intervall,,drop = F], gamma = gamma_sh, 
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
         C.theta=matrix(zeta, nrow = n, ncol = 1),
         epsilon=epsilonright,
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}

optiml <- function(k, num_cell){

    step <- ceiling( nrow(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, nrow(Y_sh))

    if(is.null(num_cell)){

      intervall <- seq.int(from = j1, to = j2)

    } else {

      intervall <- sample(x = seq.int(from = j1, to = j2), size = num_cell/children)

    }

    for (i in intervall){
      out <- optimleft_fun_nb(gamma_sh[,i],
                              W_sh[i,], Y_sh[i,] , V_sh, alpha_sh,
                              X_sh[i,], beta_sh, zeta_sh, epsilonleft)

      params <- split_params(out, eq = "left")
      gamma_sh[,i] <- params$gamma
      W_sh[i,] <- params$W
      
      
    }
    
}


optiml_delayed <- function(k, num_cell){
    
    step <- ceiling( nrow(Y_sh) / children)
    j1 <- (k-1) * step + 1
    j2 <- min(k * step, nrow(Y_sh))
    
    if(is.null(num_cell)){
      
      intervall <- seq.int(from = j1, to = j2)
      
    } else {
      
      intervall <- sample(x = seq.int(from = j1, to = j2), size = num_cell/children)
      
    }
    
    blockApply(
      x = Y_sh[intervall,],
      FUN = over_optl,
      grid = RegularArrayGrid(
        refdim = dim(Y_sh[intervall,]),
        spacings = c( 1L, ncol(Y_sh[intervall,]))),
      BPPARAM = SerialParam(),
      gamma = gamma_sh[,intervall, drop = F],  W = W_sh[intervall,, drop = F], V = V_sh,
      alpha = alpha_sh, X = X_sh[intervall,, drop = F], beta = beta_sh,  
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
         control=list(fnscale=-1,trace=0),
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

