# Markov chain monte carlo Newton Raphson algorithm
# only works for linear models currently
# requires the stan file mcml.stan
# I will add stan functionality and compiled models
# for now this will only work locally

MCNR <- function(des,
                 y,
                 start=c(0,0,1,0.1),
                 verbose=TRUE,
                 tol = 1e-2,
                 cov_tol = 1e-2,
                 m=500,
                 bSEonly=TRUE){
  
  mod <- cmdstan_model("C:/docs/glmmr-dev/mcml.stan")
  theta <- start
  thetanew <- rep(1,length(theta))
  if(verbose)cat("\nStart: ",start,"\n")
  iter <- 0
  P <- ncol(des$mean_function$X)
  parInds <- list(b = c(1:ncol(des$mean_function$X)),
                  sig = (ncol(des$mean_function$X)+1),
                  cov = c((ncol(des$mean_function$X)+2):(ncol(des$mean_function$X)+1+
                                                           length(unlist(des$covariance$parameters)))))
  
  niter <- m
  P = ncol(des$mean_function$X)
  Q = ncol(des$covariance$Z)
  #covFix <- FALSE
  
  D_func <- function(pars){
    des$covariance$parameters <- relist(des$covariance$parameters,
                                        pars)[[1]]
    des$check(verbose=FALSE)
    logdetD <- 2*sum(log(Matrix::diag(Matrix::chol(des$covariance$D))))
    invD <- Matrix::solve(des$covariance$D)
    Q <- dim(dsamps)[3]
    
    dmv <- sapply(1:dim(dsamps)[1],function(i)log_mvnd(as.matrix(dsamps[i,1,]),logdetD,invD))
    -mean(dmv)
  }
  
  while(any(abs(theta-thetanew)>tol)){
    iter <- iter + 1
    if(verbose)cat("\nIter: ",iter,": ")
    thetanew <- theta
    
    
    Xb <- Matrix::drop(des$mean_function$X %*% thetanew[parInds$b])
    
    data <- list(
      N = des$n(),
      P = P,
      Q = Q,
      Xb = Xb,
      L = as.matrix(Matrix::t(Matrix::chol(des$covariance$D))),
      Z = as.matrix(des$covariance$Z),
      y = y,
      sigma = thetanew[parInds$sig]
    )
    
    fit <- mod$sample(data = data,
                      chains = 1,
                      iter_warmup = 500,
                      iter_sampling = m,
                      refresh = 0
    )
    
    dsamps <- fit$draws("gamma")
    resids <- sapply(1:niter, function(i) y-Xb - des$covariance$Z%*%as.matrix(dsamps[i,1,]))
    sigmahat <- mean(unlist(lapply(resids,sd)))
    sigmas <- unlist(lapply(resids,sd))
    
    XtWXlist <- lapply(1:niter,function(i)Matrix::crossprod(des$mean_function$X,Matrix::diag(sigmas[i],des$n()))%*%
                         des$mean_function$X)
    EXtWX <- Reduce(`+`,XtWXlist)/niter
    Wulist <- lapply(1:niter,function(i)Matrix::diag(sigmas[i],des$n())%*%resids[[i]])
    EWu <- Reduce(`+`,Wulist)/niter
    
    theta[parInds$b] <-  theta[parInds$b] + Matrix::drop(Matrix::solve(EXtWX)%*%Matrix::t(des$mean_function$X)%*%EWu)
    theta[parInds$sig] <- sigmahat
    
    opts <-  minqa::bobyqa(theta[parInds$cov],
                           fn = D_func,
                           lower = rep(1e-6,length(parInds$cov)),
                           control = list(rhobeg = 0.02,
                                          rhoend = 0.02*tol))
    
    theta[parInds$cov] <- opts$par
    
    #if(abs(theta[parInds$cov]-thetanew[parInds$cov])<cov_tol)covFix <- TRUE
    
    if(verbose)cat("\ntheta:",theta)
  }
  
  L_func <- function(pars){
    niter <- dim(dsamps)[1]
    dlist <- sapply(1:niter,function(i)des$covariance$Z %*% Matrix(dsamps[i,1,]))
    zd <- Matrix::drop(Reduce(rbind,dlist))
    mu <- rep(Matrix::drop(des$mean_function$X %*% pars[1:P]),niter) + zd
    n <- des$n()
    lf <- log(dnorm(rep(y,niter),mu,pars[P+1]))
    lfa <- aggregate(lf,list(rep(1:niter,each=n)),sum)
    -mean(lfa$x)
  }
  
  log_lik <- function(pars){
    l1 <- L_func(c(pars[parInds$b],pars[parInds$sig]))
    l2 <- D_func(pars[parInds$cov])
    -l1-l2
  }
  
  #permutation test inference for standard errors?
  
  if(bSEonly){
    hess <- pracma::hessian(function(x)-L_func(x),theta[c(parInds$b,parInds$sig)])
    SE <- sqrt(Matrix::diag(Matrix::solve(-hess)))
    res <- data.frame(par = theta,SE=c(SE,rep(NA,length(parInds$cov))))
  } else {
    hess <- pracma::hessian(log_lik,theta)
    SE <- sqrt(Matrix::diag(Matrix::solve(-hess)))
    res <- data.frame(par = theta,SE=SE)
  }
  
  
  return(res)
}