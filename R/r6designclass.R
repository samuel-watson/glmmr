# DESIGN CLASS
Design <- R6::R6Class("Design",
                  public = list(
                    covariance = NULL,
                    mean_function = NULL,
                    exp_cond = NULL,
                    Sigma = NULL,
                    var_par = NULL,
                    initialize = function(covariance,
                                          mean.function,
                                          var_par = NULL,
                                          verbose=TRUE){
                      if(is(covariance,"R6")){
                        if(is(covariance,"Covariance")){
                          self$covariance <- covariance
                        } else {
                          stop("covariance should be Covariance class or list of appropriate arguments")
                        }
                      } else if(is(covariance,"list")){
                        self$covariance <- Covariance$new(
                          formula= covariance$formula,
                          data = covariance$data,
                          parameters = covariance$parameters,
                          verbose = verbose
                        )
                      }

                      if(is(mean.function,"R6")){
                        if(is(mean.function,"MeanFunction")){
                          self$mean_function <- mean.function
                        } else {
                          stop("mean.function should be MeanFunction class or list of appropriate arguments")
                        }
                      } else if(is(mean.function,"list")){
                        self$mean_function <- MeanFunction$new(
                          formula = mean.function$formula,
                          data = mean.function$data,
                          family = mean.function$family,
                          parameters = mean.function$parameters,
                          random_function = mean.function$random_function,
                          treat_par = mean.function$treat_par,
                          verbose = verbose
                        )
                      }



                      self$var_par <- var_par

                      self$generate()
                      private$hash <- private$hash_do()
                    },
                    print = function(){
                      cat("\n----------------------------------------\n")
                      print(self$mean_function)
                      cat("\n----------------------------------------\n")
                      print(self$covariance)
                      cat("\n----------------------------------------\n")
                    },
                    n = function(){
                      self$mean_function$n()
                    },
                    generate = function(){
                      # add check for var par with gaussian family

                      private$genW(family = self$mean_function$family,
                                   Xb = self$mean_function$.__enclos_env__$private$Xb,
                                   var_par = self$var_par)
                      private$genS(D = self$covariance$D,
                                   Z = self$covariance$Z,
                                   W = private$W)
                    },
                    power = function(par,
                                     value,
                                     alpha=0.05,
                                     type="linear",
                                     method=NULL,
                                     iter=10,
                                     skip.check = FALSE,
                                     parallel=TRUE){
                      if(!skip.check)self$check(verbose=FALSE)
                      if(missing(par)|missing(value))stop("parameter missing")
                      if(type=="linear"){
                        old_par <- self$mean_function$parameters[[par]]
                        self$mean_function$parameters[[par]] <- value
                        self$check(verbose=FALSE)
                        M <- private$information_matrix()
                        v0 <- solve(M)[par,par]
                        pwr <- pnorm(value/(sqrt(v0)) - qnorm(1-alpha/2))
                        self$mean_function$parameters[[par]] <- old_par
                      } else if(type=="sim"){
                        return(NULL)
                        # if(parallel){
                        #   cl <- parallel::makeCluster(parallel::detectCores()-2)
                        #   parallel::clusterEvalQ(cl,require(Matrix))
                        #   parallel::clusterEvalQ(cl,require(lme4))
                        #   parallel::clusterEvalQ(cl,require(minqa))
                        #  ests <- pbapply::pbreplicate(iter, private$lme_est(par=par,
                        #                                                     value=value),
                        #                               cl=cl)
                        #  parallel::stopCluster(cl)
                        # } else {
                        #   ests <- pbapply::pbreplicate(iter, private$lme_est(par=par,
                        #                                                      value=value))
                        # }
                        #  
                        #  ests <- apply(ests,2,function(x)x$par/x$SE)
                        #  
                        #  ests <- ests[par,]
                        #  ests <- (1-pnorm(abs(ests)))*2
                        #  pwr <- mean(I(ests < alpha),na.rm=TRUE)
                        #  if(any(is.na(ests)))message(paste0("failure to converge in ",length(which(is.na(ests)))," models"))
                      }
                      return(pwr)
                    },
                    subset_rows = function(index){
                      self$mean_function$subset_rows(index)
                      self$covariance$subset(index)
                    },
                    subset_cols = function(index){
                      self$mean_function$subset_cols(index)
                    },
                    plot = function(x,
                                    y,
                                    z=NULL,
                                    treat){
                      if(is.null(z)){
                        ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          geom_count(aes(color=self$mean_function$data[,treat]))+
                          theme_bw()+
                          theme(panel.grid=element_blank())+
                          scale_color_viridis_c(name=treat)+
                          scale_size_area()
                      } else {
                        ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          geom_count(aes(color=self$mean_function$data[,treat]))+
                          facet_wrap(~.data[[z]])+
                          theme_bw()+
                          theme(panel.grid=element_blank())+
                          scale_color_viridis_c(name=treat)+
                          scale_size_area()
                      }},
                    sim_data = function(type = "y",
                                        par=NULL,
                                        value=NULL){
                      if(!is.null(par)){
                        if(is.null(value))stop("set parameter value")
                        orig_par <- self$mean_function$parameters[par]
                        self$mean_function$parameters[par] <- value
                        self$mean_function$check(verbose = FALSE)
                      }
                      re <- MASS::mvrnorm(n=1,mu=rep(0,nrow(self$covariance$D)),Sigma = self$covariance$D)
                      mu <- c(drop(self$mean_function$.__enclos_env__$private$Xb)) + as.matrix(self$covariance$Z%*%re)
                      
                      f <- self$mean_function$family
                      if(f[1]=="poisson"){
                        if(f[2]=="log"){
                          y <- rpois(self$n(),exp(mu))
                        }
                        if(f[2]=="identity"){
                          y <- rpois(self$n(),mu)
                        }
                      }

                      if(f[1]=="binomial"){
                        if(f[2]=="logit"){
                          y <- rbinom(self$n(),1,exp(mu)/(1+exp(mu)))
                        }
                        if(f[2]=="log"){
                          y <- rbinom(self$n(),1,exp(mu))
                        }
                        if(f[2]=="identity"){
                          y <- rbinom(self$n(),1,mu)
                        }
                        if(f[2]=="probit"){
                          y <- rbinom(self$n(),1,pnorm(mu))
                        }
                      }

                      if(f[1]=="gaussian"){
                        if(f[2]=="identity"){
                          if(is.null(self$var_par))stop("For gaussian(link='identity') provide var_par")
                          y <- rnorm(self$n(),mu,self$var_par)
                        }
                        if(f[2]=="log"){
                          if(is.null(self$var_par))stop("For gaussian(link='log') provide var_par")
                          #CHECK THIS IS RIGHT
                          y <- rnorm(self$n(),exp(mu),self$var_par)
                        }
                      }

                      if(f[1]=="gamma"){
                        if(f[2]=="inverse"){
                          if(is.null(self$var_par))stop("For gamma(link='inverse') provide var_par")
                          #CHECK THIS IS RIGHT
                          y <- rgamma(self$n(),shape = 1/(mu*self$var_par),rate = 1/self$var_par)
                        }
                      }
                      if(!missing(par)){
                        self$mean_function$parameters[par]<-orig_par
                        self$check(verbose = FALSE)
                      }
                      if(type=="y")return(y)
                      if(type=="data.frame")return(cbind(y,self$data$data,self$covariance$location))
                    },
                    check = function(verbose=TRUE){
                      self$covariance$check(verbose=verbose)
                      self$mean_function$check(verbose = verbose)
                      if(private$hash != private$hash_do()){
                        self$generate()
                      }
                    },
                    apv = function(prior,
                                   var,
                                   prior.fun,
                                   iter,
                                   verbose=TRUE){
                      if(verbose)message("Monte Carlo integration")
                      samps <- pbapply::pbreplicate(iter,self$posterior(prior,var,do.call(prior.fun,list())))
                      summary(samps)
                    },
                    MCNR = function(y,
                                            start=c(0,0,1,0.1),
                                            verbose=TRUE,
                                            tol = 1e-2,
                                            m=100,
                                            bSEonly=TRUE,
                                            use_cmdstanr=TRUE){
                      
                      theta <- start
                      thetanew <- rep(1,length(theta))
                      if(verbose)cat("\nStart: ",start,"\n")
                      iter <- 0
                      P <- ncol(self$mean_function$X)
                      parInds <- list(b = c(1:ncol(self$mean_function$X)),
                                      sig = (ncol(self$mean_function$X)+1),
                                      cov = c((ncol(self$mean_function$X)+2):(ncol(self$mean_function$X)+1+
                                                                                length(unlist(self$covariance$parameters)))))
                      
                      niter <- m
                      P = ncol(self$mean_function$X)
                      Q = ncol(self$covariance$Z)
                      #covFix <- FALSE
                      
                      D_func <- function(pars){
                        self$covariance$parameters <- relist(self$covariance$parameters,
                                                             pars)[[1]]
                        self$check(verbose=FALSE)
                        logdetD <- 2*sum(log(Matrix::diag(Matrix::chol(self$covariance$D))))
                        invD <- Matrix::solve(self$covariance$D)
                        Q <- dim(dsamps)[3]
                        
                        dmv <- sapply(1:dim(dsamps)[1],function(i)log_mvnd(as.matrix(dsamps[i,1,]),logdetD,invD))
                        -mean(dmv)
                      }
                      
                      #parse family
                      file_type <- mcnr_family(self$mean_function$family)
                      invfunc <- self$mean_function$family$linkinv
                      
                      ## set up sampler
                      if(use_cmdstanr){
                        if(!requireNamespace("cmdstanr"))stop("cmdstanr not available")
                        model_file <- system.file("stan",
                                                  file_type$file,
                                                  package = "glmmr",
                                                  mustWork = TRUE)
                        mod <- cmdstanr::cmdstan_model(model_file)
                        
                      }
                      while(any(abs(theta-thetanew)>tol)){
                        iter <- iter + 1
                        if(verbose)cat("\nIter: ",iter,": ")
                        thetanew <- theta
                        
                        
                        Xb <- Matrix::drop(self$mean_function$X %*% thetanew[parInds$b])
                        
                        data <- list(
                          N = self$n(),
                          P = P,
                          Q = Q,
                          Xb = Xb,
                          L = as.matrix(Matrix::t(Matrix::chol(self$covariance$D))),
                          Z = as.matrix(self$covariance$Z),
                          y = y,
                          sigma = thetanew[parInds$sig],
                          type=as.numeric(file_type$type)
                        )
                        
                        if(use_cmdstanr){
                          capture.output(fit <- mod$sample(data = data,
                                                           chains = 1,
                                                           iter_warmup = 500,
                                                           iter_sampling = m,
                                                           refresh = 0),
                                         file=tempfile())
                          dsamps <- fit$draws("gamma")
                        } else {
                          capture.output(fit <- rstan::sampling(stanmodels[[gsub(".stan","",file_type$file)]],
                                                                chains = 1,
                                                                iter_warmup = 500,
                                                                iter_sampling = m))
                          dsamps <- extract(fit,"gamma")
                        }
                        
                        
                        
                        
                        resids <- sapply(1:niter, function(i) y-invfunc(Matrix::drop(Xb - self$covariance$Z%*%as.matrix(dsamps[i,1,]))))
                        
                        
                        sigmas <- c(apply(resids,2,sd))
                        sigmahat <- mean(sigmas)
                        XtWXlist <- lapply(1:niter,function(i)Matrix::crossprod(self$mean_function$X,Matrix::diag(sigmas[i],self$n()))%*%
                                             self$mean_function$X)
                        EXtWX <- Reduce(`+`,XtWXlist)/niter
                        Wulist <- lapply(1:niter,function(i)Matrix::diag(sigmas[i],self$n())%*%resids[,i])
                        EWu <- Reduce(`+`,Wulist)/niter
                        
                        theta[parInds$b] <-  theta[parInds$b] + Matrix::drop(Matrix::solve(EXtWX)%*%Matrix::t(self$mean_function$X)%*%EWu)
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
                      
                      #need to rewrite L_func for each families likelihood!
                      L_func <- function(pars){
                        niter <- dim(dsamps)[1]
                        dlist <- sapply(1:niter,function(i)self$covariance$Z %*% Matrix::Matrix(dsamps[i,1,]))
                        zd <- Matrix::drop(Reduce(rbind,dlist))
                        mu <- rep(Matrix::drop(self$mean_function$X %*% pars[1:P]),niter) + zd
                        mu <- invfunc(mu)
                        n <- self$n()
                        # change density for the model
                        if(self$mean_function$family[[1]]=="gaussian"){
                          lf <- log(dnorm(rep(y,niter),mu,theta[parInds$sig]))
                        } else if(self$mean_function$family[[1]]=="binomial"){
                          lf <- log(dbinom(rep(y,niter),1,mu))
                        }
                        lfa <- aggregate(lf,list(rep(1:niter,each=n)),sum)
                        -mean(lfa$x)
                      }
                      
                      log_lik <- function(pars){
                        l1 <- L_func(pars[parInds$b])
                        l2 <- D_func(pars[parInds$cov])
                        -l1-l2
                      }
                      
                      #FULL HESSIAN WON'T WORK FOR GLM MODELS AS SIGMA NOT A PARAMETER IN THE MODEL
                      # HAVE REMOVED FROM SE EVALUATION FOR ALL MODELS FOR NOW
                      # BUT BETTER TO REPORT ONLY FOR GAUSSIAN MODELS
                      
                      if(bSEonly){
                        hess <- pracma::hessian(function(x)-L_func(x),theta[c(parInds$b)])
                        SE <- sqrt(Matrix::diag(Matrix::solve(-hess)))
                        res <- data.frame(par = theta,SE=c(SE,rep(NA,length(parInds$cov)+1)))
                      } else {
                        hess <- pracma::hessian(log_lik,theta[c(parInds$b,parInds$cov)])
                        SE <- sqrt(Matrix::diag(Matrix::solve(-hess)))
                        res <- data.frame(par = theta[c(parInds$b,parInds$cov,parInds$sig)],SE=c(SE,NA))
                      }
                      
                      
                      return(res)
                    },
                    posterior = function(prior,
                                         var,
                                         parameters){
                      #move to private and set this as Monte Carlo integration
                      #can just request a function that outputs a new set of covariance parameters
                      R <- solve(Matrix::Matrix(diag(prior)))
                      S <- private$genS(self$covariance$sampleD(parameters),self$covariance$Z,private$W,update=FALSE)
                      M <- R + Matrix::crossprod(self$mean_function$X,solve(S))%*%self$mean_function$X
                      M <- solve(M)
                      M[var,var]
                    }
                  ),
                  private = list(
                    W = NULL,
                    Xb = NULL,
                    logit = function(x){
                      exp(x)/(1+exp(x))
                    },
                    genW = function(family,
                                    Xb,
                                    var_par=NULL){
                      # assume random effects value is at zero
                      f <- family
                      Xb <- c(Xb)
                      if(!f[1]%in%c("poisson","binomial","gaussian","gamma"))stop("family must be one of Poisson, Binomial, Gaussian, Gamma")

                      if(f[1]=="poisson"){
                        if(f[2]=="log"){
                          W <- diag(1/(exp(Xb)))
                        }
                        if(f[2]=="identity"){
                          W <- diag(exp(Xb))
                        }
                      }

                      if(f[1]=="binomial"){
                        if(f[2]=="logit"){
                          W <- diag(1/(private$logit(Xb)*(1-private$logit(Xb))))
                        }
                        if(f[2]=="log"){
                          W <- diag((1-private$logit(Xb))/(private$logit(Xb)))
                        }
                        if(f[2]=="identity"){
                          W <- diag((private$logit(Xb)*(1-private$logit(Xb))))
                        }
                        if(f[2]=="probit"){
                          W <- diag((pnorm(Xb)*(1-pnorm(Xb)))/(dnorm(Xb)))
                        }
                      }

                      if(f[1]=="gaussian"){
                        if(f[2]=="identity"){
                          if(is.null(var_par))stop("For gaussian(link='identity') provide var_par")
                          W <- var_par*diag(length(Xb))
                        }
                        if(f[2]=="log"){
                          if(is.null(var_par))stop("For gaussian(link='log') provide var_par")
                          W <- diag(var_par/exp(Xb))
                        }
                      }

                      if(f[1]=="gamma"){
                        if(f[2]=="inverse"){
                          if(is.null(var_par))stop("For gamma(link='inverse') provide var_par")
                          W <- var_par*diag(length(Xb))
                        }
                      }
                      private$W <- Matrix::Matrix(W)
                    },
                    genS = function(D,Z,W,update=TRUE){
                      if(is(D,"numeric")){
                        S <- W + D * Matrix::tcrossprod(Z)
                      } else {
                        S <- W + Z %*% Matrix::tcrossprod(D,Z)
                      }
                      if(update){
                        self$Sigma <- Matrix::Matrix(S)
                        private$hash <- private$hash_do()
                      } else {
                        return(S)
                      }

                    },
                    hash = NULL,
                    hash_do = function(){
                      digest::digest(c(self$covariance$.__enclos_env__$private$hash,
                                       self$mean_function$.__enclos_env__$private$hash))
                    },
                    lme_est = function(par=NULL,
                                       value=NULL){
                      return(NULL)
                      # lambda <-  Matrix::t(Matrix::chol(self$covariance$D))#,nrow=nrow(self$covariance$D)))
                      # ncovpar <- length(unlist(self$covariance$parameters))
                      # nmfpar <- ncol(self$mean_function$X)
                      # parInds <- list(covar = c(1:ncovpar),
                      #                 fixef = c((ncovpar+1):(ncovpar+nmfpar)))
                      # ysim <- self$sim_data(par=par,
                      #                       value=value)
                      # dat1 <- cbind(ysim,self$mean_function$data)
                      # print(summary(lme4::glmer(ysim~int+(1|cl),data=dat1,family=self$mean_function$family)))
                      # 
                      # orig_pars <- self$covariance$parameters
                      # 
                      # pp <- lme4::merPredD$new(
                      #   X = self$mean_function$X,
                      #   Zt = Matrix::t(self$covariance$Z),
                      #   Lambdat = lambda,
                      #   Lind = seq_along(lambda@x),
                      #   theta = lambda@x,
                      #   n = des$n())
                      # 
                      # resp <- lme4::glmResp$new(
                      #   y = ysim,
                      #   family = self$mean_function$family)
                      # 
                      # updateTheta <- function(pars){
                      #   # self$covariance$parameters <- relist(self$covariance$parameters,
                      #   #                           value = pars)[[1]]
                      #   # self$covariance$check(verbose = FALSE)
                      #   # newD <- self$covariance$D
                      #   # newD <- (resp$resDev()/self$n())*(pp$Lambdat%*%Matrix::chol2inv(pp$L())%*%Matrix::t(pp$Lambdat))
                      #   # cholD <- tryCatch(Matrix::chol(newD),error=function(e)NULL)
                      #   # if(is.null(cholD)){
                      #   #   newD <- Matrix::nearPD(newD)
                      #   #   cholD <- tryCatch(Matrix::chol(newD),error=function(e)print("help!"))
                      #   # }
                      #   # Matrix::t(cholD)@x
                      #   rep(pars,8)
                      # }
                      # 
                      # devfun <- function(pars) {
                      #   resp$setOffset(rep(0,self$n()))
                      #   resp$updateMu(rep(0,self$n()))
                      #   pp$setTheta(as.double(updateTheta(pars[parInds$covar])))
                      #   spars <- as.numeric(pars[parInds$fixef])
                      #   resp$setOffset( pp$X %*% spars)
                      #   p <- lme4::glmerLaplaceHandle(pp$ptr(), resp$ptr(), 1, 1e-6, 30, TRUE)
                      #   resp$updateWts()
                      #   p
                      # }
                      # 
                      # opt <- minqa::bobyqa(par = c(rep(0.01,ncovpar),rep(0.1,nmfpar)),
                      #                      fn = devfun,
                      #                      lower = c(rep(1e-5,ncovpar),rep(-Inf,nmfpar)))
                      # 
                      # 
                      # # dev1 <- pirls(X = self$mean_function$X,
                      # #               y=ysim,
                      # #               Zt = Matrix::t(self$covariance$Z),
                      # #               Lambdat = lambda,
                      # #               thfun = updateTheta,
                      # #               theta = runif(ncovpar,0.01,0.1),
                      # #               weights = rep(1,self$n()),
                      # #               offset=rep(0,self$n()),
                      # #               eta=numeric(self$n()),
                      # #               family=self$mean_function$family)
                      # sigma <- (resp$resDev()/self$n())
                      # beta <- as.numeric(opt$par[parInds$fixef])
                      # print(sqrt((resp$resDev()/self$n())*diag(pp$unsc())))
                      # print((resp$resDev()/self$n())*(pp$Lambdat%*%Matrix::crossprod(Matrix::solve(pp$L()))%*%Matrix::t(pp$Lambdat)))
                      # lambda@x <- rep(opt$par[parInds$covar],8)
                      # print(Matrix::tcrossprod(lambda))
                      # 
                      # if(is.null(opt)){
                      #   return(data.frame(b=rep(NA,nmfpar),se=rep(NA,nmfpar)))
                      # } else {
                      #   
                      #   se <- sqrt(diag(solve(Matrix::crossprod(self$mean_function$X,solve(self$Sigma))%*%self$mean_function$X )))
                      #   b <- opt$par[parInds$fixef]
                      #   theta <- c(opt$par[parInds$covar],opt$par[parInds$covar]*(resp$resDev()/self$n()))
                      #   self$covariance$parameters <- orig_pars
                      #   self$covariance$check(verbose=FALSE)
                      #   self$check(verbose=FALSE)
                      #   
                      #   return(data.frame(par= c(b,theta),SE = c(se,rep(NA,ncovpar*2))))
                      # }

                    },
                    information_matrix = function(){
                      Matrix::crossprod(self$mean_function$X,solve(self$Sigma))%*%self$mean_function$X
                    }
                  ))


relist <- function(lst,value,p=0){
  if(is(lst,"list")){
    for(i in 1:length(lst)){
      out <- Recall(lst[[i]],value,p=p)
      lst[[i]] <- out[[1]]
      p <- out[[2]]
    }
  } else {
    for(i in 1:length(lst)){
      lst[i] <- value[p+1]
      p <- p + 1
    }
  }
  return(list(lst,p))
}


