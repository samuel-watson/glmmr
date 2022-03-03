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
                    xsq = function(x){
                      re <- MASS::mvrnorm(n=1,mu=rep(0,nrow(self$covariance$D)),Sigma = self$covariance$D)
                      mu <- c(drop(as.matrix(self$mean_function$X)%*%self$mean_function$parameters)) + as.matrix(self$covariance$Z%*%re)
                      #y <- rbinom(self$n(),1,exp(mu)/(1+exp(mu)))
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
                      return(y)
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
                    analysis = function(type,
                                        iter,
                                        par,
                                        alpha = 0.05,
                                        parallel,
                                        verbose = TRUE,
                                        digits = 2,
                                        ...){
                      
                      if(type=="sim_data"&is.null(private$sim_data))stop("no simulation data saved in object")
                      if(type=="sim"){
                        if(parallel){
                          cl <- parallel::makeCluster(parallel::detectCores()-1)
                          parallel::clusterEvalQ(cl,library(Matrix))
                          #change when package built!
                          parallel::clusterEvalQ(cl,devtools::load_all())
                          # out <- parallel::parLapply(cl,
                          #                            1:10,
                          #                            function(i)self$gen_sim_data(m=m))
                          out <- pbapply::pblapply(1:iter,
                                                   function(i)private$gen_sim_data(...),
                                                   cl=cl)
                          parallel::stopCluster(cl)
                        } else {
                          out <- pbapply::pblapply(1:iter,
                                                   function(i)private$gen_sim_data(...))
                        }
                        if(verbose)message("saving simulation data")
                        private$sim_data <- out
                      }
                      if(type=="sim_data")out <- private$sim_data
                      
                      #process and generate the outputs!
                      prnt.errors(summarize.errors(out,
                                                     par = par,
                                                     true = self$mean_function$parameters[par],
                                                     alpha=alpha),
                                   digits = digits)
                      
                      prnt.stats(summarize.stats(out,
                                                   par = par),
                                  digits=2)
                      
                      return(out)
                      
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
                        #UPDATE!!
                        
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
                    sim_data = function(type = "y"){
                      # if(!is.null(par)){
                      #   if(is.null(value))stop("set parameter value")
                      #   orig_par <- self$mean_function$parameters[par]
                      #   self$mean_function$parameters[par] <- value
                      #   self$mean_function$check(verbose = FALSE)
                      # }
                      re <- MASS::mvrnorm(n=1,mu=rep(0,nrow(self$covariance$D)),Sigma = self$covariance$D)
                      mu <- c(drop(as.matrix(self$mean_function$X)%*%self$mean_function$parameters)) + as.matrix(self$covariance$Z%*%re)
                      
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
                      # if(!missing(par)){
                      #   self$mean_function$parameters[par]<-orig_par
                      #   self$check(verbose = FALSE)
                      # }
                      if(type=="data.frame")y <- cbind(y,self$data$data,self$covariance$location)
                      return(y)
                      
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
                                    start,
                                    verbose=TRUE,
                                    tol = 1e-2,
                                    m=100,
                                    max.iter = 30,
                                    options = list()){
                     
                      #set options
                      if(!is(options,"list"))stop("options should be a list")
                      b_se_only <- ifelse("b_se_only"%in%names(options),options$b_se_only,FALSE)
                      use_cmdstanr <- ifelse("use_cmdstanr"%in%names(options),options$use_cmdstanr,FALSE)
                      skip_cov_optim <- ifelse("skip_cov_optim"%in%names(options),options$skip_cov_optim,FALSE)
                      skip_se <- ifelse("skip_se"%in%names(options),options$skip_se,FALSE)
                      
                      P <- ncol(self$mean_function$X)
                      R <- length(unlist(self$covariance$parameters))
                      
                      family <- self$mean_function$family[[1]]
                      
                      parInds <- list(b = 1:P,
                                      cov = (P+1):(P+R),
                                      sig = P+R+1)
                      
                      #check starting values
                      if(family%in%c("gaussian")){
                        if(missing(start)){
                          start <- c(self$mean_function$parameters,unlist(self$covariance$parameters),self$var_par)
                        }
                        if(length(start)!=(P+R+1))stop("wrong number of starting values")
                        all_pars <- 1:(P+R+1)
                      }
                      
                      if(family%in%c("binomial","poisson")){
                        if(!missing(start)){
                          if(length(start)!=(P+R)){
                            stop("wrong number of starting values")
                          }
                        } else {
                          if(verbose)message("starting values not set, setting defaults")
                          start <- c(self$mean_function$parameters,unlist(self$covariance$parameters))
                        }
                        start <- c(start,1)
                        all_pars <- 1:(P+R)
                      }
                      
                      theta <- start
                      thetanew <- rep(1,length(theta))
                      if(verbose)cat("\nStart: ",start[all_pars],"\n")
                      iter <- 0
                      niter <- m
                      Q = ncol(self$covariance$Z)
                      
                      D_func <- function(pars){
                        self$covariance$parameters <- relist(self$covariance$parameters,
                                                             pars)[[1]]
                        self$check(verbose=FALSE)
                        logdetD <- 2*sum(log(Matrix::diag(Matrix::chol(self$covariance$D))))
                        invD <- Matrix::solve(self$covariance$D)
                        Q <- dim(dsamps)[3]
                        
                        dmv <- sapply(1:dim(dsamps)[1],function(i)log_mvnd(dsamps[i,],logdetD,invD))
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
                      while(any(abs(theta-thetanew)>tol)&iter <= max.iter){
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
                          dsamps <- matrix(dsamps[,1,],ncol=Q)
                        } else {
                          capture.output(suppressWarnings(fit <- rstan::sampling(stanmodels[[gsub(".stan","",file_type$file)]],
                                                                data = data,
                                                                chains = 1,
                                                                warmup = 500,
                                                                iter = 500+m)))
                          dsamps <- rstan::extract(fit,"gamma",permuted=FALSE)
                          dsamps <- matrix(dsamps[,1,],ncol=Q)
                        }
                        
                        resids <- sapply(1:niter, function(i) y-invfunc(Matrix::drop(Xb - self$covariance$Z%*%Matrix::Matrix(dsamps[i,]))))
                        
                        
                        sigmas <- c(apply(resids,2,sd))
                        sigmahat <- mean(sigmas)
                        XtWXlist <- lapply(1:niter,function(i)Matrix::crossprod(self$mean_function$X,Matrix::diag(sigmas[i],self$n()))%*%
                                             self$mean_function$X)
                        EXtWX <- Reduce(`+`,XtWXlist)/niter
                        Wulist <- lapply(1:niter,function(i)Matrix::diag(sigmas[i],self$n())%*%resids[,i])
                        EWu <- Reduce(`+`,Wulist)/niter
                        
                        theta[parInds$b] <-  theta[parInds$b] + Matrix::drop(Matrix::solve(EXtWX)%*%Matrix::t(self$mean_function$X)%*%EWu)
                        theta[parInds$sig] <- sigmahat
                        
                        if(!skip_cov_optim){
                          opts <-  minqa::bobyqa(theta[parInds$cov],
                                                 fn = D_func,
                                                 lower = rep(1e-6,length(parInds$cov)),
                                                 control = list(rhobeg = 0.02,
                                                                rhoend = 0.02*tol))
                          
                          theta[parInds$cov] <- opts$par
                        }
                        
                        if(verbose)cat("\ntheta:",theta[all_pars])
                      }
                      
                      #need to rewrite L_func for each families likelihood!
                      L_func <- function(pars,family){
                        niter <- dim(dsamps)[1]
                        dlist <- sapply(1:niter,function(i)self$covariance$Z %*% Matrix::Matrix(dsamps[i,]))
                        zd <- Matrix::drop(Reduce(rbind,dlist))
                        mu <- rep(Matrix::drop(self$mean_function$X %*% pars[1:P]),niter) + zd
                        mu <- invfunc(mu)
                        n <- self$n()
                        # change density for the model
                        if(family=="gaussian"){
                          lf <- log(dnorm(rep(y,niter),mu,pars[P+1]))
                        } else if(family=="binomial"){
                          lf <- log(dbinom(rep(y,niter),1,mu))
                        } else if(family=="poisson"){
                          lf <- log(dpois(rep(y,niter),mu))
                        }
                        lfa <- aggregate(lf,list(rep(1:niter,each=n)),sum)
                        -mean(lfa$x)
                      }
                      
                      log_lik <- function(pars,L,family){
                        l1 <- L_func(pars[1:L],family=family)
                        l2 <- D_func(pars[(L+1):length(pars)])
                        -l1-l2
                      }
                      
                      #FULL HESSIAN WON'T WORK FOR GLM MODELS AS SIGMA NOT A PARAMETER IN THE MODEL
                      # HAVE REMOVED FROM SE EVALUATION FOR ALL MODELS FOR NOW
                      # BUT BETTER TO REPORT ONLY FOR GAUSSIAN MODELS
                      
                      if(verbose)cat("\n\nCalculating standard errors...\n")
                      
                      if(family%in%c("gaussian")){
                        mf_pars <- theta[c(parInds$b,parInds$sig)]
                        mf_pars_names <- c(paste0("b",1:P),"sigma")
                      } else {
                        mf_pars <- theta[c(parInds$b)]
                        mf_pars_names <- c(paste0("b",1:P))
                      }
                      
                      cov_pars_names <- paste0("cov",1:R)
                      
                      if(!skip_se){
                        if(b_se_only){
                          hess <- tryCatch(pracma::hessian(function(x,family)-L_func(x,family),
                                                           mf_pars,
                                                           family=family),error=function(e)NULL)
                          if(!is.null(hess)){
                            SE <- tryCatch(sqrt(Matrix::diag(Matrix::solve(-hess))),error=function(e)rep(NA,length(mf_pars)))
                          } else {
                            SE <- rep(NA,length(mf_pars))
                          }
                          
                          res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                            est = c(mf_pars,theta[parInds$cov]),
                                            SE=c(SE,rep(NA,length(parInds$cov) )))
                        } else {
                          hess <- tryCatch(pracma::hessian(log_lik,
                                                           c(mf_pars,theta[parInds$cov]),
                                                           L = length(mf_pars),
                                                           family = family),error=function(e)NULL)
                          if(!is.null(hess)){
                            SE <- tryCatch(sqrt(Matrix::diag(Matrix::solve(-hess))),error=function(e)rep(NA,length(mf_pars)+length(cov_pars_names)))
                          } else {
                            SE <- rep(NA,length(mf_pars)+length(cov_pars_names))
                          }
                          res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                            est = c(mf_pars,theta[parInds$cov]),
                                            SE=SE)
                        }
                        
                        if(any(is.na(res$SE[1:P]))){
                          warning("Hessian was not positive definite, using approximation")
                          self$check(verbose=FALSE)
                          res$SE[1:P] <- sqrt(Matrix::diag(Matrix::solve(private$information_matrix())))
                        }
                          
                      } else {
                        res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                          est = c(mf_pars,theta[parInds$cov]),
                                          SE=NA)
                      }
                      
                      if(iter >= max.iter|any(abs(theta-thetanew)>tol))warning("not converged")
                      
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
                    saved_sim_data = NULL,
                    gen_sim_data = function(...){
                      ysim <- self$sim_data()
                      res <- do.call(self$MCNR,list(y=ysim,verbose=FALSE,...))
                      return(res[,c('est','SE')])
                    },
                    hash = NULL,
                    hash_do = function(){
                      digest::digest(c(self$covariance$.__enclos_env__$private$hash,
                                       self$mean_function$.__enclos_env__$private$hash))
                    },
                    information_matrix = function(){
                      Matrix::crossprod(self$mean_function$X,solve(self$Sigma))%*%self$mean_function$X
                    }
                  ))


