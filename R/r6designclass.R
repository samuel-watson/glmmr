#' R6 Class representing a GLMM Design
#' 
#' For the generalised linear mixed model 
#' 
#' \deqn{Y \sim F(\mu,\sigma)}
#' \deqn{\mu = h^-1(X\beta + Z\gamma)}
#' \deqn{\gamma \sim MVN(0,D)}
#' 
#' where h is the link function. A Design in comprised of a \link[glmmr]{MeanFunction} object, which defines the family F, 
#' link function h, and fixed effects design matrix X, and a \link[glmmr]{Covariance} object, which defines Z and D. The class provides
#' methods for analysis and simulation with these models.
#' @importFrom Matrix Matrix
#' @export 
Design <- R6::R6Class("Design",
                  public = list(
                    #' @field covariance A \link[glmmr]{Covariance} object defining the random effects covariance.
                    covariance = NULL,
                    #' @field mean_function A \link[glmmr]{MeanFunction} object, defining the mean function for the model, including the data and covariate design matrix X.
                    mean_function = NULL,
                    #' @field exp_condition A vector indicting the unique experimental conditions for each observation, see Details.
                    exp_condition = NULL,
                    #' @field Sigma The overall covariance matrix for the observations. Calculated and updated automatically as \eqn{W^{-1} + ZDZ^T} where W is an n x n 
                    #' diagonal matrix with elements on the diagonal equal to the GLM iterated weights. See Details.
                    Sigma = NULL,
                    #' @field var_par Scale parameter required for some distributions (Gaussian, Gamma, Beta).
                    var_par = NULL,
                    #' @description 
                    #' Return predicted values based on the currently stored parameter values in `mean_function`
                    #' @param type One of either "`link`" for values on the scale of the link function, or "`response`" 
                    #' for values on the scale of the response
                    #' @return A \link[Matrix]{Matrix} class object containing the predicted values
                    fitted = function(type="link"){
                      Xb <- Matrix::drop(self$mean_function$X %*% self$mean_function$parameters)
                      if(type=="response"){
                        Xb <- self$mean_function$family$linkinv(Xb)
                      }
                      return(Xb)
                    },
                    #' @description 
                    #' Create a new Design object
                    #' @param covariance Either a \link[glmmr]{Covariance} object, or an equivalent list of arguments
                    #' that can be passed to `Covariance` to create a new object.
                    #' @param mean.function Either a \link[glmmr]{MeanFunction} object, or an equivalent list of arguments
                    #' that can be passed to `MeanFunction` to create a new object.
                    #' @param var_par Scale parameter required for some distributions, including Gaussian. Default is NULL.
                    #' @param verbose Logical indicating whether to provide detailed output
                    #' @return A new Design class object
                    initialize = function(covariance,
                                          mean.function,
                                          var_par = NULL,
                                          verbose=TRUE,
                                          skip.sigma = FALSE){
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

                      if(!skip.sigma)self$generate()
                      private$hash <- private$hash_do()
                    },
                    #' @description 
                    #' Print method for `Design` class
                    #' @details 
                    #' Calls the respective print methods of the linked covariance and mean function objects.
                    #' @param ... ignored
                    print = function(){
                      cat("\n----------------------------------------\n")
                      print(self$mean_function)
                      cat("\n----------------------------------------\n")
                      print(self$covariance)
                      cat("\n----------------------------------------\n")
                    },
                    #' @description 
                    #' Returns the number of observations in the model
                    #' @details 
                    #' The matrices X and Z both have n rows, where n is the number of observations in the model/design.
                    #' @param ... ignored
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
                    analysis = function(type,
                                        iter,
                                        par,
                                        alpha = 0.05,
                                        sim_design,
                                        parallel,
                                        verbose = TRUE,
                                        digits = 2,
                                        ...){
                      
                      if(!missing(sim_design)){
                        f1 <- sim_design$sim_data
                      } else {
                        f1 <- self$sim_data
                      }
                      if(type=="sim_data"&is.null(private$saved_sim_data))stop("no simulation data saved in object")
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
                                                   function(i){
                                                     ysim <- f1()
                                                     private$gen_sim_data(par=par,
                                                                          ysim = ysim,...)},
                                                   cl=cl)
                          parallel::stopCluster(cl)
                        } else {
                          out <- pbapply::pblapply(1:iter,
                                                   function(i){
                                                     ysim <- f1()
                                                     private$gen_sim_data(par=par,
                                                                          ysim = ysim,...)})
                        }
                        
                        res <- list(
                          coefficients = lapply(out,function(i)i[[1]]$coefficients),
                          #dfbeta = lapply(out,function(i)i[[2]]),
                          sim_method = "full.sim",
                          mcml_method = out[[1]][[1]]$method,
                          convergence = unlist(lapply(out,function(i)i[[1]]$converged)),
                          m = out[[1]][[1]]$m,
                          tol =out[[1]][[1]]$tol,
                          nsim = iter,
                          alpha = alpha,
                          b_parameters = self$mean_function$parameters,
                          cov_parameters = self$covariance$parameters,
                          mean_formula = self$mean_function$formula,
                          cov_formula = self$covariance$formula,
                          family = self$mean_function$family,
                          n = self$n(),
                          par = par
                        )
                        
                        class(res) <- "glmmr.sim"
                        
                        if(verbose)message("saving simulation data")
                        private$saved_sim_data <- res
                      }
                      if(type == "sim_approx"){
                        if(parallel){
                          cl <- parallel::makeCluster(parallel::detectCores()-1)
                          parallel::clusterEvalQ(cl,library(Matrix))
                          #change when package built!
                          parallel::clusterEvalQ(cl,devtools::load_all())
                          # out <- parallel::parLapply(cl,
                          #                            1:10,
                          #                            function(i)self$gen_sim_data(m=m))
                          out <- pbapply::pblapply(1:iter,
                                                   function(i){
                                                     ysim <- f1()
                                                     private$gen_sim_data_approx(par=par,
                                                                                 ysim=ysim)},
                                                   cl=cl)
                          parallel::stopCluster(cl)
                        } else {
                          out <- pbapply::pblapply(1:iter,
                                                   function(i){
                                                     ysim <- f1()
                                                     private$gen_sim_data_approx(par=par,
                                                                                 ysim=ysim)})
                        }
                        
                        res <- list(
                          coefficients = lapply(out,function(i)i[[1]]),
                          #dfbeta = lapply(out,function(i)i[[2]]),
                          sim_method = "approx.sim",
                          mcml_method = NA,
                          convergence = NA,
                          m = NA,
                          tol = NA,
                          nsim = iter,
                          alpha = alpha,
                          b_parameters = self$mean_function$parameters,
                          cov_parameters = self$covariance$parameters,
                          mean_formula = self$mean_function$formula,
                          cov_formula = self$covariance$formula,
                          family = self$mean_function$family,
                          n = self$n(),
                          par = par
                        )
                        
                        class(res) <- "glmmr.sim"
                        
                        if(verbose)message("saving simulation data")
                        private$saved_sim_data <- res
                      }
                      if(type=="sim_data")res <- private$saved_sim_data
                      
                      invisible(res)
                      
                    },
                    power = function(par,
                                     value,
                                     alpha=0.05,
                                     method=NULL,
                                     iter=10,
                                     skip.check = FALSE,
                                     parallel=TRUE){
                      if(!skip.check)self$check(verbose=FALSE)
                      if(missing(par)|missing(value))stop("parameter missing")
                      old_par <- self$mean_function$parameters[[par]]
                      self$mean_function$parameters[[par]] <- value
                      self$check(verbose=FALSE)
                      M <- private$information_matrix()
                      v0 <- solve(M)[par,par]
                      pwr <- pnorm(value/(sqrt(v0)) - qnorm(1-alpha/2))
                      self$mean_function$parameters[[par]] <- old_par
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
                        ggplot2::ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          ggplot2::geom_count(aes(color=self$mean_function$data[,treat]))+
                          ggplot2::theme_bw()+
                          ggplot2::theme(panel.grid=ggplot2::element_blank())+
                          ggplot2::scale_color_viridis_c(name=treat)+
                          ggplot2::scale_size_area()
                      } else {
                        ggplot2::ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          ggplot2::geom_count(aes(color=self$mean_function$data[,treat]))+
                          ggplot2::facet_wrap(~.data[[z]])+
                          ggplot2::theme_bw()+
                          ggplot2::theme(panel.grid=ggplot2::element_blank())+
                          ggplot2::scale_color_viridis_c(name=treat)+
                          ggplot2::scale_size_area()
                      }},
                    sim_data = function(type = "y"){
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
                    MCML = function(y,
                                    start,
                                    se.method = "lik",
                                    permutation.par,
                                    verbose=TRUE,
                                    tol = 1e-2,
                                    m=100,
                                    max.iter = 30,
                                    options = list()){
                     
                      # checks
                      if(!se.method%in%c("perm","lik"))stop("se.method should be 'perm' or 'lik'")
                      if(se.method=="perm" & missing(permutation.par))stop("if using permutational based
inference, set permuation.par")
                      if(se.method=="perm" & is.null(self$mean_function$randomise))stop("random allocations
are created using the function in self$mean_function$randomise, but this has not been set. Please see help(MeanFunction)
for more details")
                      #set options
                      if(!is(options,"list"))stop("options should be a list")
                      b_se_only <- ifelse("b_se_only"%in%names(options),options$b_se_only,FALSE)
                      use_cmdstanr <- ifelse("use_cmdstanr"%in%names(options),options$use_cmdstanr,FALSE)
                      skip_cov_optim <- ifelse("skip_cov_optim"%in%names(options),options$skip_cov_optim,FALSE)
                      #skip_se <- ifelse("skip_se"%in%names(options),options$skip_se,FALSE)
                      method <- ifelse("method"%in%names(options),options$method,"mcnr")
                      sim_lik_step <- ifelse("sim_lik_step"%in%names(options),options$sim_lik_step,FALSE)
                      no_warnings <- ifelse("no_warnings"%in%names(options),options$no_warnings,FALSE)
                      perm_type <- ifelse("perm_type"%in%names(options),options$perm_type,"cov")
                      perm_iter <- ifelse("perm_iter"%in%names(options),options$perm_iter,100)
                      perm_parallel <- ifelse("perm_iter"%in%names(options),options$perm_iter,TRUE)
                      warmup_iter <- ifelse("warmup_iter"%in%names(options),options$warmup_iter,500)
                      perm_ci_steps <- ifelse("warmup_iter"%in%names(options),options$perm_ci_steps,1000)
                      #theta_in <- ifelse("theta_in"%in%names(options),options$theta_in,NULL)
                      
                      P <- ncol(self$mean_function$X)
                      R <- length(unlist(self$covariance$parameters))
                      
                      family <- self$mean_function$family[[1]]
                      
                      parInds <- list(b = 1:P,
                                      cov = (P+1):(P+R),
                                      sig = P+R+1)
                      
                      if(family%in%c("gaussian")){
                        mf_parInd <- c(parInds$b,parInds$sig)
                      } else {
                        mf_parInd <- c(parInds$b)
                      }
                      
                      orig_par_b <- self$mean_function$parameters
                      orig_par_cov <- self$covariance$parameters
                      
                      #check starting values
                      if(family%in%c("gaussian")){
                        if(missing(start)){
                          if(verbose)message("starting values not set, setting defaults")
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
                      
                      if(verbose)message(paste0("using method: ",method))
                      if(verbose)cat("\nStart: ",start[all_pars],"\n")
                      
                      iter <- 0
                      niter <- m
                      Q = ncol(self$covariance$Z)
                      
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
                      
                      ## ALGORITHMS
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
                                                           iter_warmup = warmup_iter,
                                                           iter_sampling = m,
                                                           refresh = 0),
                                         file=tempfile())
                          dsamps <- fit$draws("gamma")
                          dsamps <- matrix(dsamps[,1,],ncol=Q)
                        } else {
                          capture.output(suppressWarnings(fit <- rstan::sampling(stanmodels[[gsub(".stan","",file_type$file)]],
                                                                data = data,
                                                                chains = 1,
                                                                warmup = warmup_iter,
                                                                iter = warmup_iter+m)))
                          dsamps <- rstan::extract(fit,"gamma",permuted=FALSE)
                          dsamps <- matrix(dsamps[,1,],ncol=Q)
                        }
                        
                        dsamps <<- dsamps
                        # BETA PARAMETERS STEP
                        if(method == "mcnr"){
                          beta_step <- mcnr_step(y,
                                                 as.matrix(self$mean_function$X),
                                                 as.matrix(self$covariance$Z),
                                                 theta[parInds$b],
                                                 dsamps,
                                                 self$mean_function$family[[2]])
                          
                          theta[parInds$b] <-  theta[parInds$b] + beta_step$beta_step
                          theta[parInds$sig] <- beta_step$sigmahat
                          
                          
                        } else if(method == "mcem"){
                          theta[mf_parInd] <- drop(l_lik_optim(as.matrix(self$covariance$Z),
                                                          as.matrix(self$mean_function$X),
                                                          y,
                                                          dsamps,
                                                          family=self$mean_function$family[[1]],
                                                          link=self$mean_function$family[[2]],
                                                          start = theta[mf_parInd],
                                                          lower = rep(-Inf,length(mf_parInd)),
                                                          upper = rep(Inf,length(mf_parInd)),
                                                          trace= 0))
                          
                        }
                        
                        
                        # COVARIANCE PARAMETERS STEP
                        if(!skip_cov_optim){
                          theta[parInds$cov] <- drop(d_lik_optim(self$covariance$.__enclos_env__$private$Funclist,
                                                      self$covariance$.__enclos_env__$private$Distlist,
                                                      dsamps,
                                                      start = c(theta[parInds$cov]),
                                                      lower= rep(0,length(parInds$cov)),
                                                      upper= rep(Inf,length(parInds$cov)),
                                                      trace=0))
                          
                          
                        }
                        
                        if(verbose)cat("\ntheta:",theta[all_pars])
                      }
                      
                      not_conv <- iter >= max.iter|any(abs(theta-thetanew)>tol)
                      if(not_conv&!no_warnings)warning(paste0("algorithm not converged. Max. difference between iterations :",max(abs(theta-thetanew)),". Suggest 
                                                 increasing m, or trying a different algorithm."))
                      
                      if(sim_lik_step){
                        if(verbose)cat("\n\n")
                        if(verbose)message("Optimising simulated likelihood")
                        theta[all_pars] <- f_lik_optim(self$covariance$.__enclos_env__$private$Funclist,
                                           self$covariance$.__enclos_env__$private$Distlist,
                                           as.matrix(self$covariance$Z),
                                           as.matrix(self$mean_function$X),
                                           y,
                                           dsamps,
                                           theta[parInds$cov],
                                           family=self$mean_function$family[[1]],
                                           link=self$mean_function$family[[2]],
                                           start = theta[all_pars],
                                           lower = c(rep(-Inf,P),rep(1e-5,length(all_pars)-P)),
                                           importance = TRUE)
                        
                      }
                      
                      if(verbose)cat("\n\nCalculating standard errors...")
                      
                      if(family%in%c("gaussian")){
                        mf_pars <- theta[c(parInds$b,parInds$sig)]
                        mf_pars_names <- c(colnames(self$mean_function$X),"sigma")
                      } else {
                        mf_pars <- theta[c(parInds$b)]
                        mf_pars_names <- colnames(self$mean_function$X)
                      }
                      
                      cov_pars_names <- paste0("cov",1:R)
                      permutation = FALSE
                      if(se.method=="lik"){
                        if(verbose)cat("using Hessian method\n")
                          hess <- tryCatch(f_lik_optim(self$covariance$.__enclos_env__$private$Funclist,
                                                       self$covariance$.__enclos_env__$private$Distlist,
                                                       as.matrix(self$covariance$Z),
                                                       as.matrix(self$mean_function$X),
                                                       y,
                                                       dsamps,
                                                       theta[parInds$cov],
                                                       family=self$mean_function$family[[1]],
                                                       link=self$mean_function$family[[2]],
                                                       start = theta[all_pars],
                                                       lower = c(rep(-Inf,P),rep(1e-5,length(all_pars)-P)),
                                                       importance = FALSE),error=function(e)NULL)
                          
                          hessused <- TRUE
                          if(!is.null(hess)){
                            SE <- tryCatch(sqrt(Matrix::diag(Matrix::solve(hess))),error=function(e)rep(NA,length(mf_pars)+length(cov_pars_names)))
                          } else {
                            SE <- rep(NA,length(mf_pars)+length(cov_pars_names))
                          }
                          res <- data.frame(par = c(mf_pars_names,cov_pars_names,paste0("d",1:Q)),
                                            est = c(mf_pars,theta[parInds$cov],colMeans(dsamps)),
                                            SE=c(SE,apply(dsamps,2,sd)))
                        
                        
                        if(any(is.na(res$SE[1:P]))){
                          if(!no_warnings)warning("Hessian was not positive definite, using approximation")
                          hessused <- FALSE
                          self$check(verbose=FALSE)
                          res$SE[1:P] <- sqrt(Matrix::diag(Matrix::solve(private$information_matrix())))
                        }
                          
                        res$lower <- res$est - qnorm(1-0.05/2)*res$SE
                        res$upper <- res$est + qnorm(1-0.05/2)*res$SE
                          
                      } else if(se.method=="perm") {
                        if(verbose)cat("using permutational method\n")
                        permutation = TRUE
                        #get null model
                        # use parameters from fit above rather than null marginal model
                        perm_out <- self$perumtation_test(permutation.par,
                                                          start = theta[parInds$b][permutation.par],
                                                          nsteps = perm_ci_steps,
                                                          type = perm_type,
                                                          verbose= verbose)
                        tval <- qnorm(1-perm_out$p/2)
                        par <- theta[parInds$b][permutation.par]
                        se <- abs(par/tval)
                        se1 <- rep(NA,length(mf_pars))
                        se1[permutation.par] <- se
                        se2 <- rep(NA,length(parInds$cov))
                        ci1l <- ci1u <- rep(NA,length(mf_pars))
                        ci2l <- ci2u <- rep(NA,length(parInds$cov))
                        ci1l[permutation.par] <- perm_out$lower
                        ci1u[permutation.par] <- perm_out$upper
                        
                        res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                          est = c(mf_pars,theta[parInds$cov]),
                                          SE=c(se1,se2),
                                          lower=c(ci1l,ci2l),
                                          upper=c(ci1u,ci2u))
                      hessused <- FALSE
                      } else {
                        res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                          est = c(mf_pars,theta[parInds$cov]),
                                          SE=NA,
                                          lower = NA,
                                          upper =NA)
                      }
                      
                     colnames(dsamps) <- Reduce(c,rev(cov1$.__enclos_env__$private$flistlabs))
                     out <- list(coefficients = res,
                                 converged = !not_conv,
                                 method = method,
                                 hessian = hessused,
                                 permutation = permutation,
                                 m = m,
                                 tol = tol,
                                 sim_lik = sim_lik_step,
                                 re.samps = dsamps)
                     
                     class(out) <- "mcml"
                      
                      self$mean_function$parameters <- orig_par_b 
                      self$covariance$parameters <- orig_par_cov
                      #self$check(verbose=FALSE)
                      
                      return(out)
                    },
                    # dfbeta = function(y,
                    #                   par,
                    #                   alpha=0.05,
                    #                   b_hat ,
                    #                   verbose = TRUE){
                    #   
                    #   orig_par_b <- self$mean_function$parameters
                    #   orig_par_cov <- self$covariance$parameters
                    #   
                    #   iter <- 0
                    #   sig <- NULL
                    #   sign <- NULL
                    #   sigsign <- NULL
                    #   if(!missing(b_hat)){
                    #     ests <- b_hat
                    #     self$mean_function$parameters <- b_hat
                    #   } else {
                    #     if(verbose)message("parameter estimates not provided, using values stored in mean function")
                    #     ests <- self$mean_function$parameters
                    #   }
                    #   self$check(verbose = verbose)
                    #   newest <- ests[par]
                    #   n <- self$n()
                    #   nid <- 1:n
                    #   
                    #   
                    #   #fix the indexes to remove
                    #   
                    #   while(is.null(sig)|is.null(sign)|is.null(sigsign)){
                    #     invS <- Matrix::solve(self$Sigma[nid,nid])
                    #     invSX <- invS %*% self$mean_function$X[nid,]
                    #     invM <- Matrix::solve(Matrix::crossprod(self$mean_function$X[nid,],invSX))
                    #     B <-  invM %*% Matrix::t(invSX)
                    #     Q <- invS - invSX %*% B
                    #     dQ <- 1/Matrix::diag(Q)
                    #     e <- dQ*Matrix::t(Q)%*%y[nid]
                    #     Be <- B%*%Matrix::diag(Matrix::drop(e))
                    #     SE <- sqrt(Matrix::diag(invM))
                    #     tstat <- newest/SE
                    #     pval <- 2*(1-pnorm(abs(tstat)))
                    #     maxid <- order(Be[par,],decreasing = sign(ests[par])==1)[1]
                    #     newest <- newest - Be[par,maxid]
                    #     self$mean_function$parameters <- self$mean_function$parameters - Be[,maxid]
                    #     self$check(verbose = FALSE)
                    #     
                    #     nid <- nid[-maxid]
                    #     iter <- iter+1
                    #     # check significance
                    #     if(pval[par] >= alpha & is.null(sig)){
                    #       sig <- iter
                    #       attr(sig,"id") <- which(!c(1:n)%in%nid)
                    #     }
                    #     if(sign(ests[par]) != sign( newest ) & is.null(sign)){
                    #       sign <- iter
                    #       attr(sign,"id") <- which(!c(1:n)%in%nid)
                    #     }
                    #     if(sign(ests[par]) != sign( newest ) & pval[par] < alpha &  is.null(sigsign)){
                    #       sigsign <- iter
                    #       attr(sigsign,"id") <- which(!c(1:n)%in%nid)
                    #     }
                    #     
                    #     
                    #     if(verbose)cat("\rIteration: ",iter)
                    #     
                    #   }
                    #   
                    #   self$mean_function$parameters <- orig_par_b 
                    #   self$covariance$parameters <- orig_par_cov
                    #   self$check(verbose=FALSE)
                    #   
                    #   
                    #   
                    #   return(list(sig,sign,sigsign))
                    # },
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
                    },
                    permutation_test = function(y,
                                                permutation.par,
                                           start,
                                           iter = 1000,
                                           nsteps=1000,
                                           type="cov",
                                           parallel = TRUE,
                                           verbose=TRUE){
                      if(is.null(self$mean_function$randomise))stop("random allocations
are created using the function in self$mean_function$randomise, but this has not been set. Please see help(MeanFunction)
for more details")
                      Xnull <- as.matrix(self$mean_function$X)
                      Xnull <- Xnull[,-permutation.par]
                      null_fit <- stats::glm.fit(Xnull,y,family=self$mean_function$family)
                      xb <- null_fit$linear.predictors
                      
                      tr <- self$mean_function$X[,permutation.par]
                      if(any(!tr%in%c(0,1)))stop("permuational inference only available for dichotomous treatments")
                      tr[tr==0] <- -1
                      
                      if(verbose&type=="cov")message("using covariance weighted statistic, to change permutation statistic set option perm_type, see details in help(Design)")
                      if(verbose&type=="unw")message("using unweighted statistic, to change permutation statistic set option perm_type, see details in help(Design)")
                      w.opt <- type=="cov"
                      invS <- ifelse(w.opt,Matrix::solve(self$Sigma),1)
                      qstat <- private$qscore(y,tr,xb,permutation.par,invS,w.opt)
                      
                      if(verbose)cat("Starting permutations\n")
                      if(parallel){
                        cl <- parallel::makeCluster(parallel::detectCores()-1)
                        parallel::clusterEvalQ(cl,library(Matrix))
                        #change when package built!
                        parallel::clusterEvalQ(cl,devtools::load_all())
                        qtest <- pbapply::pbsapply(1:iter,function(i){
                          new_tr <- self$mean_function$randomise()
                          new_tr[new_tr==0] <- -1
                          private$qscore(y,new_tr,xb,permutation.par,invS,w.opt)
                        }, cl = cl)
                        parallel::stopCluster(cl)
                      } else {
                        qtest <- pbapply::pbsapply(1:iter,function(i){
                          new_tr <- self$mean_function$randomise()
                          new_tr[new_tr==0] <- -1
                          private$qscore(y,new_tr,xb,permutation.par,invS,w.opt)
                        })
                      }
                      
                      #permutation confidence intervals
                      if(verbose)cat("Starting permutational confidence intervals\n")
                      pval <- length(qtest[qtest>qstat])/iter
                      #print(pval)
                      if(pval==0)pval <- 0.5/iter
                      tval <- qnorm(1-pval/2)
                      par <- start#theta[parInds$b][permutation.par]
                      se <- abs(par/tval)
                      if(verbose)cat("Lower\n")
                      lower <- private$confint_search(y,
                                                      permutation.par,
                                                      start = par - 2*se,
                                                      b = start,
                                                      w.opt = w.opt,
                                                      nsteps = nsteps)
                      if(verbose)cat("\nUpper\n")
                      upper <- private$confint_search(y,
                                                      permutation.par,
                                                      start = par + 2*se,
                                                      b = start,
                                                      w.opt = w.opt,
                                                      nsteps = nsteps)
                      return(list(p=pval,lower=lower,upper=upper))
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
                    gen_sim_data = function(par,
                                            ysim,
                                            ...){
                      
                      #ysim <- self$sim_data()
                      # choose mcem for glmm and mcnr for lmm
                      res <- do.call(self$MCML,list(y=ysim,
                                                    verbose=TRUE,
                                                    options= list(method="mcem",
                                                                  no_warnings=TRUE),...))
                      dfb <- NULL#self$dfbeta(y=ysim,
                                  #       par = par,
                                   #      b_hat = res$coefficients[grepl("b",res$coefficients$par),'est'])
                      
                      return(list(res,dfb))
                    },
                    gen_sim_data_approx = function(par,
                                                   ysim,
                                            ...){
                      
                      #ysim <- self$sim_data()
                      # choose mcem for glmm and mcnr for lmm
                      invM <- Matrix::solve(private$information_matrix())
                      b <- Matrix::drop(invM %*% Matrix::crossprod(self$mean_function$X,Matrix::solve(self$Sigma))%*%ysim)
                      res <- data.frame(par = paste0("b",1:length(b)),est=b,SE= sqrt(Matrix::drop(Matrix::diag(invM))))
                      dfb <- NULL#self$dfbeta(y=ysim,
                                #         par = par,
                                 #        b_hat = res[grepl("b",res$par),'est'])
                      
                      return(list(res,dfb))
                    },
                    hash = NULL,
                    hash_do = function(){
                      digest::digest(c(self$covariance$.__enclos_env__$private$hash,
                                       self$mean_function$.__enclos_env__$private$hash))
                    },
                    information_matrix = function(){
                      Matrix::crossprod(self$mean_function$X,solve(self$Sigma))%*%self$mean_function$X
                    },
                    qscore = function(y,
                                      tr,
                                      xb,
                                      permutation.par,
                                      invS,
                                      weight=TRUE){
                      
                      #xb <- self$mean_function$X[,-permutation.par] %*% b + offset
                      ypred <- self$mean_function$family$linkinv(Matrix::drop(xb))
                      resids <- Matrix::Matrix(y-ypred)
                      
                      if(weight){
                        tr_mat <- Matrix::diag(tr)
                        g <- Matrix::t(private$get_G(xb))
                        q <- (g%*%invS)%*%(tr_mat%*%resids)
                      } else {
                        tr_mat <- Matrix::Matrix(tr,nrow=1)
                        q <- (tr_mat%*%resids)
                      }
                      return(abs(Matrix::drop(q)))
                    },
                    get_G = function(x){
                      family = self$mean_function$family
                      
                      if(family[[2]] == "identity"){
                        dx <- rep(1,length(x))
                      } else if(family[[2]] == "log"){
                        dx <- exp(x)
                      } else if(family[[2]] == "logit"){
                        dx <- (exp(x)/(1+exp(x)))*(1-exp(x)/(1+exp(x)))
                      } else if(family[[2]] == "probit"){
                        dx <- -1/dnorm(x)
                      }
                      return(dx)
                    },
                    confint_search = function(y,
                                              permutation.par,
                                              start,
                                              b,
                                              nsteps = 1000,
                                              w.opt = TRUE,
                                              alpha=0.05,
                                              verbose=TRUE){
                      bound <- start
                      Xnull <- as.matrix(self$mean_function$X)
                      Xnull <- Xnull[,-permutation.par]
                      tr <- as.matrix(self$mean_function$X)[,permutation.par]
                      dtr <- tr
                      dtr[dtr==0] <- -1
                      k <- 2/(pnorm(1-alpha)*((2*pi)^-0.5)*exp((-pnorm(1-alpha)^2)/2))
                      invS <- ifelse(w.opt,Matrix::solve(self$Sigma),1)
                      
                      for(i in 1:nsteps){
                        null_fit <- stats::glm.fit(Xnull,y,family=self$mean_function$family,offset = tr*bound)
                        #pars <- null_fit$coefficients
                        xb <- null_fit$linear.predictors
                        qstat <- private$qscore(y,dtr,xb,permutation.par,invS,w.opt)
                        
                        new_tr <- self$mean_function$randomise()
                        new_tr[new_tr==0] <- -1
                        qtest <- private$qscore(y,new_tr,xb,permutation.par,invS,w.opt)
                        rjct <- qstat > qtest
                        #cat("\r bound: ",bound," qscore: ",qstat," qtest: ",qtest," iter: ",i)
                        
                        step <- k*(b - bound)
                        if(rjct){
                          bound <- bound + step*alpha/i
                        } else {
                          bound <- bound - step*(1-alpha)/i
                        }
                        
                        if(verbose)cat("\r",progress_bar(i,nsteps))
                      }
                      
                      return(bound)
                      
                    }
                  ))


