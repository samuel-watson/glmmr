#' A GLMM Design 
#' 
#' An R6 class representing a GLMM and study design
#' @details
#' For the generalised linear mixed model 
#' 
#' \deqn{Y \sim F(\mu,\sigma)}
#' \deqn{\mu = h^-1(X\beta + Z\gamma)}
#' \deqn{\gamma \sim MVN(0,D)}
#' 
#' where h is the link function. A Design in comprised of a \link[glmmr]{MeanFunction} object, which defines the family F, 
#' link function h, and fixed effects design matrix X, and a \link[glmmr]{Covariance} object, which defines Z and D. The class provides
#' methods for analysis and simulation with these models.
#' 
#' This class provides methods for: data simulation (`sim_data()` and `fitted()`), model fitting using Markov Chain 
#' Monte Carlo Maximum Likelihood (MCML) methods (`MCML()`), design analysis via simulation including power (`analysis()`),
#' deletion diagnostics (`dfbeta()`), and permutation tests including p-values and confidence intervals (`permutation()`).
#' 
#' The class by default calculates the covariance matrix of the observations as:
#' 
#' \deqn{\Sigma = W^{-1} + ZDZ^T}
#' 
#' where _W_ is a diagonal matrix with the WLS iterated weights for each observation equal
#' to, for individual _i_ \eqn{\phi a_i v(\mu_i)[h'(\mu_i)]^2} (see Table 2.1 in McCullagh 
#' and Nelder (1989) <ISBN:9780412317606>). For very large designs, this can be disabled as
#' the memory requirements can be prohibitive.
#' @references 
#' Braun and Feng
#' McCullagh
#' Stan
#' McCullagh and Nelder
#' Approx GLMMs paper
#' Watson confidence interval
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
                    #' @seealso \link[glmmr]{nelder}, \link[glmmr]{MeanFunction}, \link[glmmr]{Covariance}
                    #' @examples 
                    #' #create a data frame describing a cross-sectional parallel cluster
                    #' #randomised trial
                    #' df <- nelder(~(cl(10)*t(5)) > ind(10))
                    #' df$int <- 0
                    #' df[df$cl > 5, 'int'] <- 1
                    #' 
                    #' mf1 <- MeanFunction$new(
                    #'   formula = ~ factor(t) + int - 1,
                    #'   data=df,
                    #'   parameters = c(rep(0,5),0.6),
                    #'   family = gaussian()
                    #' )
                    #' cov1 <- Covariance$new(
                    #'   data = df,
                    #'   formula = ~ (1|gr(cl)) + (1|gr(cl*t)),
                    #'   parameters = c(0.25,0.1)
                    #' )
                    #' des <- Design$new(
                    #'   covariance = cov1,
                    #'   mean.function = mf1,
                    #'   var_par = 1
                    #' )
                    #' 
                    #' #alternatively we can pass the data directly to Design
                    #' #here we will specify a cohort study
                    #' df <- nelder(~ind(20) > t(6))
                    #' df$int <- 0
                    #' df[df$t > 3, 'int'] <- 1
                    #' 
                    #' des <- Design$new(
                    #' covariance = list(
                    #'   data=df,
                    #'   formula = ~ (1|pexp(t)*gr(ind)),
                    #'   parameters = c(0.8,1)),
                    #' mean.function = list(
                    #'   formula = ~int + factor(t),
                    #'   data=df,
                    #'   parameters = rep(0,7),
                    #'   family = poisson()))
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

                      if(!skip.sigma)private$generate()
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
                    n = function(...){
                      self$mean_function$n()
                    },
                    #' @description
                    #' Returns the number of clusters at each level
                    #' @details
                    #' Returns a data frame describing the number of independent clusters or groups at each level in the design. For example,
                    #' if there were cluster-periods nested in clusters, then the top level would be clusters, and the second level would be 
                    #' cluster periods.
                    #' @param ... ignored
                    #' @return A data frame with the level, number of clusters, and variables describing each level.
                    #' @examples 
                    #' ...
                    n_cluster = function(){
                      flist <- rev(self$covariance$.__enclos_env__$private$flistvars)
                      gr_var <- unlist(lapply(flist,function(x)"gr"%in%x$funs))
                      gr_count <- unlist(rev(self$covariance$.__enclos_env__$private$flistcount))
                      gr_cov_var <- lapply(flist,function(x)x$rhs)
                      if(any(gr_var)){
                        dfncl <- data.frame(Level = 1:sum(gr_var),"N.clusters"=sort(gr_count[gr_var]),"Variables"=unlist(lapply(gr_cov_var,paste0,collapse=" "))[gr_var][order(gr_count[gr_var])])
                      } else {
                        dfncl <- data.frame(Level = 1,"N.clusters"=1,"Variables"=paste0(unlist(gr_cov_var)[!duplicated(unlist(gr_cov_var))],collapse=" "))
                      }
                      return(dfncl)
                      
                    },
                    #' @description 
                    #' Run a design analysis of the model via simulation
                    #' @details 
                    #' The analysis function conducts a detailed design analysis using the analysis
                    #' model specified by the object. Data are simulated either using the same
                    #' data generating process, or using a different Design object specified by 
                    #' the user to allow for model misspecification. On each iteration the model
                    #' is estimated with the simulated data _y_ using either the `MCML` function or
                    #' approximate parameters and standard errors using generalised least squares. 
                    #' MCML is an exact maximum likelihood algorithm, and can be slow, so results
                    #' from previous simulations are saved in the design object and can be recalled
                    #' later. Deletion diagnostics are also calculated to calculate influential parts of 
                    #' the design, although these are typically not useful for balanced 
                    #' experimental designs.
                    #' 
                    #' The function returns an `glmmr.sim` object, which estimates and summarises: 
                    #' 
                    #' **Model fitting and simulation diagnostics** 
                    #' Convergence failure percent and coverage. Maximum likelihood estimators
                    #' for GLMMs can fail to converge or reach the MLE. GLMMs can also be 
                    #' susceptible to small sample biases where the relevant sample size is
                    #' at the level of clustering and correlation. 
                    #' 
                    #' **Error rates** 
                    #' Type 2 (power), Type S (significance), and Type M (magnitude) errors are
                    #' reported.
                    #' 
                    #' **p-value and confidence interval width distributions**
                    #' 
                    #' **Deletion diagnostics**
                    #' For unbalanced designs and under model misspecifications, certain parts of 
                    #' the design may have more influence than others over the estimate of interest,
                    #' or have a larger than desired effect. A summary of the DFBETA diagnostic 
                    #' is provided.
                    #'  
                    #' @param type One of either `sim_data` (recalls saved data from a previous
                    #' call to this function), `sim` (full simulation using MCML), or `sim_approx` (
                    #' uses GLS to approximate the MLE and standard errors)
                    #' @param iter Integer. The number of iterations of the simulation to run
                    #' @param par Integer. The parameter of interest for which design analysis 
                    #' statistics should be calculated. Refers to the column of X.
                    #' @param alpha Numeric. The type I error rate.
                    #' @param sim_design Optional. A different `Design` object that will be used to 
                    #' simulate data to allow for model misspecification.
                    #' @param parallel Logical indicating whether to run the simulations in parallel
                    #' @param verbose Logical indicating whether to report detailed output. Defaults to TRUE.
                    #' @param ... Additional arguments passed to `MCML`, see below.
                    #' @return A `glmmr.sim` object containing the estimates from all the simulations, including
                    #' standard errors, deletion diagnostic statistics, and details about the simulation.
                    #' @seealso \link[glmmr]{print.glmmr.sim}
                    #' @examples 
                    #' ...
                    analysis = function(type, 
                                        iter,
                                        par,
                                        alpha = 0.05,
                                        sim_design,
                                        parallel,
                                        verbose = TRUE,
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
                        
                        
                        sim_mean_formula <- ifelse(missing(sim_design),NA,as.character(sim_design$mean_function$formula))
                        sim_cov_formula <- ifelse(missing(sim_design),NA,as.character(sim_design$covariance$formula))
                        
                        
                        res <- list(
                          coefficients = lapply(out,function(i)i[[1]]$coefficients),
                          dfbeta = lapply(out,function(i)i[[2]]),
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
                          sim_mean_formula = sim_mean_formula,
                          sim_cov_formula = sim_cov_formula,
                          family = self$mean_function$family,
                          aic = lapply(out,function(i)i[[1]]$aic),
                          Rsq = lapply(out,function(i)i[[1]]$Rsq),
                          n = self$n(),
                          par = par,
                          type = 1
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
                          dfbeta = lapply(out,function(i)i[[2]]),
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
                          aic = lapply(out,function(i)i[[1]]$aic),
                          Rsq = lapply(out,function(i)i[[1]]$Rsq),
                          n = self$n(),
                          par = par,
                          type = 1
                        )
                        
                        class(res) <- "glmmr.sim"
                        
                        if(verbose)message("saving simulation data")
                        private$saved_sim_data <- res
                      }
                      if(type=="sim_data")res <- private$saved_sim_data
                      
                      invisible(res)
                    },
                    #' @description 
                    #' Approximate power of the design using the GLS variance
                    #' @details 
                    #' Calculates the approximate power of the design using the square root
                    #' of the relevant element of the GLS variance matrix:
                    #' 
                    #'  \deqn{(X^T\Sigma^{-1}X)^{-1}}
                    #'  
                    #' Note that this is equivalent to using the "design effect" for many
                    #' models.
                    #' @param par Integer indicating which parameter of the design the power
                    #' should be calculated for. Refers to the order of parameters and column
                    #' of X
                    #' @param value Numeric specifying the value of the parameter to calculate
                    #' the power at
                    #' @param alpha Numeric between zero and one indicating the type I error rate. 
                    #' Default of 0.05.
                    #' @return A value between zero and one indicating the approximate power of the
                    #' design.
                    #' @examples 
                    #' ...
                    power = function(par,
                                     value,
                                     alpha=0.05){
                      self$check(verbose=FALSE)
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
                    #' @description
                    #' Subsets the design by removing the specified observations
                    #' 
                    #' Given a vector of row indices, the corresponding rows will be removed from the 
                    #' mean function and covariance
                    #' @param index Integer or vector integers listing the rows to remove from the design
                    #' @return The function updates the object and nothing is returned
                    #' @examples
                    #' ...
                    subset_rows = function(index){
                      self$mean_function$subset_rows(index)
                      self$covariance$subset(index)
                    },
                    #' @description 
                    #' Removes a column or columns from the X matrix 
                    #' 
                    #' Removes the specified columns from the linked mean function object's X matrix
                    #' @param index Integer or vector of integers specifying the indexes of the columns to remove from X
                    #' @return The function updates the object and nothing is returned
                    #' @examples
                    #' ...
                    subset_cols = function(index){
                      self$mean_function$subset_cols(index)
                    },
                    #'@description 
                    #'Generates a plot of the design
                    #'
                    #'Generates a 'bubble' plot of the design with respect to two or three variables 
                    #'in which the size of the points at each location are scaled by the number of observations at that
                    #'location. For example, for a cluster randomised trial the user might specify
                    #'time period on the x-axis and cluster ID on the y-axis. For a geospatial
                    #'sampling design the x and y axes might represent spatial dimensions.
                    #'@param x String naming a column in the data frame in the linked covariance object (self$covariance$data) 
                    #'for the x-axis
                    #'@param y  String naming a column in the data frame in the linked covariance object (self$covariance$data) 
                    #'for the y-axis
                    #'@param z Optional. String naming a column in the data frame in the linked covariance object (self$covariance$data) 
                    #'for a 'third axis' used for faceting
                    #'@param treat String naming a column in the data frame in the linked mean function
                    #'object (self$mean_function$data) that identifies the treatment status of the observations
                    #'at each location, used to set the colour of the points in the plot
                    #'@return A \link[ggplot2]{ggplot2} plot
                    #'@examples
                    #' ...
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
                    #'@description 
                    #'Generates a realisation of the design
                    #'
                    #'Generates a single vector of outcome data based upon the 
                    #'specified GLMM design
                    #'@param type Either 'y' to return just the outcome data, or 'data'
                    #' to return a data frame with the simulated outcome data alongside the model data 
                    #' @return Either a vector or a data frame
                    #' @examples
                    #' ...
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
                    #'@description
                    #'Checks for any changes in linked objects and updates
                    #'@param verbose Logical indicating whether to report if any updates are made, defaults to TRUE
                    #'@return Linked objects are updated by nothing is returned
                    #'@examples
                    #'...
                    check = function(verbose=TRUE){
                      self$covariance$check(verbose=verbose)
                      self$mean_function$check(verbose = verbose)
                      if(private$hash != private$hash_do()){
                        private$generate()
                      }
                    },
                    # apv = function(prior,
                    #                var,
                    #                prior.fun,
                    #                iter,
                    #                verbose=TRUE){
                    #   if(verbose)message("Monte Carlo integration")
                    #   samps <- pbapply::pbreplicate(iter,self$posterior(prior,var,do.call(prior.fun,list())))
                    #   summary(samps)
                    # },
                    #'@description
                    #'Markov Chain Monte Carlo Maximum Likelihood  model fitting
                    #'
                    #'@details
                    #'Fits generalised linear mixed models using one of three algorithms: Markov Chain Newton
                    #'Raphson (MCNR), Markov Chain Expectation Maximisation (MCEM), or Maximum simulated
                    #'likelihood (MSL). All the algorithms are described by McCullagh (1997). For each iteration
                    #'of the algorithm the unobserved random effect terms (\eqn{\gamma}) are simulated
                    #'using Markov Chain Monte Carlo (MCMC) methods (we use Hamiltonian Monte Carlo through Stan),
                    #'and then these values are conditioned on in the subsequent steps to estimate the covariance
                    #'parameters and the mean function parameters (\eqn{\beta}). For all the algorithms, 
                    #'the covariance parameter estimates are updated using an expectation maximisation step.
                    #'For the mean function parameters you can either use a Newton Raphson step (MCNR) or
                    #'an expectation maximisation step (MCEM). A simulated likelihood step can be added at the 
                    #'end of either MCNR or MCEM, which uses an importance sampling technique to refine the 
                    #'parameter estimates. 
                    #'
                    #'The accuracy of the algorithm depends on the user specified tolerance. For higher levels of
                    #'tolerance, larger numbers of MCMC samples are likely need to sufficiently reduce Monte Carlo error.
                    #'
                    #'The function also offers different methods of obtaining standard errors. First, one can generate
                    #'estimates from the estimated Hessian matrix (`se.method = 'lik'`). Second, there are robust standard 
                    #'errors using a sandwich estimator based on White (1982) (`se.method = 'robust'`). 
                    #'Third, there are use approximate generalised least squares estimates based on the maximum likelihood 
                    #'estimates of the covariance
                    #'parameters (`se.method = 'approx'`), or use a permutation test approach (`se.method = 'perm'`).
                    #'Note that the permutation test can be accessed separately with the function `permutation_test()`.
                    #'
                    #'There are several options that can be specified to the function using the `options` argument. 
                    #'The options should be provided as a list, e.g. `options = list(method="mcnr")`. The possible options are:
                    #'* `b_se_only` TRUE (calculate standard errors of the mean function parameters only) or FALSE (calculate
                    #'all standard errors), default it FALSE.
                    #'* `use_cmdstanr` TRUE (uses `cmdstanr` for the MCMC sampling, requires cmdstanr), or FALSE (uses `rstan`). Default is FALSE.
                    #'* `skip_cov_optim` TRUE (skips the covariance parameter estimation step, and uses the values covariance$parameters), or 
                    #'FALSE (run the whole algorithm)], default is FALSE
                    #'* `method` One of either `mcnr` or `mcem`, see above. Default is `mcnr`.
                    #'* `sim_lik_step` TRUE (conduct a simulated likelihood step at the end of the algorithm), or FALSE (does
                    #'not do this step), defaults to FALSE.
                    #'* `no_warnings` TRUE (do not report any warnings) or FALSE (report warnings), default to FALSE
                    #'* `perm_type` Either `cov` (use weighted test statistic in permutation test) or `unw` (use unweighted
                    #' test statistic), defaults to `cov`. See `permutation_test()`.
                    #' * `perm_iter` Number of iterations for the permutation test, default is 100.
                    #' * `perm_parallel` TRUE (run permuation test in parallel) or FALSE (runs on a single thread), default to TRUE
                    #' * `warmup_iter` Number of warmup iterations on each iteration for the MCMC sampler, default is 500
                    #' * `perm_ci_steps` Number of steps for the confidence interval search procedure if using the permutation
                    #' test, default is 1000. See `permutation_test()`.
                    #' * `fd_tol` The tolerance of the first difference method to estimate the Hessian and Gradient, default 
                    #' is 1e-4.
                    #'
                    #'@param y A numeric vector of outcome data
                    #'@param start Optional. A numeric vector indicating starting values for the MCML algorithm iterations. 
                    #'If this is not specified then the parameter values stored in the linked mean function object will be used.
                    #'@param se.method One of either `'lik'`, `'approx'`, `'perm'`, or `'none'`, see Details.
                    #'@param verbose Logical indicating whether to provide detailed output, defaults to TRUE.
                    #'@param tol Numeric value, tolerance of the MCML algorithm, the maximum difference in parameter estimates 
                    #'between iterations at which to stop the algorithm.
                    #'@param m Integer. The number of MCMC draws of the random effects on each iteration.
                    #'@param max.iter Integer. The maximum number of iterations of the MCML algorithm.
                    #'@param options An optional list providing options to the algorithm, see Details.
                    #'@return A `mcml` object
                    #'@examples
                    #'...
                    #'@md
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
                      if(!se.method%in%c("perm","lik","none","robust"))stop("se.method should be 'perm', 'lik', 'robust', or 'none'")
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
                      method <- ifelse("method"%in%names(options),options$method,"mcnr")
                      sim_lik_step <- ifelse("sim_lik_step"%in%names(options),options$sim_lik_step,FALSE)
                      no_warnings <- ifelse("no_warnings"%in%names(options),options$no_warnings,FALSE)
                      perm_type <- ifelse("perm_type"%in%names(options),options$perm_type,"cov")
                      perm_iter <- ifelse("perm_iter"%in%names(options),options$perm_iter,100)
                      perm_parallel <- ifelse("perm_parallel"%in%names(options),options$perm_iter,TRUE)
                      warmup_iter <- ifelse("warmup_iter"%in%names(options),options$warmup_iter,500)
                      perm_ci_steps <- ifelse("perm_ci_steps"%in%names(options),options$perm_ci_steps,1000)
                      fd_tol <- ifelse("fd_tol"%in%names(options),options$fd_tol,1e-4)
                      block <- ifelse("block"%in%names(options),options$block,TRUE)
                      
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
                          dsamps <- t(dsamps)
                        }
                        
                        
                        # BETA PARAMETERS STEP
                        if(method == "mcnr"){
                          beta_step <- mcnr_step(y = y,
                                                 X= as.matrix(self$mean_function$X),
                                                 Z = as.matrix(self$covariance$Z),
                                                 beta = theta[parInds$b],
                                                 u = dsamps,
                                                 family = self$mean_function$family[[1]],
                                                 link = self$mean_function$family[[2]])
                          
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
                          newtheta <- do.call(d_lik_optim,append(self$covariance$.__enclos_env__$private$D_data,
                                                                 list(u = dsamps,
                                                                      start = c(theta[parInds$cov]),
                                                                      lower= rep(1e-6,length(parInds$cov)),
                                                                      upper= rep(Inf,length(parInds$cov)),
                                                                      trace=0,
                                                                      block = block)))
                          theta[parInds$cov] <- drop(newtheta)
                        }
                        
                        if(verbose)cat("\ntheta:",theta[all_pars])
                      }
                      
                      not_conv <- iter >= max.iter|any(abs(theta-thetanew)>tol)
                      if(not_conv&!no_warnings)warning(paste0("algorithm not converged. Max. difference between iterations :",max(abs(theta-thetanew)),". Suggest 
                                                 increasing m, or trying a different algorithm."))
                      
                      if(sim_lik_step){
                        if(verbose)cat("\n\n")
                        if(verbose)message("Optimising simulated likelihood")
                        newtheta <- do.call(f_lik_optim,append(self$covariance$.__enclos_env__$private$D_data,
                                                               list(as.matrix(self$covariance$Z),
                                                                    as.matrix(self$mean_function$X),
                                                                    y,
                                                                    dsamps,
                                                                    theta[parInds$cov],
                                                                    family=self$mean_function$family[[1]],
                                                                    link=self$mean_function$family[[2]],
                                                                    start = theta[all_pars],
                                                                    lower = c(rep(-Inf,P),rep(1e-5,length(all_pars)-P)),
                                                                    importance = TRUE)))
                        theta[all_pars] <- newtheta
                      }
                      
                      if(verbose)cat("\n\nCalculating standard errors...")
                      
                      if(family%in%c("gaussian")){
                        mf_pars <- theta[c(parInds$b,parInds$sig)]
                        mf_pars_names <- c(colnames(self$mean_function$X),"sigma")
                      } else {
                        mf_pars <- theta[c(parInds$b)]
                        mf_pars_names <- colnames(self$mean_function$X)
                      }
                      
                      cov_pars_names <- rep(as.character(unlist(rev(self$covariance$.__enclos_env__$private$flist))),
                                            drop(self$covariance$.__enclos_env__$private$D_data$N_par))#paste0("cov",1:R)
                      permutation <- FALSE
                      robust <- FALSE
                      if(se.method=="lik"|se.method=="robust"){
                        if(verbose&!robust)cat("using Hessian\n")
                        if(verbose&robust)cat("using robust sandwich estimator\n")
                        
                          hess <- tryCatch(do.call(f_lik_hess,append(self$covariance$.__enclos_env__$private$D_data,
                                                                     list(as.matrix(self$covariance$Z),
                                                                          as.matrix(self$mean_function$X),
                                                                          y,
                                                                          dsamps,
                                                                          theta[parInds$cov],
                                                                          family=self$mean_function$family[[1]],
                                                                          link=self$mean_function$family[[2]],
                                                                          start = theta[all_pars],
                                                                          lower = c(rep(-Inf,P),rep(1e-5,length(all_pars)-P)),
                                                                          upper = c(rep(Inf,P),rep(Inf,length(all_pars)-P)),
                                                                          tol=fd_tol))),
                                           error=function(e)NULL)
                          
                          hessused <- TRUE
                          
                          semat <- tryCatch(Matrix::solve(hess),error=function(e)NULL)
                      
                          if(se.method == "robust"&!is.null(semat)){
                            hlist <- list()
                            #identify the clustering and sum over independent clusters
                            D_data <- self$covariance$.__enclos_env__$private$D_data
                            gr_var <- apply(D_data$func_def,1,function(x)any(x==1))
                            gr_count <- D_data$N_dim
                            gr_id <- which(gr_count == min(gr_count[gr_var]))
                            gr_cov_var <- D_data$cov_data[1:D_data$N_dim[gr_id],
                                                                   1:D_data$N_var_func[gr_id,which(D_data$func_def[gr_id,]==1)],gr_id,drop=FALSE]
                            gr_cov_var <- as.data.frame(gr_cov_var)
                            colnames(gr_cov_var) <- all.vars(rev(self$covariance$.__enclos_env__$private$flist)[[gr_id]])
                            Z_in <- match_rows(self$covariance$data,as.data.frame(gr_cov_var),by=colnames(gr_cov_var))
                            
                            for(i in 1:ncol(Z_in)){
                              id_in <- which(Z_in[,i]==1)
                              g1 <- matrix(0,nrow=length(all_pars),ncol=1)
                              g1 <- do.call(f_lik_grad,append(self$covariance$.__enclos_env__$private$D_data,
                                                              list(as.matrix(self$covariance$Z)[id_in,,drop=FALSE],
                                                                   as.matrix(self$mean_function$X)[id_in,,drop=FALSE],
                                                                   y[id_in],
                                                                   dsamps,
                                                                   theta[parInds$cov],
                                                                   family=self$mean_function$family[[1]],
                                                                   link=self$mean_function$family[[2]],
                                                                   start = theta[all_pars],
                                                                   lower = c(rep(-Inf,P),rep(1e-5,length(all_pars)-P)),
                                                                   upper = c(rep(Inf,P),rep(Inf,length(all_pars)-P)),
                                                                   tol=fd_tol)))
                              
                              hlist[[i]] <- g1%*%t(g1)
                            }
                            g0 <- Reduce('+',hlist)
                            semat <- semat%*%g0%*%semat
                            robust <- TRUE
                          }
                          
                          if(!is.null(semat)){
                            SE <- tryCatch(sqrt(Matrix::diag(semat)),
                                           error=function(e)rep(NA,length(mf_pars)+length(cov_pars_names)))
                          } else {
                            SE <- rep(NA,length(mf_pars)+length(cov_pars_names))
                          }
                          
                          res <- data.frame(par = c(mf_pars_names,cov_pars_names,paste0("d",1:Q)),
                                            est = c(mf_pars,theta[parInds$cov],rowMeans(dsamps)),
                                            SE=c(SE,apply(dsamps,1,sd)))
                        
                        
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
                      robust <- FALSE
                      } else {
                        res <- data.frame(par = c(mf_pars_names,cov_pars_names),
                                          est = c(mf_pars,theta[parInds$cov]),
                                          SE=NA,
                                          lower = NA,
                                          upper =NA)
                        hessused <- FALSE
                        robust <- FALSE
                      }
                      
                     rownames(dsamps) <- Reduce(c,rev(cov1$.__enclos_env__$private$flistlabs))
                     
                     ## model summary statistics
                     aic_data <- append(list(Z = as.matrix(self$covariance$Z),
                                             X = as.matrix(self$mean_function$X),
                                             y = y,
                                             u = dsamps,
                                             family = self$mean_function$family[[1]],
                                             link=self$mean_function$family[[2]]), 
                                        self$covariance$.__enclos_env__$private$D_data)
                     aic <- do.call(aic_mcml,append(aic_data,list(beta_par = mf_pars,
                                                                  cov_par = theta[parInds$cov])))
                     
                     xb <- self$mean_function$X %*% theta[parInds$b]
                     zd <- self$covariance$Z %*% rowMeans(dsamps)
                     
                     wdiag <- gen_dhdmu(Matrix::drop(xb),
                                        family=self$mean_function$family[[1]],
                                        link = self$mean_function$family[[2]])
                     
                     if(self$mean_function$family[[1]]%in%c("gaussian","gamma")){
                       wdiag <- theta[parInds$sig] * wdiag
                     }
                     
                     total_var <- var(Matrix::drop(xb)) + var(Matrix::drop(zd)) + mean(wdiag)
                     condR2 <- (var(Matrix::drop(xb)) + var(Matrix::drop(zd)))/total_var
                     margR2 <- var(Matrix::drop(xb))/total_var
                     
                     
                     out <- list(coefficients = res,
                                 converged = !not_conv,
                                 method = method,
                                 hessian = hessused,
                                 robust = robust,
                                 permutation = permutation,
                                 m = m,
                                 tol = tol,
                                 sim_lik = sim_lik_step,
                                 aic = aic,
                                 Rsq = c(cond = condR2,marg=margR2),
                                 mean_form = as.character(self$mean_function$formula),
                                 cov_form = as.character(self$covariance$formula),
                                 family = self$mean_function$family[[1]],
                                 link = self$mean_function$family[[2]],
                                 re.samps = dsamps)
                     
                     class(out) <- "mcml"
                      
                      self$mean_function$parameters <- orig_par_b 
                      self$covariance$parameters <- orig_par_cov
                      #self$check(verbose=FALSE)
                      
                      return(out)
                    },
                    #'@description
                    #'Calculates DFBETA deletion diagnostic values
                    #'
                    #'@details
                    #'Add more description...
                    #'@param y Numeric vector of outcomes data
                    #'@param par Integer indicating which parameter, as a column of X, to report the 
                    #'deletion diagnostics for.
                    #'@return A vector of length self$n() with the value of DFBETA for each observation for
                    #'the specified parameter
                    #'@examples
                    #'...
                    dfbeta = function(y,
                                      par){
                        dfbeta_stat(as.matrix(self$Sigma),
                                    as.matrix(self$mean_function$X),
                                    y,
                                    par)
                      
                    },
                    # approx_posterior = function(prior,
                    #                      var,
                    #                      parameters){
                    #   #move to private and set this as Monte Carlo integration
                    #   #can just request a function that outputs a new set of covariance parameters
                    #   R <- solve(Matrix::Matrix(diag(prior)))
                    #   S <- private$genS(self$covariance$sampleD(parameters),self$covariance$Z,private$W,update=FALSE)
                    #   M <- R + Matrix::crossprod(self$mean_function$X,solve(S))%*%self$mean_function$X
                    #   M <- solve(M)
                    #   M[var,var]
                    # },
                    #'@description
                    #'Conducts a permuation test
                    #'
                    #' Estimates p-values and confidence intervals using a permutation test
                    #'
                    #'@details
                    #' If the user provided a re-randomisation function to the linked mean function object (see \link[glmmr]{MeanFunction}),
                    #' then a permuation test can be conducted. A new random assignment is generated on each iteration of the permutation test.
                    #' The test statistic can be either a quasi-score statistic, weighting the observations using the covariance matrix (`type="cov"`),
                    #' or an unweighted statistic that weights each observation in each cluster equally (`type="unw"`). The 1-alpha% 
                    #' confidence interval limits are estimated using an efficient iterative stochastic search procedure. On each step of the algorithm
                    #' a single permuation and test statistic is generated, and the current estimate of the confidence interval limit either 
                    #' increased or decreased depedning on its value. The procedure converges in probability to the true limits, see Watson et al (2021)
                    #' and Garthwaite (1996).
                    #' 
                    #'@param y Numeric vector of outcome data
                    #'@param permutation.par Integer indicator which parameter to conduct a permutation
                    #'test for. Refers to a column of the X matrix.
                    #'@param start Value of the parameter. Used both as a starting value for the algorithms
                    #'and as a best estimate for the confidence interval search procedure.
                    #'@param iter Integer. Number of iterations of the permuation test to conduct
                    #'@param nsteps Integer. Number of steps of the confidence interval search procedure
                    #'@param type Either `cov` for a test statistic weighted by the covariance matrix, or 
                    #'`unw` for an unweighted test statistic. See Details.
                    #'@param parallel Logical indicating whether to run the tests in parallel
                    #'@param verbose Logical indicating whether to report detailed output
                    #'@return A list with the estimated p-value and the estimated lower and upper 95% confidence interval
                    #' @references 
                    #' Watson et al. Arxiv
                    #' Braun and Feng
                    #' Gail
                    #' Garthwaite
                    #'@examples
                    #'...
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
                      #invS <- ifelse(type=="cov",Matrix::solve(self$Sigma),1)
                      if(w.opt){
                        invS <- Matrix::solve(self$Sigma)
                      } else {
                        invS <- 1
                      }
                      
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
                    },
                    #' @description 
                    #' Fit the GLMM using MCMC
                    #' 
                    #' Fit the GLMM using MCMC. The function calls a Stan program to draw posterior samples.
                    #' @details 
                    #' Draws samples from the posterior distribution of the model parameters using Stan. Priors are specified using the `priors` argument.
                    #' Currently, only Gaussian (or half-Gaussian for covariance parameters) prior distributions are supported. The argument `priors`
                    #' accepts a list with three or four elements, `prior_b_mean`, `prior_b_sd` which are vectors specifying the prior mean and 
                    #' standard deviation for the mean function parameters \eqn{\beta} in the model; `prior_g_sd` specifying the prior standard deviation
                    #' for the half-Gaussian prior for the covariance parameters, and optionally `sigma_sd` for the half-Gaussian prior for the scale 
                    #' terms in models that have a scale parameter (Gaussian and Gamma currently). By default the function uses `rstan`, however the function
                    #' can optionally call `cmdstanr` instead if it is installed, using the option `use_cmdstanr=TRUE`. For further details on the use and
                    #' arguments that can be used, see \link[rstan]{sampling}.
                    #' @param y Numeric vector providing the outcome data
                    #' @param priors A named list specifying the prior mean and standard deviation for the Gaussian prior distributions, see Details.
                    #' @param warmup_iter A positive integer specifying the number of warmup iterations for each MCMC chain
                    #' @param sampling_iter A positive integer specifying the number of sampling iterations for each MCMC chain
                    #' @param chain A positive integer specifying the number of MCMC chains to run
                    #' @param use_cmdstanr Logical indicating whether to use `cmdstanr`, the default is FALSE, which will use `rstan`
                    #' @param ... Additional arguments to pass to \link[rstan]{sampling}.
                    #' @return Either an S4 `stanfit` object returned by \link[rstan]{sampling}, or a `cmdstanr` environment, depending on the sampler used.
                    #' @seealso \link[rstan]{sampling}, \link[rstan]{stan}
                    #' @examples 
                    #' ...
                    MCMC = function(y,
                                    priors,
                                    warmup_iter = 1000,
                                    sampling_iter = 1000,
                                    chains = 4,
                                    use_cmdstanr=FALSE,...){
                      
                      fn1 <- Reduce(cbind,self$covariance$.__enclos_env__$private$Funclist)
                      if(any(fn1[1,]==5))warning("Matern not implemented exactly in Stan, order one is assumed as 
fixed for the modified Bessel function of the second kind.")
                      
                      if(missing(priors) || !is("priors",list)){
                        warning("priors not specified properly, using defaults")
                        priors <- list(
                          prior_b_mean = rep(0,ncol(self$mean_function$X)),
                          prior_b_sd = rep(1,ncol(self$mean_function$X)),
                          prior_g_sd = rep(1,length(self$covariance$parameters)),
                          sigma_sd = 1
                        )
                      } else {
                        if(!"prior_b_mean" %in%names(priors)){
                          message("prior_b_mean not specified, setting to 0")
                          priors$prior_b_mean <- rep(0,ncol(self$mean_function$X))
                        }
                        if(!"prior_b_sd" %in%names(priors)){
                          message("prior_b_sd not specified, setting to 1")
                          priors$prior_b_sd <- rep(1,ncol(self$mean_function$X))
                        }
                        if(!"prior_g_sd" %in%names(priors)){
                          message("prior_g_sd not specified, setting to 1")
                          priors$prior_g_sd <- rep(1,length(self$covariance$parameters))
                        }
                        if(self$mean_function$family[[1]] %in% c("gaussian","gamma")){
                          if(!"sigma_sd" %in%names(priors)){
                            message("sigma_sd not specified, setting to 1")
                            priors$sigma_sd <- 1
                          }
                        }
                      }
                      
                      if(use_cmdstanr){
                        if(!requireNamespace("cmdstanr"))stop("cmdstanr not available")
                        model_file <- system.file("stan",
                                                  paste0(self$mean_function$family[[1]],".stan"),
                                                  package = "glmmr",
                                                  mustWork = TRUE)
                        mod <- cmdstanr::cmdstan_model(model_file)
                      }
                      
                      ## get data
                      stan_data <- private$gen_stan_data(type = "y",
                                                         y=y)
                      stan_data$prior_b_mean <- priors$prior_b_mean
                      stan_data$prior_b_sd <- priors$prior_b_sd
                      stan_data$prior_g_sd <- priors$prior_g_sd
                      if(self$mean_function$family[[1]] %in% c("gaussian","gamma")){
                        stan_data$sigma_sd <- priors$sigma_sd
                      }
                      
                      
                      if(use_cmdstanr){
                        fit <- mod$sample(data = stan_data,
                                          chains = chains,
                                          iter_warmup = warmup_iter,
                                          iter_sampling = sampling_iter,
                                          ...)
                      } else {
                        fit <- rstan::sampling(stanmodels[[self$mean_function$family[[1]]]],
                                               data = stan_data,
                                               chains = chains,
                                               warmup = warmup_iter,
                                               iter = warmup_iter+sampling_iter,...)
                      }
                      return(fit)
                    },
                    #' @description 
                    #' Bayesian design analysis
                    #' 
                    #' Runs a Bayesian design analysis using simulated data
                    #' @details 
                    #' The Bayesian design analysis conducts multiple iterations of simulating data from a GLMM and then fitting a GLMM to analyse
                    #' the properties of the posterior distribution and model calibration to inform study design. Data are either simulated from the same
                    #' model as the analysis model, i.e. there is correct model specification (nothing is specified for the option `sim_design`), or data
                    #' are simulated from a different model (a `Design` object passed to the `sim_design` argument of this function) and then fit using
                    #' the `Design` object from which this function was called, i.e. with model misspecification. The function analyses four related statistics 
                    #' summarising the design. First, the expected posterior variance; second, the variance of the posterior variance; third, the probability the
                    #' parameter of interest will be greater than some threshold; and fourth, simulation-based calibration (see Talts et al, 2021).
                    #' 
                    #' Model fitting is done using Stan's sampler (see \link[rstan]{sampling} and the `MCMC()` function in this class).
                    #' @param iter A positive integer specifying the number of simulation iterations to run
                    #' @param par A positive interger specifying the index of the parameter in the mean function to analyse, refers specifically to the column 
                    #' of X
                    #' @param priors A named list specifying the priors, see Details in the `MCMC()` function in this class
                    #' @param threshold A number specifying the threshold. The probability that the parameter of interest is greater than this threshold will be calculated.
                    #' @param sim_design. Optional. A different `Design` object that should be used to simulate the data for the simulation. 
                    #' @param priors_sim Optional. A named list of the same structure as `priors` the specifies the priors from which to simulate data, if they are different
                    #' to the priors for model fitting.
                    #' @param warmup_iter A positive integer specifying the number of warmup iterations for each MCMC chain when fitting the model on each simulation
                    #' iteration.
                    #' @param sampling_iter A positive integer specifying the number of sampling iterations for each MCMC chain when fitting the model on each simulation
                    #' iteration.
                    #' @param chains A positive integer specifying the number of MCMC chains when fitting the model on each simulation
                    #' iteration.
                    #' @param verbose Logical indicating whether to provide detailed output of the progress of simulations.
                    #' @param use_cmdstanr Logical indicating whether to use `cmdstanr` instead of `rstan`, the default is FALSE (i.e. use `rstan`)
                    #' @return A `glmmr.sim` object containing the samples of posterior variance and relevant probabilities to summarise the analysis
                    #' @references 
                    #' ...
                    #' @examples 
                    #' ...
                    analysis_bayesian = function(iter,
                                                 par,
                                                 priors,
                                                 threshold = 0,
                                                 sim_design,
                                                 priors_sim,
                                                 warmup_iter = 200,
                                                 sampling_iter = 200,
                                                 chains=1,
                                                 verbose=TRUE,
                                                 use_cmdstanr=FALSE,...){
                      if(missing(priors) || !is("priors",list)){
                        warning("priors not specified properly, using defaults")
                        priors <- list(
                          prior_b_mean = rep(0,ncol(self$mean_function$X)),
                          prior_b_sd = rep(1,ncol(self$mean_function$X)),
                          prior_g_sd = rep(1,length(self$covariance$parameters)),
                          sigma_sd = 1
                        )
                      } else {
                        if(!"prior_b_mean" %in%names(priors)){
                          if(verbose)message("prior_b_mean not specified, setting to 0")
                          priors$prior_b_mean <- rep(0,ncol(self$mean_function$X))
                        }
                        if(!"prior_b_sd" %in%names(priors)){
                          if(verbose)message("prior_b_sd not specified, setting to 1")
                          priors$prior_b_sd <- rep(1,ncol(self$mean_function$X))
                        }
                        if(!"prior_g_sd" %in%names(priors)){
                          if(verbose)message("prior_g_sd not specified, setting to 1")
                          priors$prior_g_sd <- rep(1,length(self$covariance$parameters))
                        }
                        if(self$mean_function$family[[1]] %in% c("gaussian","gamma")){
                          if(!"sigma_sd" %in%names(priors)){
                            message("sigma_sd not specified, setting to 1")
                            priors$sigma_sd <- 1
                          }
                        }
                      }
                      
                      f_str <- ifelse(missing(sim_design),"_sim","_sim_misspec")
                      
                      if(use_cmdstanr){
                        if(!requireNamespace("cmdstanr"))stop("cmdstanr not available")
                        
                        model_file <- system.file("stan",
                                                  paste0(self$mean_function$family[[1]],f_str,".stan"),
                                                  package = "glmmr",
                                                  mustWork = TRUE)
                        mod <- cmdstanr::cmdstan_model(model_file)
                      }
                      
                      stan_data <- private$gen_stan_data(type = "s")
                      stan_data$prior_b_mean <- priors$prior_b_mean
                      stan_data$prior_b_sd <- priors$prior_b_sd
                      stan_data$prior_g_sd <- priors$prior_g_sd
                      if(self$mean_function$family[[1]] %in% c("gaussian","gamma")){
                        stan_data$sigma_sd <- priors$sigma_sd
                      }
                      
                      if(!missing(sim_design)){
                        if(!is(sim_design,"Design"))stop("sim_design must be another design")
                        
                        if(missing(priors_sim)|!is("priors",list)){
                          warning("priors_sim not specified properly, using defaults")
                          priors_sim <- list(
                            prior_b_mean = rep(0,ncol(sim_design$mean_function$X)),
                            prior_b_sd = rep(1,ncol(sim_design$mean_function$X)),
                            prior_g_sd = rep(1,length(sim_design$covariance$parameters)),
                            sigma_sd = 1
                          )
                        } else {
                          if(!"prior_b_mean" %in%names(priors)){
                            message("prior_b_mean not specified in priors_sim, setting to 0")
                            priors_sim$prior_b_mean <- rep(0,ncol(sim_design$mean_function$X))
                          }
                          if(!"prior_b_sd" %in%names(priors)){
                            message("prior_b_sd not specified in priors_sim, setting to 1")
                            priors_sim$prior_b_sd <- rep(1,ncol(sim_design$mean_function$X))
                          }
                          if(!"prior_g_sd" %in%names(priors)){
                            message("prior_g_sd not specified in priors_sim, setting to 1")
                            priors_sim$prior_g_sd <- rep(1,length(sim_design$covariance$parameters))
                          }
                          if(self$mean_function$family[[1]] %in% c("gaussian","gamma")){
                            if(!"sigma_sd" %in%names(priors)){
                              message("sigma_sd not specified in priors_sim, setting to 1")
                              priors_sim$sigma_sd <- 1
                            }
                          }
                        }
                        
                        stan_data_m <- sim_design$.__enclos_env__$private$gen_stan_data(type = "m")
                        stan_data_m$prior_b_mean_m <- priors_sim$prior_b_mean
                        stan_data_m$prior_b_sd_m <- priors_sim$prior_b_sd
                        stan_data_m$prior_g_sd_m <- priors_sim$prior_g_sd
                        if(sim_design$mean_function$family[[1]] %in% c("gaussian","gamma")){
                          stan_data_m$sigma_sd <- priors_sim$sigma_sd
                        }
                        
                        stan_data <- append(stan_data,stan_data_m)
                      }
                      
                      stan_data$par_ind <- par
                      stan_data$threshold <- threshold
                      
                      posterior_var <- c()
                      sbc_ranks <- c()
                      posterior_threshold <- c()
                      if(verbose)cat("\nStarting simuations and model fitting...\n")
                      if(verbose)cat(progress_bar(0,iter))
                      for(i in 1:iter){
                        if(use_cmdstanr){
                          fit <- mod$sample(data = stan_data,
                                            chains = chains,
                                            iter_warmup = warmup_iter,
                                            iter_sampling = sampling_iter,
                                            ...)
                        } else {
                          fit <- rstan::sampling(stanmodels[[paste0(self$mean_function$family[[1]],f_str)]],
                                                 data = stan_data,
                                                 chains = chains,
                                                 warmup = warmup_iter,
                                                 iter = warmup_iter+sampling_iter,...)
                          
                          b1 <- rstan::extract(fit,"beta")
                          b1 <- b1$beta[,par]
                          posterior_var[i] <- var(b1)
                          b2 <- rstan::extract(fit,"betaIn")
                          sbc_ranks[i] <- sum(b2$betaIn)
                          b3 <- rstan::extract(fit,"betaThresh")
                          posterior_threshold[i] <- mean(b3$betaThresh)
                        }
                        
                        if(verbose)cat("\rSimulations progress: ",progress_bar(i,iter))
                      }
                      
                      sim_mean_formula <- ifelse(missing(sim_design),NA,as.character(sim_design$mean_function$formula))
                      sim_cov_formula <- ifelse(missing(sim_design),NA,as.character(sim_design$covariance$formula))
                      psim <- ifelse(missing(priors_sim),NA,priors_sim)
                      
                      res <- list(posterior_var = posterior_var,
                                  sbc_ranks=sbc_ranks,
                                  posterior_threshold=posterior_threshold,
                                  iter = iter,
                                  priors = priors,
                                  priors_sim = psim,
                                  mean_formula = self$mean_function$formula,
                                  cov_formula = self$covariance$formula,
                                  sim_mean_formula = sim_mean_formula,
                                  sim_cov_formula = sim_cov_formula,
                                  family = self$mean_function$family,
                                  type = 2)
                      
                      class(res) <- "glmmr.sim"
                      ## add on other return, include model specification and misspecification
                      return(res)
                    }
                  ),
                  private = list(
                    W = NULL,
                    Xb = NULL,
                    logit = function(x){
                      exp(x)/(1+exp(x))
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
                    genW = function(family,
                                    Xb,
                                    var_par=NULL){
                      # assume random effects value is at zero
                      if(!family[[1]]%in%c("poisson","binomial","gaussian","gamma"))stop("family must be one of Poisson, Binomial, Gaussian, Gamma")
                      
                      wdiag <- gen_dhdmu(c(Xb),
                                         family=family[[1]],
                                         link = family[[2]])
                      
                      if(family[[1]]%in%c("gaussian","gamma")){
                        wdiag <- var_par * wdiag
                      }
                      W <- diag(drop(wdiag))
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
                      dfb <- self$dfbeta(y=ysim,
                                         par = par)
                        
                        
                      
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
                      dfb <- self$dfbeta( y=ysim,
                                         par = par)
                      
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
                      if(w.opt){
                        invS <- Matrix::solve(self$Sigma)
                      } else {
                        invS <- 1
                      }
                      
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
                      
                    },
                    gen_stan_data = function(type="s",
                                             y = NULL){
                      B <- length(self$covariance$.__enclos_env__$private$flist)
                      N_dim <- c(unlist(rev(self$covariance$.__enclos_env__$private$flistcount)),0)
                      N_func <-  c(unlist(lapply(self$covariance$.__enclos_env__$private$Funclist,function(x)ncol(x))),0)
                      func_def <- matrix(0,nrow=B,ncol=max(N_func))
                      for(b in 1:B)func_def[b,1:N_func[b]] <- self$covariance$.__enclos_env__$private$Funclist[[b]][1,]
                      fvar <- lapply(rev(self$covariance$.__enclos_env__$private$flistvars),function(x)x$groups)
                      nvar <- lapply(lapply(fvar,table),as.vector)
                      N_var_func <- matrix(0,ncol=max(N_func),nrow=B)
                      for(b in 1:B)N_var_func[b,1:N_func[b]] <- nvar[[b]]
                      col_id <- array(0,dim=c(B,max(N_func),max(N_var_func)))
                      
                      for(b in 1:B){
                        for(k in 1:N_func[b]){
                          vnames <- rev(self$covariance$.__enclos_env__$private$flistvars)[[b]]
                          vnames <- vnames$rhs[vnames$groups==(N_func[b]+1-k)]
                          col_id[b,k,1:N_var_func[b,k]] <- match(vnames,colnames(self$covariance$.__enclos_env__$private$Distlist[[b]]))
                        }
                      }
                      n_par <- matrix(0,nrow=B,ncol=max(N_func))
                      for(b in 1:B)for(k in 1:N_func[b])n_par[b,k] <- sum(self$covariance$.__enclos_env__$private$Funclist[[b]][6:9,k] > -1)
                      cov_data <- array(0,dim=c(B,max(N_dim),max(rowSums(N_var_func))))
                      for(b in 1:B)cov_data[b,1:N_dim[b],1:ncol(self$covariance$.__enclos_env__$private$Distlist[[b]])] <- self$covariance$.__enclos_env__$private$Distlist[[b]]
                      
                      stan_data <- list(
                        B = B,
                        N_dim = N_dim,
                        max_N_dim = max(N_dim),
                        N_func = N_func,
                        max_N_func = max(N_func),
                        func_def = func_def,
                        N_var_func = N_var_func,
                        max_N_var = max(rowSums(N_var_func)),
                        max_N_var_func = max(N_var_func),
                        col_id = col_id,
                        N_par = n_par,
                        sum_N_par = sum(n_par),
                        cov_data = cov_data,
                        N = self$n(),
                        P = ncol(self$mean_function$X),
                        Q = ncol(self$covariance$Z),
                        X = as.matrix(self$mean_function$X),
                        Z = as.matrix(self$covariance$Z)
                      )
                      if(type=="m"){
                        for(i in 1:length(stan_data))names(stan_data)[i] <- paste0(names(stan_data)[i],"_m")
                      } else if(type=="y"){
                        if(is.null(y))stop("Must specify y for type y stan data")
                        stan_data$y <- y
                      }
                      return(stan_data)
                    }
                  ))


