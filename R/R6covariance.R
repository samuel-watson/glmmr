#' R6 Class representing a covariance function and data
#' 
#' For the generalised linear mixed model 
#' 
#' \deqn{Y \sim F(\mu,\sigma)}
#' \deqn{\mu = h^-1(X\beta + Z\gamma)}
#' \deqn{\gamma \sim MVN(0,D)}
#' 
#' where h is the link function, this class defines Z and D. The covariance is defined by a covariance function, data, and parameters.
#' A new instance can be generated with $new(). The class will generate the 
#' relevant matrices Z and D automatically.  
#' @export
Covariance <- R6::R6Class("Covariance",
                      public = list(
                        #' @field data Data frame with data required to build covariance
                        data=NULL,
                        #' @field formula Covariance function formula. See `help(Covariance$new())` for details.
                        formula = NULL,
                        #' @field parameters List of lists holding the model parameters. See `help(Covariance$new())` for details.
                        parameters = NULL,
                        #' @field Z Design matrix
                        Z = NULL,
                        #' @field D Covariance matrix of the random effects
                        D = NULL,
                        #' @description 
                        #' Return the size of the design
                        #' @return Scalar 
                        n= function(){
                          nrow(self$Z)
                        },
                        #' @description 
                        #' Create a new Covariance object
                        #' @param formula Formula describing the covariance function. See Details
                        #' @param data Data frame with data required for constructing the covariance.
                        #' @param parameters List of lists with parameter values for the functions in the model
                        #' formula. See Details.
                        #' @param verbose Logical whether to provide detailed output.
                        #' @details A covariance function is specified as an additive formula made up of 
                        #' components with structure \code{(1|f(j))}. The left side of the vertical bar 
                        #' specifies the covariates in the model that have a random effects structure. 
                        #' The right side of the vertical bar specify the covariance function `f` for 
                        #' that term using variable named in the data `j`. If there are multiple 
                        #' covariates on the left side, it is assumed their random effects are 
                        #' correlated, e.g. \code{(1+x|f(j))}. Additive functions are assumed to be 
                        #' independent, for example, \code{(1|f(j))+(x|f(j))} would create random effects 
                        #' with zero correlation for the intercept and the parameter on covariate \code{x}. 
                        #' Covariance functions on the right side of the vertical bar are multiplied 
                        #' together, i.e. \code{(1|f(j)*g(t))}. 
                        #' 
                        #' There are several common functions included for a named variable in data \code{x},
                        #' see respective help files for the functions for parameterisation:
                        #' * \code{gr(x)}: Indicator function   
                        #' * \code{fexp(x)}: Exponential function
                        #' * \code{pexp(x)}: Power function
                        #' 
                        #' One can add other functions by specifying a function that takes a list as an 
                        #' argument with first element called data that contains the data, and a second 
                        #' element called pars that contains the parameters as a vector.
                        #' 
                        #' Parameters are provided to the covariance function as a list of lists. 
                        #' The elements of the list should be lists corresponding to the additive elements of the 
                        #' covariance function. Each of those lists should have elements that are vectors or scalars 
                        #' providing the values of the parameters for each function in the order they are written. 
                        #' For example,
                        #' * Formula: `~(1|gr(j))+(1|gr(j)*gr(t))`; parameters: `c(0.25,0.1)`
                        #' * Formula: `~(1|gr(j)*fexp(t))`; parameters: `c(0.25,0.5)`
                        #' @return A Covariance object
                        #' @seealso \code{gr}, \code{fexp}, \code{pexp}
                        #' @examples 
                        #' df <- nelder(~(cl(5)*t(5)) > ind(5))
                        #' cov <- Covariance$new(formula = ~(1|gr(j)*pexp(t)),
                        #'                       parameters = c(0.25,0.7),
                        #'                       data= df)
                        initialize = function(formula=NULL,
                                              data = NULL,
                                              parameters= NULL,
                                              verbose=TRUE){
                          if(any(is.null(data),is.null(formula),is.null(parameters))){
                            cat("not all attributes set. call check() when all attributes set.")
                          } else {
                            self$data = data
                            self$formula = as.formula(formula, env=.GlobalEnv)
                            self$parameters = parameters

                            private$cov_form()
                          }
                        },
                        #' @description 
                        #' Check if anything has changed and update matrices if so.
                        #' @param verbose Logical whether to report if any changes detected.
                        #' @return NULL
                        #' @examples 
                        #' df <- nelder(~(cl(5)*t(5)) > ind(5))
                        #' cov <- Covariance$new(formula = ~(1|gr(j)*pexp(t)),
                        #'                       parameters = list(list(0.05,0.8)),
                        #'                       data= df)
                        #' cov$parameters <- list(list(0.01,0.1))
                        #' cov$check(verbose=FALSE)
                        check = function(verbose=TRUE){
                          new_hash <- private$hash_do()
                          if(private$hash[1] != new_hash[1]){
                            if(verbose)message("changes found, updating Z")
                            private$cov_form()
                          } else if(private$hash[2] != new_hash[2]){
                            if(verbose)message("changes found, updating D")
                            private$genD()
                          }

                          invisible(self)
                        },
                        #' @description 
                        #' Show details of Covariance object
                        #' @param ... ignored
                        #' @examples
                        #' df <- nelder(~(cl(5)*t(5)) > ind(5))
                        #' Covariance$new(formula = ~(1|gr(j)*pexp(t)),
                        #'                       parameters = list(list(0.05,0.8)),
                        #'                       data= df)
                        print = function(){
                          # MAKE CLEARER ABOUT FUNCTIONS AND PARAMETERS?

                          cat("Covariance\n")
                          print(self$formula)
                          cat("Parameters:")
                          print(unlist(self$parameters))
                          #print(head(self$data))
                        },
                        #' @description 
                        #' Keep specified indices and removes the rest
                        #' @param index vector of indices to keep
                        #' @examples 
                        #' df <- nelder(~(cl(10)*t(5)) > ind(10))
                        #' cov <- Covariance$new(formula = ~(1|gr(j)*pexp(t)),
                        #'                       parameters = list(list(0.05,0.8)),
                        #'                       data= df)
                        #' cov$subset(1:100)                     
                        subset = function(index){
                          self$data <- self$data[index,]
                          self$check()
                        },
                        #' @description 
                        #' Generate a new D matrix
                        #' 
                        #' @details 
                        #' D is the covariance matrix of the random effects terms in the generalised linear mixed
                        #' model. This function will return a matrix D for a given set of parameters.
                        #' @param parameters list of lists, see initialize()
                        #' @return matrix 
                        #' @examples 
                        #' df <- nelder(~(cl(10)*t(5)) > ind(10))
                        #' cov <- Covariance$new(formula = ~(1|gr(j)*pexp(t)),
                        #'                       parameters = list(list(0.05,0.8)),
                        #'                       data= df)
                        #' cov$sampleD(list(list(0.01,0.1)))
                        sampleD = function(parameters){
                          return(private$genD(update = FALSE,
                                              new_pars = parameters))
                        }
                      ),
                      private = list(
                        details = list(), # update to contain details of covariance for printing later
                        Distlist = NULL,
                        flist = NULL,
                        flistvars = NULL,
                        Funclist = NULL,
                        Zlist = NULL,
                        cov_functions = function(arg,x,pars){
                          if(arg=="exponential"){
                            return(pars[1]*exp(-x/pars[2]))
                          }
                          if(arg == "indicator"){
                            return(pars[1]*I(x==0))
                          }
                          if(arg == "exp_power"){
                            return(pars[1]^(x))
                          }
                        },
                        hash = NULL,
                        hash_do = function(){
                          c(digest::digest(c(self$data)),digest::digest(c(self$formula,self$parameters)))
                        },
                        cov_form = function(){
                          #1. split into independent components that can be combined in block diagonal form
                          flist <- list()
                          flistvars <- list()
                          formExt <- TRUE
                          count <- 0
                          form0 <- self$formula[[2]]
                          while(formExt){
                            count <- count + 1
                            formExt <- I(form0[[1]]=="+")
                            if(formExt){
                              flist[[count]] <- form0[[3]][[2]]
                              form0 <- form0[[2]]
                            } else {
                              flist[[count]] <- form0[[2]]
                            }
                            if("+"%in%all.names(flist[[count]][[3]]))stop("covariance only products")
                            rhsvars <- all.vars(flist[[count]][[3]])
                            funlist <- list()
                            rhsvargroup <- rep(NA,length(rhsvars))
                            formMult <- TRUE
                            countMult <- 0
                            form3 <- flist[[count]][[3]]
                            while(formMult){
                              countMult <- countMult + 1
                              formMult <- I(form3[[1]] == "*")
                              if(formMult){
                                funlist[[countMult]] <- form3[[3]][[1]]
                                rhsvargroup[match(all.vars(form3[[3]]),rhsvars)] <- countMult
                                form3 <- form3[[2]]
                              } else {
                                funlist[[countMult]] <- form3[[1]]
                                rhsvargroup[match(all.vars(form3),rhsvars)] <- countMult
                              }

                            }
                            flistvars[[count]] <- list(lhs=all.vars(flist[[count]][[2]]),
                                                       rhs = rhsvars,
                                                       funs = funlist,
                                                       groups = rhsvargroup)

                          }

                          # build each Z matrix and cbind
                          Zlist <- list()
                          Distlist <- list()
                          Funclist <- list()
                          for(i in 1:length(flist)){
                            data_nodup <- self$data[!duplicated(self$data[,flistvars[[i]]$rhs]),flistvars[[i]]$rhs]
                            if(!is(data_nodup,"data.frame")){
                              data_nodup <- data.frame(data_nodup)
                              colnames(data_nodup) <- flistvars[[i]]$rhs
                            }
                            Distlist[[i]] <- as.matrix(data_nodup)
                            zdim2 <- nrow(data_nodup)
                            Zlist[[i]] <- match_rows(self$data,data_nodup,by=flistvars[[i]]$rhs)
                            if(length(flistvars[[i]]$lhs)>0){
                              ZlistNew <- list()
                              for(j in 1:length(flistvars[[i]]$lhs)){
                                ZlistNew[[j]] <- Zlist[[i]]*df[,flistvars[[i]]$lhs[j]]
                              }
                              Zlist[[i]] <- Reduce(cbind,ZlistNew)
                            }
                            # Dist1 <- list()
                            # for(j in 1:length(unique(flistvars[[i]]$groups))){
                            #   Dist1[[j]] <- as.matrix(dist(data_nodup[,flistvars[[i]]$rhs[flistvars[[i]]$groups==j]],upper = TRUE, diag=TRUE))
                            # }
                            # Distlist[[i]] <- Dist1
                            
                            #create matrix for functions
                            
                          }
                          
                          fl <- rev(flistvars)
                          fnames <- c("gr","fexp","pexp")
                          fnpar <- c(1,2,1)
                          parcount <- 0
                          Funclist <- list()
                          
                          for(i in 1:length(fl)){
                            nfunc <- length(fl[[i]]$funs)
                            f1 <- matrix(-1,nrow=9,ncol=nfunc)
                            for(k in 1:nfunc){
                              f1[1,k] <- which(fnames==rev(fl[[i]]$funs)[[k]])
                              for(l in 1:sum(fl[[i]]$groups==(nfunc +1 - k))){
                                f1[1+l,k] <- which(fl[[i]]$groups==(nfunc +1 - k))[l] - 1
                              }
                              # deal with the parameters
                              npar1 <- fnpar[f1[1,k]]
                              f1[6:(5+npar1),k] <- ((parcount + 1):(parcount + npar1)) - 1
                              parcount <- parcount + npar1
                            }
                            
                            Funclist[[i]] <- f1
                          }
                          
                          Z <- Reduce(cbind,Zlist)
                          Z <- Matrix::Matrix(Z)
                          #if(ncol(Z)>nrow(Z))warning("Model underidentified")
                          if(nrow(Z) < ncol(Z))warning("More random effects than observations")
                          self$Z <- Z
                          private$Distlist <- rev(Distlist)
                          private$flistvars <- flistvars
                          private$flist <- flist
                          for(i in 1:length(Zlist))Zlist[[i]] <- Matrix::Matrix(Zlist[[i]])
                          private$Zlist <- Zlist
                          private$Funclist <- Funclist
                          
                          private$genD()

                        },
                        genD = function(update=TRUE,
                                        new_pars=NULL){
                          
                          D <- genD(private$Funclist,
                                    private$Distlist,
                                    self$parameters)
                          
                          # D.sublist <- list()
                          # for(d in 1:length(private$flist)){
                          #   ngroup <- length(unique(private$flistvars[[d]]$groups))
                          #   D.sublist[[d]] <- matrix(1, nrow=ncol(private$Zlist[[d]]),
                          #                            ncol=ncol(private$Zlist[[d]]))
                          #   if(update){
                          #      pars <-rev(self$parameters[[d]])[1:ngroup]
                          #   } else {
                          #     pars <- rev(new_pars[[d]])
                          #   }
                          #   for(j in 1:ngroup){
                          #     D.sublist[[d]] <- D.sublist[[d]]*
                          #       do.call(as.character(private$flistvars[[d]]$funs[[j]]),
                          #               list(list(data=private$Distlist[[d]][[j]],
                          #                         pars=pars[[j]])))
                          #   }
                          #   #remove AsIs class
                          #   class(D.sublist[[d]]) <- "matrix"
                          # 
                          #   if(length(private$flistvars[[d]]$lhs)>0){
                          #     #fix cov parameter here
                          #     Dmatlist <- list(D.sublist[[d]])
                          #     for(j in 2:length(private$flistvars[[d]]$lhs)){
                          #       Dmatlist[[j]] <- diag(zdim2)*self$parameters[[d]][[ngroup-1+j]]
                          #     }
                          #     D.sublist[[d]] <- do.call("blockmat",Dmatlist)
                          #   }
                          # }
                          # 
                          # #finally bdiag to combine them all
                          # D <- do.call(Matrix::bdiag,D.sublist)

                          #add warning if number of re > n
                          

                          if(update){
                            self$D <- Matrix::Matrix(D)
                            private$hash <- private$hash_do()
                          } else {
                            return(D)
                          }

                        }
                      ))

