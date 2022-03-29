#' A GLMM Design Space
#' 
#' A class-based representation of a "design space" that contains one or more \link{glmmr}[Design] objects.
#' @details 
#' An experimental study is comprised of a collection of experimental conditions, which are one or more observations made a pre-specified locations/values
#' of covariates. A design space represents the collection of all possible experimental conditions for the study design and plausible models describing
#' the data generating process. The main purpose of this class is to identify optimal study designs, that is the set of `n` experimental conditions 
#' from all possible experimental conditions that minimise the variance of a parameter of interest across the specified GLMMs.
#' 
#' A DesignSpace object is intialised using one or more Design objects. Design objects can be added or removed from the collection. 
#' All designs must have the same number of rows in their design matrices (X and Z) and the same number of experimental conditions. 
#' The DesignSpace functions can modify the linked design objects.
DesignSpace <- R6::R6Class("DesignSpace",
                 public = list(
                   #' @field 
                   #' A vector denoting the prior weighting of each Design in the design space. Required if robust optimisation is used based on a 
                   #' weighted average variance over the linked designs. If it is not specified in the call to `new()` then designs are assumed
                   #' to have equal weighting.
                   weights = NULL,
                   #' @field 
                   #' A vector indicating the unique identifier of the experimental condition for each observation/row in the matrices X and Z.
                   experimental_condition = NULL,
                   #' @description 
                   #' Create a new Design Space
                   #' 
                   #' Creates a new design space from one or more glmmr designs.
                   #' @details 
                   #' The experimental condition refers to the smallest "unit" of the study design that could be included in the design. For example, in a
                   #' cluster randomised trial, the experimental condition may be single individuals such that we can observed any number of individuals 
                   #' in any cluster period (including none at all). In this case the experimental condition would be equivalent to row number. Alternatively,
                   #' we may have to observe whole cluster periods, and we need to choose which cluster periods to observe, in which case the each observation 
                   #' in a different cluster-period would have the same experimental condition identifier. Finally, we may determine that the whole cluster in 
                   #' all periods (a "sequence") is either observed or not.
                   #' @param ... One or more glmmr \link{glmmr}[Design] objects. The designs must have an equal number of observations.
                   #' @param weights Optional. A numeric vector of values between 0 and 1 indicating the prior weights to assign to each of the designs. The weights
                   #' are required for optimisation, if a weighted average variance is used across the designs. If not specified then designs are assumed 
                   #' to have equal weighting.
                   #' @param experimental_condition Optional. A vector of the same length as the number of observations in each design indicating the unique
                   #' identifier of the experimental condition that observation belongs to, see Details. If not provided, then it is assumed that all observations
                   #' are separate experimental conditions.
                   #' @return A `DesignSpace` object
                   #' @examples
                   #' ...
                   initialize = function(...,
                                         weights=NULL,
                                         experimental_condition = NULL) {
                     samp.size <- c()
                     i <- 0
                     for (item in list(...)) {
                       i <- i + 1
                       if(!is(item,"Design"))stop(paste0(item," is not a Design"))
                       samp.size[i] <- item$n()
                     }
                     #print(samp.size)
                     if(length(samp.size) > 1 && !all.equal(samp.size))stop("designs are of different sizes")
                     samp.size <- unique(samp.size)
                     
                     #if the weights are null assign equal weighting
                     if(!is.null(weights)){
                       if(length(weights)!=length(list(...)))stop("weights not same length as designs")
                       self$weights <- weights
                     } else {
                       message("weights not provided, assigning equal weighting. weights can be changed manually in self$weights")
                       self$weights <- rep(1/length(list(...)),length(list(...)))
                     }
                     
                     #if the experimental condition is null assign all separate
                     if(!is.null(experimental_condition)){
                       if(length(experimental_condition)!=samp.size)stop("experimental condition not same size as designs")
                       self$experimental_condition <- experimental_condition
                     } else {
                       message("experimental condition not provided, assuming each observation is a separate experimental condition. experimental condition can be changed manually in self$experimental_condition")
                       self$experimental_condition <- 1:samp.size
                     }
                     
                     for (item in list(...)) {
                       self$add(item)
                     }

                   },
                   #' @description 
                   #' Add a design to the design space
                   #' @param x A `Design` to add to the design space
                   #' @return Nothing
                   #' @examples 
                   #' ...
                   add = function(x) {
                     if(x$n()!=private$designs[[1]]$n())stop("New design is not same size as designs in this design space.")
                     private$designs <- append(private$designs, list(x))
                     self$weights <- rep(1/length(private$designs),length(private$designs))
                     invisible(self)
                   },
                   #' @description 
                   #' Removes a design from the design space
                   #' @param index Index of the design to remove
                   #' @return Nothing
                   #' @examples 
                   #' ...
                   remove = function(index) {
                     if (private$length() == 0) return(NULL)
                     private$designs <- private$designs[-index]
                   },
                   #' @description 
                   #' Print method for the Design Space
                   #' @param ... ignored
                   #' @return Prints to the console all the designs in the design space
                   #' @examples 
                   #' ...
                   print = function(){
                     cat(paste0("Design space with ",self$n()," design(s): \n"))
                     for(i in 1:length(private$designs)){
                       cat(paste0("=========================================================\nDESIGN ",i,"(weight ",self$weights[i],"):\n"))
                       print(private$designs[[i]])
                     }
                   },
                   # power = function(par,
                   #                  value,
                   #                  alpha=0.05){
                   #   
                   #   #change/remove this function
                   #   pwr <- c()
                   #   if(self$n()>1 && length(par)==1)par <- rep(par,self$n())
                   #   if(self$n()>1 && length(value)==1)value <- rep(value,self$n())
                   #   for(i in 1:self$n()){
                   #     pwr[i] <- private$designs[[i]]$power(par[i],value[i],alpha)
                   #   }
                   #   return(data.frame(Design = 1:self$n(),power = pwr))
                   # },
                   #' @description 
                   #' Returns the size of the design space and number of observations
                   #' @examples 
                   #' ...
                   n = function(){
                     c("n.designs"=length(private$designs),"n" = private$designs[[1]]$n())
                   },
                   #' @description 
                   #' Identify a c-optimal design of size m
                   #' 
                   #' Algorithms to identify a c-optimal design of size m within the design space.
                   #' @details 
                   #' The algorithm identifies a c-optimal design of size m from the design space with N designs each with n observations. The objective
                   #' function is
                   #' 
                   #' \deqn{C^TM^{-1}C}
                   #' 
                   #' where M is the information matrix and C is a vector. Typically C will be a vector of zeros with a single 1 in the position of the
                   #' parameter of interest. For example, if the columns of X in the design are an interept, the treatment indicator, and then time 
                   #' period indicators, the vector C may be `c(0,1,0,0,...)`, such that the objective function is the variance of that parameter. 
                   #' If there are multiple designs in the design space, the C vectors do 
                   #' not have to be the same as the columns of X in each design might differ, in which case a list of vectors can be provided.
                   #' 
                   #' If the experimental conditions are correlated with one another, then a hill climbing algorithm is used to find the optimal 
                   #' design by using the convexity of the objective function to "climb the hill" towards the optimal design. 
                   #' If the experimental conditional are uncorrelated (but there is correlation between observations within the same
                   #' experimental condition) then optionally a fast algorithm can be used to approximate the optimal design using a second-order 
                   #' cone program (see Sangol 2015 and van Dette). The approximate algorithm will return weights for each unique experimental condition representing
                   #' the "proportion of effort" to spend on each design condition. There are different ways to translate these weights into integer
                   #' values. Use of the approximate optimal design algorithm can be disabled used `force_hill=TRUE`
                   #' 
                   #' In some cases the optimal design will not be full rank with respect to the design matrix X of the design space. This will result
                   #' in a non-positive definite information matrix, and an error. The program will indicate which columns of X are likely "empty" in the optimal
                   #' design. The user can then optionally remove these columns in the algorithm using the `rm_cols` argument, which will delete the
                   #' specified columns and linked observations before starting the algorithm. 
                   #' 
                   #' The algorithm will also identify robust optimal designs if there are multiple designs in the design space. 
                   #' There are two options for robust optimisation. First, a weighted average of objective functions, where the weights are specified 
                   #' by the `weights` field in the design space (`robust_function = "weighted"`). The weights may represent the prior probability or plausibility of each design, 
                   #' for example. Second, a minimax approach can be used, where the function identifies the design that minimises the maximum objective
                   #' function across all designs (`robust_function = "minimax"`).
                   #' @param m A positive integer specifying the number of experimental conditions to include.
                   #' @param C Either a vector or a list of vectors of the same length as the number of designs, see Details.
                   #' @param rm_cols Optional. A list of vectors indicating columns of X to remove from each design, see Details.
                   #' @param keep Logical indicating whether to "keep" the optimal design in the linked design objects and remove any experimental
                   #' conditions and columns that are not part of the optimal design. Irreversible, so that these observations will be lost from the 
                   #' linked design objects. Defaults to FALSE.
                   #' @param verbose Logical indicating whether to reported detailed output on the progress of the algorithm. Default is TRUE.
                   #' @param robust_function Either "weighted" or "minimax". Specifies the objective function to use in robust optimisation, see Details.
                   #' @param force_hill Logical. If the experimental conditions are uncorrelated, if this option is TRUE then the hill climbing 
                   #' algorithm will be used, otherwise if it is FALSE, then a fast approximate alternative will be used. See Details
                   #' @return A vector indicating the identifiers of the experimental conditions in the optimal design, or a vector indicating the
                   #' weights if the approximate algorithm is used. Optionally the linked designs are also modified (see option `keep`).
                   #' @examples
                   #' ...
                   optimal = function(m,
                                      C,
                                      rm_cols=NULL,
                                      keep=FALSE,
                                      verbose=TRUE,
                                      robust_function = "weighted",
                                      force_hill=FALSE){
                     if(keep&verbose)message("linked design objects will be overwritten with the new design")
                     
                     ## add checks
                     
                     # dispatch to correct algorithm
                     # check if the experimental conditions are correlated or not
                     #loop through each sigma
                     if(verbose)message("Checking experimental condition correlations...")
                     if(length(self$experimental_condition)!=private$designs[[1]]$n())stop("experimental condition not the same length as design")
                     uncorr <- TRUE
                     unique_exp_cond <- unique(self$experimental_condition)
                     for(i in 1:self$n()){
                       for(j in unique_exp_cond){
                         uncorr <- all(private$designs[[i]]$Sigma[which(self$experimental_condition==j),which(self$experimental_condition!=j)]==0)
                         if(!uncorr)break
                       }
                       if(!uncorr)break
                     }
                     ## need to detect if the experimental conditions are duplicated
                     ## can update this but currently only provides a warning to the user
                     if(uncorr&!force_hill){
                       datahashes <- c()
                       for(j in unique_exp_cond){
                         datalist <- list()
                         for(k in 1:self$n()){
                           datalist[[k]] <- list(des$mean_function$X[self$experimental_condition==j,],
                                                 des$Sigma[self$experimental_condition==j,self$experimental_condition==j])
                         }
                         datahashes <- c(datahashes, digest::digest(datalist))
                       }
                       
                       if(any(duplicated(datahashes))){
                         unique_hash <- unique(datahashes)
                         n_unique_hash <- length(unique_hash)
                         datahashes <- match(datahashes,unique_hash)
                         message(paste0("Duplicated experimental conditions in the design space, ",n_unique_hash," unique 
experimental conditions, which are uncorrelated. 
force_hill=FALSE so weights will be calculated for each experimental condition separately. Sum of weights for
each condition will be reported below."))
                       }
                     }
                     
                     if(!is(C,"list")){
                       C_list <- list()
                       for(i in 1:self$n()){
                         C_list[[i]] <- matrix(C,ncol=1)
                       }
                     } else {
                       C_list <- C
                     }
                     
                     if(verbose&uncorr&!force_hill)message("Experimental conditions uncorrelated, using second-order cone program")
                     if(verbose&uncorr&force_hill)message("Experimental conditions uncorrelated, but using hill climbing algorithm")
                     if(verbose&!uncorr)message("Experimental conditions correlated, using hill climbing algorithm")
                     
                     if(uncorr&!force_hill){
                       # this returns the experimental designs to keep
                       idx_out <- private$socp_roptimal(C_list,m)
                       idx_out <- drop(idx_out)
                       if(verbose)cat("Weights for each experimental condition in the optimal design: ", idx_out)
                       if(any(duplicated(datahashes))){
                         agg_weights <- aggregate(idx_out,list(datahashes),sum)
                         cat("\nSum of weights for unique experimental conditions: ",agg_weights$x)
                         idx_out <- list(weights = idx_out, unique_weights = agg_weights$x)
                       }
                       return(invisible(idx_out))
                     } else {
                       #initialise from random starting index
                       N <- private$designs[[1]]$mean_function$n()
                       X_list <- private$genXlist()
                       sig_list <- private$genSlist()
                       weights <- self$weights
                       exp_cond <- as.numeric(as.factor(self$experimental_condition))
                       ncond <- length(unique(exp_cond))
                       idx_in <- sort(sample(1:ncond,m,replace=FALSE))
                       rdmode <- ifelse(robust_function=="weighted",1,0)
                       
                       if(!is.null(rm_cols))
                       {
                         if(!is(rm_cols,"list"))stop("rm_cols should be a list")
                         idx_original <- list()
                         zero_idx <- c()
                         idx_original <- 1:nrow(X_list[[1]])
                         
                         # find all the entries with non-zero values of the given columns in each design
                         for(i in 1:length(rm_cols))
                         {
                           if(!is.null(rm_cols[[i]])){
                             for(j in 1:length(rm_cols[[i]]))
                             {
                               zero_idx <- c(zero_idx,which(X_list[[i]][,rm_cols[[i]][j]]!=0))
                             }
                           }
                         }
                         zero_idx <- sort(unique(zero_idx))
                         idx_original <- idx_original[-zero_idx]
                         idx_in <- match(idx_in,idx_original)
                         
                         if(verbose)message(paste0("removing ",length(zero_idx)," observations"))
                         
                         #update the matrices
                         for(i in 1:length(rm_cols))
                         {
                           X_list[[i]] <- X_list[[i]][-zero_idx,-rm_cols[[i]]]
                           C_list[[i]] <- matrix(C_list[[i]][-rm_cols[[i]]],ncol=1)
                           sig_list[[i]] <- sig_list[[i]][-zero_idx,-zero_idx]
                         }
                       }
                       
                       out_list <- GradRobustStep(N = N,
                                                  idx_in = idx_in, 
                                                  C_list = C_list, 
                                                  X_list = X_list, 
                                                  sig_list = sig_list,
                                                  exp_cond = exp_cond,
                                                  nfix = 0,
                                                  weights = weights, 
                                                  rd_mode=rdmode,
                                                  trace=verbose)
                       
                       # old algorithm - it doesn't really work!
                       #   out_list <- GradRobustAlg1(N = N,
                       #                              idx_in = idx_in, 
                       #                              C_list = C_list, 
                       #                              X_list = X_list, 
                       #                              sig_list = sig_list,
                       #                              exp_cond = exp_cond,
                       #                              nfix = 0,
                       #                              weights = weights, 
                       #                              rd_mode=rdmode,
                       #                              trace=verbose)
                       
                       idx_out <- drop(out_list[["idx_in"]] )
                       idx_out <- sort(idx_out)
                       idx_out_exp <- unique(self$experimental_condition)[idx_out]
                       rows_to_keep <- which(self$experimental_condition %in% idx_out_exp)
                       
                       
                       if(keep){
                         for(i in 1:self$n()){
                           private$designs[[i]]$subset_rows(rows_to_keep)
                           ncol <- 1:ncol(private$designs[[i]]$mean_function$X)
                           if(!is.null(rm_cols))private$designs[[i]]$subset_cols(ncol[-rm_cols[[i]]])
                         }
                       }
                       if(verbose)cat("Experimental conditions in the optimal design: ", idx_out_exp)
                       return(invisible(idx_out_exp))
                     }
                     
                     
                     
                     
                     
                   },
                   #' @description 
                   #' Returns a linked design
                   #' @param i Index of the design to return
                   #' @examples 
                   show = function(i){
                     return(private$designs[[i]])
                   }
                 ),
                 private = list(
                   designs = list(),
                   genXlist = function(){
                     X_list <- list()
                     for(i in 1:self$n()){
                       X_list[[i]] <- as.matrix(private$designs[[i]]$mean_function$X)
                     }
                     return(X_list)
                   },
                   genSlist = function(){
                     S_list <- list()
                     for(i in 1:self$n()){
                       S_list[[i]] <- as.matrix(private$designs[[i]]$Sigma)
                     }
                     return(S_list)
                   },
                   socp_roptimal = function(C,
                                            m){
                     X <- private$genXlist()
                     sig <- private$genSlist()
                     weights <- self$weights
                     exp_cond <- self$experimental_condition
                     
                     n_r <- length(sig)
                     constr <- list()
                     n <- nrow(X[[1]])
                     n_ex <- unique(exp_cond)
                     mu <- CVXR::Variable(length(n_ex))
                     z <- CVXR::Variable(n*n_r)
                     
                     for(i in 1:n_r){
                       constr[[i]] <- t(X[[i]])%*%Matrix::t(Matrix::chol(Matrix::solve(sig[[i]])))%*%z[c(1:n + n*(i-1))] == C[[i]]
                     }
                     
                     for(i in n_ex){
                       #build expression
                       cons_str <- "weights[1]*CVXR::p_norm(z[which(exp_cond==i)])"
                       if(n_r > 1){
                         for(j in 1:(n_r-1)){
                           cons_str <- paste0(cons_str," + weights[",j+1,"]*p_norm(z[(which(exp_cond==i)+",j,"*n)])")
                         }
                       }
                       cons_str <- paste0(cons_str, " <= mu[i]")
                       
                       constr[[(length(constr)+1)]] <- eval(str2lang(cons_str))
                     }
                     obj <- sum(mu)
                     prob <- CVXR::Problem(CVXR::Minimize(obj),constr)
                     stopifnot(CVXR::is_dcp(prob))
                     res <- CVXR::solve(prob)
                     weights <- res$getValue(mu)
                     # choose the m biggest to keep
                     weights/sum(weights)
                     #order(weights,decreasing = TRUE)[1:m]
                   }
                 )
)








