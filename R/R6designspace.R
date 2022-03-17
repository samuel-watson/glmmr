# DESIGN SPACE CLASS

DesignSpace <- R6::R6Class("DesignSpace",
                 public = list(
                   weights = NULL,
                   experimental_condition = NULL,
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
                   add = function(x) {
                     private$designs <- append(private$designs, list(x))
                     self$weights <- rep(1/length(private$designs),length(private$designs))
                     invisible(self)
                   },
                   remove = function(index) {
                     if (private$length() == 0) return(NULL)
                     private$designs <- private$designs[-index]
                   },
                   print = function(){
                     cat(paste0("Design space with ",self$n()," design(s): \n"))
                     for(i in 1:length(private$designs)){
                       cat(paste0("=========================================================\nDESIGN ",i,"(weight ",self$weights[i],"):\n"))
                       print(private$designs[[i]])
                     }
                   },
                   power = function(par,
                                    value,
                                    alpha=0.05){
                     
                     #change/remove this function
                     pwr <- c()
                     if(self$n()>1 && length(par)==1)par <- rep(par,self$n())
                     if(self$n()>1 && length(value)==1)value <- rep(value,self$n())
                     for(i in 1:self$n()){
                       pwr[i] <- private$designs[[i]]$power(par[i],value[i],alpha)
                     }
                     return(data.frame(Design = 1:self$n(),power = pwr))
                   },
                   n = function(){
                     length(private$designs)
                   },
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
                     for(i in 1:self$n()){
                       for(j in unique(self$experimental_condition)){
                         uncorr <- all(private$designs[[i]]$Sigma[which(self$experimental_condition==j),which(self$experimental_condition!=j)]==0)
                         if(!uncorr)break
                       }
                       if(!uncorr)break
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
                       idx_out <- drop(out_list[["idx_in"]] )
                     }
                     
                     
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
                   },
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
                     order(weights,decreasing = TRUE)[1:m]
                   }
                 )
)








