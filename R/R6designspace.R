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
                                      verbose=TRUE){
                     if(keep&verbose)message("linked design objects will be overwritten with the new design")
                     
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
                     
                     if(verbose&uncorr)message("Experimental conditions uncorrelated, using second-order cone program")
                     if(verbose&!uncorr)message("Experimental conditions correlated, using hill climbing algorithm")
                     
                     if(uncorr){
                       # this returns the experimental designs to keep
                       idx_out <- private$socp_roptimal(C,m)
                       
                     } else {
                       #initialise from random starting index
                       n <- private$designs[[1]]$mean_function$n()
                       idx_in <- sort(sample(1:n,m,replace=FALSE))
                       if(!is(C,"list")){
                         Clist <- list()
                         for(i in 1:self$n()){
                           Clist[[i]] <- matrix(C,ncol=1)
                         }
                       }
                       
                       ## this now expects that idx_out are the experimental conditions to keep and NOT the row numbers
                       ## see rows_to_keep for the rows to keep
                       
                       idx_out <- grad_robust2_step(
                         idx_in,
                         C_list = Clist,
                         X_list = private$genXlist(),
                         sig_list = private$genSlist(),
                         rm_cols = rm_cols,
                         trace=verbose
                       )
                       
                     }
                     
                     
                     idx_out <- sort(idx_out)
                     rows_to_keep <- which(self$experimental_condition %in% idx_out)
                     
                     
                     if(keep){
                       for(i in 1:self$n()){
                         private$designs[[i]]$subset_rows(rows_to_keep)
                         ncol <- 1:ncol(private$designs[[i]]$mean_function$X)
                         if(!is.null(rm_cols))private$designs[[i]]$subset_cols(ncol[-rm_cols[[i]]])
                       }
                     }
                     
                     if(verbose)cat("Experimental conditions in the optimal design: ", idx_out)
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
                       constr[[i]] <- t(X[[i]])%*%Matrix::t(Matrix::chol(Matrix::solve(sig[[i]])))%*%z[c(1:n + n*(i-1))] == C
                     }
                     
                     for(i in n_ex){
                       #build expression
                       cons_str <- "weights[1]*p_norm(z[which(exp_cond==i)])"
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








