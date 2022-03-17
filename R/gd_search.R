check_psd <- function(M){
  cholM <- tryCatch(chol(M),error=function(e)NA)
  if(!is(cholM,"matrix"))
  {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

## COMMENTING OUT ALL THE BELOW AS UNUSED


# optim_fun <- function(C,X,S){
#   if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
#   M <- t(X) %*% solve(S) %*% X
#   val <- t(C) %*% solve(M) %*% C
#   return(val)
# }
# 
# optim_fun2 <- function(C,X,S){
#   if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
#   M <- t(X) %*% solve(S) %*% X
#   val <- diag(solve(M))[c(C) != 0]
#   return(val)
# }
# 
# M_fun <- function(C,X,S){
#   M <- t(X) %*% solve(S) %*% X
#   return(M)
# }
# 
# hill_climb <-function(idx_in, 
#                       C_list, 
#                       X_list, 
#                       sig_list, 
#                       exp_cond, 
#                       w=NULL, 
#                       tol=1e-9,
#                       trace = TRUE, 
#                       rm_cols = NULL, 
#                       nfix = 0, 
#                       rd_mode = 1){
#   
#   if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
#   if(sum(w)!=1)w <- w/sum(w)
#   if(!is(w,"matrix"))w <- matrix(w,ncol=1)
#   if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
#   if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
#   if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
#   if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")
# 
#   #define sample size
#   N <- nrow(X_list[[1]])
# 
#   # added this function to give the user the option to remove columns from particular
#   # designs quickly if the algorithm previously stopped and said to remove
#   if(!is.null(rm_cols))
#   {
#     if(!is(rm_cols,"list"))stop("rm_cols should be a list")
#     idx_original <- list()
#     zero_idx <- c()
#     idx_original <- 1:nrow(X_list[[1]])
# 
#     # find all the entries with non-zero values of the given columns in each design
#     for(i in 1:length(rm_cols))
#     {
#       if(!is.null(rm_cols[[i]])){
#         for(j in 1:length(rm_cols[[i]]))
#         {
#           zero_idx <- c(zero_idx,which(X_list[[i]][,rm_cols[[i]][j]]!=0))
#         }
#       }
#     }
#     zero_idx <- sort(unique(zero_idx))
#     idx_original <- idx_original[-zero_idx]
#     idx_in <- match(idx_in,idx_original)
# 
#     if(trace)message(paste0("removing ",length(zero_idx)," observations"))
# 
#     #update the matrices
#     for(i in 1:length(rm_cols))
#     {
#       X_list[[i]] <- X_list[[i]][-zero_idx,-rm_cols[[i]]]
#       C_list[[i]] <- matrix(C_list[[i]][-rm_cols[[i]]],ncol=1)
#       sig_list[[i]] <- sig_list[[i]][-zero_idx,-zero_idx]
#     }
# 
#     if(any(is.na(idx_in)))
#     {
#       if(trace)message("generating new random starting point")
#       idx_in <- sample(1:nrow(X_list[[1]]),n,replace=FALSE)
#     }
#   }
# 
#   #MAIN BODY OF THE FUNCTION
#   
#   out_list <- GradRobustStep(N,
#                              idx_in, 
#                              C_list, 
#                              X_list, 
#                              sig_list,
#                              exp_cond = exp_cond,
#                              nfix = 0,
#                              weights = c(1), 
#                              rd_mode=1)
#   
#   out_list <- GradRobustStep(idx_in -1, do.call(rbind, C_list), do.call(rbind, X_list), do.call(rbind, sig_list), weights = w, nfix, rd_mode)
#   idx_in <- out_list[["idx_in"]] + 1
#   idx_out <- out_list[["idx_out"]] + 1
#   best_val_vec <- out_list[["best_val_vec"]]
# 
#   # as a check, see if the next swap would mean that M is not positive semi-definite, in which
#   # case the algorithm has terminated early and the solution is not in this design
#   # this code is copied from the choose_swap_robust function, so could be streamlined into its own function perhaps
# 
#   idx_test <- c(idx_in[2:n],idx_out[1])
#   val_out_mat <- matrix(NA,nrow=length(idx_test),ncol=length(A_list))
#   for(idx in 1:length(A_list)){
#     val_out_mat[,idx] <- sapply(1:length(idx_test),function(i)
#       remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_test]))
#   }
#   val_out <- as.numeric(val_out_mat %*% w)
#   idx_test <- idx_test[-which.max(val_out)]
#   rm1A <- list()
#   for(idx in 1:length(A_list)){
#     rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
#   }
#   for(j in 1:length(sig_list))
#   {
#     M_list[[j]] <- gen_m(X_list[[j]][idx_test,],rm1A[[j]])
# 
#     if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
#                                            " of design ",j," is not part of an optimal design."))
#   }
# 
#   ## if columns were removed then return the index to the indexing of the original X matrix
#   if(!is.null(rm_cols))
#   {
#     idx_in <- idx_original[idx_in]
#   }
# 
#   #return variance
#   if(trace){
#     cat("\nVariance for individual model(s):\n")
#     print(c(best_val_vec))
#     if(length(A_list)>1){
#       cat("\n Weighted average variance:\n")
#       print(sum(best_val_vec*c(w)))
#     }
#   }
# 
#   return(idx_in)
# }

## HAS THERE BEEN ANY DEVELOPMENT ON THE BELOW?

# For a given m find the optimal power vector
max_var <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)

  if (length(idx_in) != nrow(X_list[[1]]))
    idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  v0
}



# For a given m find the optimal power vector
max_power <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)
  #idx_in <- (1:m)*round(nrow(X_list[[1]])/m,0)

  if (length(idx_in) != nrow(X_list[[1]]))
    idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  pow <- pnorm(sqrt(theta[unlist(C_list)!=0]/sqrt(v0)) - qnorm(1-alpha/2))

  pow
}

sample_size <- function(theta, alpha, pwr_target, m, C_list, X_list, sig_list, w) {
  iter <- 0
  pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w)
  while (!all(pwr_new - pwr_target > 0) & m < nrow(X_list[[1]])) {
    iter <- iter + 1
    m <- m + 1
    pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE)
    cat("\nm = ", m)
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }
  return(m)
}

sample_size2 <- function(theta, alpha, pwr_target, C_list, X_list, sig_list, w) {

  cat("\nTarget power = ", pwr_target)

  lo <- max(unlist(lapply(C_list,function(i)length(unlist(i)))))*3
  hi <- nrow(X_list[[1]])
  pwr_new_lo <-NULL
  while(is.null(pwr_new_lo)){
    cat("\nlo = ", lo)
    pwr_new_lo <- tryCatch(
      max_power(theta, alpha, lo, C_list, X_list, sig_list, w, trace = FALSE),
      error=function(i)NULL)
    lo <- lo+10
  }
  pwr_new_hi <- max_power(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\nmin power = ", min(pwr_new_lo))
  cat("\nmax power = ", min(pwr_new_hi))

  if (min(pwr_new_hi) < pwr_target | min(pwr_new_lo) > pwr_target)
    stop("\ntarget power is not in range of ", min(pwr_new_lo) , " and ", min(pwr_new_hi))

  v_hi    <- max_var(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)
  v_target<- (theta[unlist(C_list)!=0]/(qnorm(pwr_target) + qnorm(1-alpha/2))^2)^2
  guess   <- round(max(v_hi / v_target * hi))
  pwr_new_guess <- max_power(theta, alpha, guess, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\ninitial guess = ", guess, " with power = ", min(pwr_new_guess))

  if (min(pwr_new_guess) < pwr_target) lo <- guess
  if (min(pwr_new_guess) > pwr_target) hi <- guess

  while (lo <= hi) {
    mid <- lo + round((hi - lo) / 2)
    cat("\nlo = ", lo)
    cat("  hi = ", hi)
    cat(" mid = ", mid)
    pwr_new <- max_power(theta, alpha, mid, C_list, X_list, sig_list, w, trace = FALSE)
    if (pwr_target < min(pwr_new)) hi = mid - 1
    if (pwr_target > min(pwr_new)) lo = mid + 1
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }

  return(mid)
}
