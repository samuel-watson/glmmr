#' Prints an mcml fit output
#' 
#' Print method for class "`mcml`"
#' 
#' @param x an object of class "`mcml`" as a result of a call to MCML, see \link[glmmr]{Design}
#' @param digits Number of digits to print
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.mcml` tries to replicate the output of other regression functions, such
#' as `lm` and `lmer` reporting parameters, standard errors, and z- and p- statistics.
#' The z- and p- statistics should be interpreted cautiously however, as generalised
#' linear mixed models can suffer from severe small sample biases where the effective
#' sample size relates more to the higher levels of clustering than individual observations.
#' 
#' Parameters `b` are the mean function beta parameters, parameters `cov` are the
#' covariance function parameters in the same order as `$covariance$parameters`, and
#' parameters `d` are the estimated random effects.
#' @return TBC
#' @export
print.mcml <- function(x, digits =2, ...){
  cat("Markov chain Monte Carlo Maximum Likelihood Estimation\nAlgorithm: ",
      ifelse(x$method=="mcem","Markov Chain Expectation Maximisation",
             "Markov Chain Newton-Raphson"),
      ifelse(x$sim_step," with simulated likelihood step\n","\n"))
  
  cat("\nNumber of Monte Carlo simulations per iteration: ",x$m," with tolerance ",x$tol,"\n")
  semethod <- ifelse(x$permutation,"permutation test",ifelse(x$hessian,"hessian","approx"))
  cat("P-value and confidence interval method: ",semethod,"\n")
  pars <- x$coefficients[!grepl("d",x$coefficients$par),c('est','SE','lower','upper')]
  z <- pars$est/pars$SE
  pars <- cbind(pars[,1:2],z=z,p=2*(1-pnorm(abs(z))),pars[,3:4])
  colnames(pars) <- c("Estimate","Std. Err.","z value","p value","2.5% CI","97.5% CI")
  rownames(pars) <- x$coefficients$par[!grepl("d",x$coefficients$par)]
  pars <- apply(pars,2,round,digits = digits)
  print(pars)
  
  #messages
  if(x$permutation)message("Permutation test used for one parameter, other SEs are not reported. SEs and Z values
are approximate based on the p-value, and assume normality.")
  if(!x$hessian&!x$permutation)warning("Hessian was not positive definite, standard errors are approximate")
  if(!x$converged)warning("Algorithm did not converge")
  return(invisible(pars))
}

#' Summarises an mcml fit output
#' 
#' Summary method for class "`mcml`"
#' 
#' @param x an object of class "`mcml`" as a result of a call to MCML, see \link[glmmr]{Design}
#' @param digits Number of digits to print
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.mcml` tries to replicate the output of other regression functions, such
#' as `lm` and `lmer` reporting parameters, standard errors, and z- and p- statistics.
#' The z- and p- statistics should be interpreted cautiously however, as generalised
#' linear mixed models can suffer from severe small sample biases where the effective
#' sample size relates more to the higher levels of clustering than individual observations.
#' TBC!!
#' 
#' Parameters `b` are the mean function beta parameters, parameters `cov` are the
#' covariance function parameters in the same order as `$covariance$parameters`, and
#' parameters `d` are the estimated random effects.
#' @return TBC
#' @export
summary.mcml <- function(x,digits=2,...){
  pars <- print(x)
  ## summarise random effects
  dfre <- data.frame(Mean = round(apply(x$re.samps,2,mean),digits = digits), 
                     lower = round(apply(x$re.samps,2,function(i)quantile(i,0.025)),digits = digits),
                     upper = round(apply(x$re.samps,2,function(i)quantile(i,0.975)),digits = digits))
  colnames(dfre) <- c("Estimate","2.5% CI","97.5% CI")
  cat("Random effects estimates\n")
  print(dfre)
  ## add in model fit statistics
  return(invisible(list(coefficients = pars,re.terms = dfre)))
}

#' Prints a glmmr simuation output
#' 
#' Print method for class "`glmmr.sim`"
#' 
#' @param x an object of class "`mcml`" as a result of a call to MCML, see \link[glmmr]{Design}
#' @param digits Number of digits to print
#' @param ... Further arguments passed from other methods
#' @details 
#' `print.glmmr.sim` calculates multiple statistics summarising the design analysis.
#' 
#'  Simulation diagnostics. TBC
#' @return TBC
#' @export
print.glmmr.sim <- function(x, digits = 2,...){
  ## sim summary
  cat("glmmr simulation-based analysis\n",paste0(rep("-",31),collapse = ""),"\n")
  cat("Number of iterations: ",x$nsim,"\n")
  cat("Simulation method: ",x$sim_method,"\n")
  cat("For model with family",x$family[[1]],", link function",x$family[[2]],", ",
      x$n,"observations and \nMean function: ",as.character(x$mean_formula),"\nCovariance function: ",
      as.character(x$cov_formula),"\nTrue beta parameters: ",x$b_parameters,"\nCovariance parameters: ",unlist(x$cov_parameters),"\n")
  cat("\nSimulation diagnostics \n",paste0(rep("-",31),collapse = ""),
      "\nSimulation algorithm: ",x$mcml_method,"\n")
  conv <- mean(x$convergence)
  ## get coverage
  thresh <- qnorm(1-x$alpha/2)
  nbeta <- length(x$b_parameters)
  rows_to_include <- 1:nbeta
  cover <- Reduce(rbind,lapply(x$coefficients,function(i){
    (i$est[rows_to_include] - thresh*i$SE[rows_to_include]) <= x$b_parameters & 
      (i$est[rows_to_include] + thresh*i$SE[rows_to_include]) >= x$b_parameters
  }))
  cover <- colMeans(cover)
  cat("MCML algorithm convergence: ",round(conv*100,1),"%\nalpha: ",paste0(x$alpha*100,"%"),
      "\nCI coverage (beta):",paste0(round(cover*100,1),"%"))
  
  ## errors 
  cat("\n\n Errors\n",paste0(rep("-",31),collapse = ""),"\n")
  
  errdf <- sapply(1:nbeta,function(i)summarize.errors(x$coefficients,
                                  par = i,
                                  true = x$b_parameters[i],
                                  alpha = x$alpha))
  rownames(errdf) <- c("Type 2 (Power)","Type M (Exaggeration ratio)","Type S1 (Wrong sign)","Type S2 (Significant & wrong sign)")
  colnames(errdf) <- paste0("b",1:nbeta)
  print(apply(errdf,2,round,digits=digits))
  
  ## statistics
  cat("\n\n Distribution of statistics\n",paste0(rep("-",31),collapse = ""),"\np-values\n")
  statdf <- sapply(1:nbeta,function(i)summarize.stats(x$coefficients,
                                                      par = i,
                                                      alpha = x$alpha))
  pvdf <- apply(statdf,2,function(i)i$ptot)
  rownames(pvdf) <- c("0.00 - 0.01","0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.25", "0.25 - 0.50", "0.50 - 1.00")
  colnames(pvdf) <- paste0("b",1:nbeta)
  print(apply(pvdf,2,round,digits = digits))
  cat("\nConfidence interval half-width (+/-) quantiles\n")
  cidf <- apply(statdf,2,function(i)i$citot)
  colnames(cidf) <- paste0("b",1:nbeta)
  print(apply(cidf,2,round,digits = digits))
  
  ### dfbeta
  
  
  
  ##robustness
  cat("\n\n DFBETA for parameter: ",x$coefficients[[1]]$par[x$par],"\n",paste0(rep("-",41),collapse = ""),"\n")
  
  
   
  # dfb <- summarize.dfbeta(x$dfbeta,n=x$n)
  # cat("Mean minimum number of observations required to: \n\n")
  # dfbdf <- data.frame(x=c("Make estimate not significant","Change the sign of the estimate","Create wrong sign and significant estimate"),
  #                   Number = round(c(mean(
  # dfb[[1]]),mean(dfb[[3]]),mean(dfb[[5]])),digits = digits),
  #                   Proportion = round(c(mean(dfb[[2]]),mean(dfb[[4]]),mean(dfb[[6]])),digits = digits))
  # print(dfbdf)
}

#' Method to summarise errors 
#' 
#' Method to summarise errors from the simulation output of `$analysis()` from the Design class.
#' Not generally required by the user.
#' 
#' @param out list of mcml model outputs
#' @param par the index of the parameter to summarise
#' @param true true value of the parameter to summarise
#' @param alpha the type I error value
#' @return A vector with errors of type 2, M, S1, and S2
#' @importFrom stats qnorm setNames
summarize.errors <- function(out,
                              par,
                              true,
                              alpha){
  # generate t-stats and p-vals
  tstats <- lapply(out,function(i)i$est/i$SE)
  thresh <- abs(qnorm(alpha/2))
  pstats <- lapply(tstats,function(i)abs(i)>thresh)
  tstats <- Reduce(rbind,tstats)
  pstats <- Reduce(rbind,pstats)

  bests <- Reduce(rbind,lapply(out,function(i)i$est))

  # power
  pwr <- mean(pstats[,par])

  # type M
  ss.ests <- bests[,par][pstats[,par]]
  m.err <- NA
  if(length(ss.ests)>0){
    ss.ests.sgn <- ss.ests[sign(ss.ests)==sign(true)]
    if(length(ss.ests.sgn)>0){
      m.err <- mean(abs(ss.ests.sgn))/abs(true)
    }
  }

  # type S2
  s.err1 <- NA
  s.err1 <- mean(sign(bests[,par])!=sign(true))
  
  # type S2
  s.err2 <- NA
  if(length(ss.ests)>0){
    s.err2 <- mean(sign(ss.ests)!=sign(true))
  }

  setNames(c(pwr,m.err,s.err1,s.err2),c("p","m","s1","s2"))

}


#' Method to summarise statistics 
#' 
#' Method to summarise statistics from the simulation output of `$analysis()` from the Design class.
#' Not generally required by the user.
#' 
#' @param out list of mcml model outputs
#' @param par the index of the parameter to summarise
#' @param alpha the type I error value
#' @return A list containing two names data frames (`ptot` and `citot`) summarising 
#' the values of p-statistics from the model fits, and the quantiles of 1-alpha% confidence
#' interval half-widths, respectively.
#' @importFrom stats pnorm qnorm setNames
summarize.stats <- function(out,
                             par,
                            alpha=0.05){

  tstats <- lapply(out,function(i)i$est/i$SE)
  pstats <- lapply(tstats,function(i)(1-pnorm(abs(i)))*2)
  pstats <- Reduce(rbind,pstats)
  pstats <- pstats[,par]
  thresh <- qnorm(1-alpha/2)
  ses <- Reduce(rbind,lapply(out,function(i)i$SE))
  cis <- ses*thresh
  cis <- cis[,par]

  pgr <- c(0,0.01,0.05,0.1,0.25,0.5,1)
  ptot <- rep(NA,6)
  for(i in 1:6)ptot[i] <- mean((pstats >= pgr[i] & pstats < pgr[i+1]))

  citot <- quantile(cis,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))

  return(list(ptot = ptot,citot = citot))
}


#' Method to summarise DFBETA output
#' 
#' Method to summarise statistics from the simulation output of `$dfbeta()` from the Design class.
#' Not generally required by the user.
#' 
#' @param out list of mcml model outputs
#' @param n Total number of observations in the model
#' @return A list containing the number of observations and proportion of observations
#' required to change significance, sign, and significant sign from the models.
summarize.dfbeta <- function(out,
                             n){
  
  #significance change
  n.sig <- drop(Reduce(rbind,lapply(out,function(i)i[[1]])))
  p.sig <- n.sig/n
  
  # sign change
  n.sign <- drop(Reduce(rbind,lapply(out,function(i)i[[2]])))
  p.sign <- n.sign/n
  
  # sigsign change
  n.sigsign <- drop(Reduce(rbind,lapply(out,function(i)i[[3]])))
  p.sigsign <- n.sigsign/n
  
  return(list(n.sig, p.sig, n.sign, p.sign, n.sigsign, p.sigsign))
}

