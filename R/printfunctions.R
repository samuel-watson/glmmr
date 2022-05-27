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
  
  cat("\nFixed effects formula :",x$mean_form)
  cat("\nCovariance function formula: ",x$cov_form)
  cat("\nFamily: ",x$family,", Link function:",x$link,"\n")
  cat("\nNumber of Monte Carlo simulations per iteration: ",x$m," with tolerance ",x$tol,"\n")
  semethod <- ifelse(x$permutation,"permutation test",
                     ifelse(x$robust,"robust",ifelse(x$hessian,"hessian","approx")))
  cat("P-value and confidence interval method: ",semethod,"\n\n")
  pars <- x$coefficients[!grepl("d",x$coefficients$par),c('est','SE','lower','upper')]
  z <- pars$est/pars$SE
  pars <- cbind(pars[,1:2],z=z,p=2*(1-pnorm(abs(z))),pars[,3:4])
  colnames(pars) <- c("Estimate","Std. Err.","z value","p value","2.5% CI","97.5% CI")
  rnames <- x$coefficients$par[!grepl("d",x$coefficients$par)]
  if(any(duplicated(rnames))){
    did <- unique(rnames[duplicated(rnames)])
    for(i in unique(did)){
      rnames[rnames==i] <- paste0(rnames[rnames==i],".",1:length(rnames[rnames==i]))
    }
  }
  rownames(pars) <- rnames
  pars <- apply(pars,2,round,digits = digits)
  print(pars)
  
  cat("\ncAIC: ",round(x$aic,digits))
  cat("\nApproximate R-squared: Conditional: ",round(x$Rsq[1],digits)," Marginal: ",round(x$Rsq[2],digits))
  
  #messages
  if(x$permutation)message("Permutation test used for one parameter, other SEs are not reported. SEs and Z values
are approximate based on the p-value, and assume normality.")
  #if(!x$hessian&!x$permutation)warning("Hessian was not positive definite, standard errors are approximate")
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
  if(all(is.na(x$sim_mean_formula))){
    cat("\nSimulation and analysis model:\nFamily",x$family[[1]],", link function",x$family[[2]],", ",
        x$n,"observations and \nMean function: ",as.character(x$mean_formula),"\nCovariance function: ",
        as.character(x$cov_formula))
  } else {
    cat("\nSimulation model:\nFamily",x$sim_family[[1]],", link function",x$sim_family[[2]],", ",
        x$n,"observations and \nMean function: ",as.character(x$sim_mean_formula),"\nCovariance function: ",
        as.character(x$sim_cov_formula))
    cat("\n\nAnalysis model:\nFamily",x$family[[1]],", link function",x$family[[2]],", ",
        x$n,"observations and \nMean function: ",as.character(x$mean_formula),"\nCovariance function: ",
        as.character(x$cov_formula))
  }
  
  if(x$sim_method!="bayesian"){
    cat("\n\nTrue beta parameters: ",x$b_parameters,"\nTrue covariance parameters: ",
        unlist(x$cov_parameters),"\n")
  } else {
    if(any(!is.na(x$priors_sim))){
      sim_priors <- x$priors_sim
      cat("\n\nSimulation model priors:\nBeta: ",paste0("~N([",paste0(sim_priors$prior_b_mean,collapse=","),
                                                                         "],[",paste0(sim_priors$prior_b_sd,collapse=","),"])^2)"))
      cat("\nCovariance parameters: ",paste0("~N_+([",paste0(rep(0,length(sim_priors$prior_g_sd)),collapse=","),
                                                   "],[",paste0(sim_priors$prior_g_sd,collapse=","),"]^2)"))
      if(x$family[[1]]=="gaussian"){
        cat("\nScale parameter: ",paste0("~N_+(0,",sim_priors$sigma_sd,"^2)"))
      }
    }
    sim_priors <- x$priors
    cat("\n\nAnalysis model priors:\nBeta: ",paste0("~N([",paste0(sim_priors$prior_b_mean,collapse=","),
                                             "],[",paste0(sim_priors$prior_b_sd,collapse=","),"])^2)"))
    cat("\nCovariance parameters: ",paste0("~N_+([",paste0(rep(0,length(sim_priors$prior_g_sd)),collapse=","),
                                                 "],[",paste0(sim_priors$prior_g_sd,collapse=","),"]^2)"))
    if(x$family[[1]]=="gaussian"){
      cat("\nScale parameter: ",paste0("~N_+(0,",sim_priors$sigma_sd,"^2)"))
    }
  }
  
  if(x$sim_method!="bayesian"){
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
    if(x$sim_method=="full.sim"){
      rsq <- Reduce(rbind,x$Rsq)
      raic <- paste0(round(range(unlist(x$aic)),digits = digits),collpase=",")
      rmsq <- paste0(round(range(rsq[,1]),digits = digits),collpase=",")
      rcsq <- paste0(round(range(rsq[,2]),digits = digits),collpase=",")
    } else {
      rsq <- NA
      raic <- NA
      rmsq <- NA
      rcsq <- NA
    }
    
    cat("\nRange of cAIC: ",raic,
        "\nRange of: conditional R-squared",rmsq,
        " marginal R-squared: ",rcsq)
    
    ## errors 
    cat("\n\n Errors\n",paste0(rep("-",31),collapse = ""),"\n")
    
    errdf <- sapply(1:nbeta,function(i)summarize.errors(x$coefficients,
                                                        par = i,
                                                        true = x$b_parameters[i],
                                                        alpha = x$alpha))
    rownames(errdf) <- c("Type 2 (Power)","Type M (Exaggeration ratio)",
                         "Type S1 (Wrong sign)","Type S2 (Significant & wrong sign)","Bias")
    colnames(errdf) <- x$coefficients[[1]]$par[1:nbeta]#paste0("b",1:nbeta)
    print(apply(errdf,2,round,digits=digits))
    
    ## statistics
    cat("\n\n Distribution of statistics\n",paste0(rep("-",31),collapse = ""),"\np-values\n")
    statdf <- sapply(1:nbeta,function(i)summarize.stats(x$coefficients,
                                                        par = i,
                                                        alpha = x$alpha))
    pvdf <- apply(statdf,2,function(i)i$ptot)
    rownames(pvdf) <- c("0.00 - 0.01","0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.25", "0.25 - 0.50", "0.50 - 1.00")
    colnames(pvdf) <- x$coefficients[[1]]$par[1:nbeta]
    print(apply(pvdf,2,round,digits = digits))
    cat("\nConfidence interval half-width (+/-) quantiles\n")
    cidf <- apply(statdf,2,function(i)i$citot)
    colnames(cidf) <- x$coefficients[[1]]$par[1:nbeta]
    print(apply(cidf,2,round,digits = digits))
    
    ### dfbeta
    
    
    
    ##robustness
    cat("\n\n Deletion diagnostics for parameter: ",x$coefficients[[1]]$par[x$par],"\n",paste0(rep("-",41),collapse = ""),"\n")
    
    dfb <- summarize.dfbeta(x$dfbeta)
    cat("Mean maximum DFBETA: ",signif(dfb$maxb,digits = digits),
        "\nRange of mean DFBETA for each observation: ",signif(dfb$dfbrange,digits=digits),
        "\nObservations with largest DFBETA: ",dfb$maxobs)
    # cat("Mean minimum number of observations required to: \n\n")
    # dfbdf <- data.frame(x=c("Make estimate not significant","Change the sign of the estimate","Create wrong sign and significant estimate"),
    #                   Number = round(c(mean(
    # dfb[[1]]),mean(dfb[[3]]),mean(dfb[[5]])),digits = digits),
    #                   Proportion = round(c(mean(dfb[[2]]),mean(dfb[[4]]),mean(dfb[[6]])),digits = digits))
    # print(dfbdf)
  } else {
    cat("\n\nSimulation diagnostics \n",paste0(rep("-",31),collapse = ""))
    cat("\nUse plot() to view the simulation based calibration plot\n")
    cat("\nQuantiles of the posterior variance:\n")
    print(round(quantile(x$posterior_var,seq(0,1,by=0.1)),digits = digits))
    cat("\nProbability of the probability that the parameter will be greater than ",x$threshold,":\n")
    print(round(quantile(x$posterior_threshold,seq(0,1,by=0.1)),digits = digits))
  }
}

#' Plotting method for glmmr.sim
#' 
#' Plots a glmmr.sim object
#' 
#' For a Bayesian simulation analysis, this will plot the simulation based calibration ranks and the 
#' distribution of the posterior variance. Otherwise, it will plot the distribution of confidence interval 
#' half widths.
#' @param x A `glmmr.sim` object
#' @param par Integer indicating the index of the parameter to plot if not a Bayesian analysis
#' @param alpha Numeric indictaing the type I error rate for non-Bayesian analysis
#' @return A `ggplot2` plot
#' @examples 
#' ...
#' @export
plot.glmmr.sim <- function(x,
                          par,
                          alpha=0.05){
  if(x$type==1){
    out <- summarize.stats(x$coefficients,par,alpha)
    dfp <- data.frame(stat = rep(c("p-values","CI half-width"),each=x$iter),
                      value = c(out$pstats,out$cis))
    ggplot2::ggplot(aes(x=value),data=dfp)+
      ggplot2::facet_wrap(~stat,scales = "free")+
      ggplot2::geom_histogram()+
      ggplot2::labs(x="Value",y="Count")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid=ggplot2::element_blank())
  } else {
    dfp <- data.frame(stat = rep(c("Posterior variance","Rank"),each=x$iter),
                      value = c(x$posterior_var,x$sbc_ranks))
    ggplot2::ggplot(aes(x=value),data=dfp)+
      ggplot2::facet_wrap(~stat,scales = "free")+
      ggplot2::geom_histogram()+
      ggplot2::labs(x="Value",y="Count")+
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid=ggplot2::element_blank())
  }
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
  
  #bias 
  bias <- mean(bests[,par]) - true

  setNames(c(pwr,m.err,s.err1,s.err2,bias),c("p","m","s1","s2","bias"))

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

  return(list(ptot = ptot,citot = citot,pstats=pstats,cis=cis))
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
summarize.dfbeta <- function(out){
  
  # max values
  maxb <- unlist(lapply(out,function(i)max(abs(i))))
  
  #individual obs - do any have consistent leverage?
  dfball <- colMeans(Reduce(rbind,out))
  dfbrange <- range(dfball)
  maxobs <- order(abs(dfball))[1:4]
  return(list(maxb=mean(maxb),dfbrange=dfbrange,maxobs=maxobs))
  
  # #significance change
  # n.sig <- drop(Reduce(rbind,lapply(out,function(i)i[[1]])))
  # p.sig <- n.sig/n
  # 
  # # sign change
  # n.sign <- drop(Reduce(rbind,lapply(out,function(i)i[[2]])))
  # p.sign <- n.sign/n
  # 
  # # sigsign change
  # n.sigsign <- drop(Reduce(rbind,lapply(out,function(i)i[[3]])))
  # p.sigsign <- n.sigsign/n
  # 
  # return(list(n.sig, p.sig, n.sign, p.sign, n.sigsign, p.sigsign))
}

