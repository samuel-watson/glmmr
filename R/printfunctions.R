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

prnt.errors <- function(errors,digits=2){
  cat("Errors\n----------------------------------\n")

  errdf <- data.frame(Value = round(errors,digits))
  rownames(errdf) <- c("Type 2 (Power)","Type M (Exaggeration ratio)","Type S1 (Wrong sign)","Type S2 (Significant & wrong sign)")
  print(errdf)
}

summarize.stats <- function(out,
                             par,
                            alpha=0.05){

  tstats <- lapply(out,function(i)i$est/i$SE)
  pstats <- lapply(tstats,function(i)(1-pnorm(abs(i)))*2)
  pstats <- Reduce(rbind,pstats)
  pstats <- pstats[,par]
  thresh <- qnorm(1-alpha/2)
  ses <- Reduce(rbind,lapply(out,function(i)i$SE))
  cis <- ses*thresh*2
  cis <- cis[,par]

  pgr <- c(0,0.01,0.05,0.1,0.25,0.5,1)
  ptot <- rep(NA,6)
  for(i in 1:6)ptot[i] <- mean((pstats >= pgr[i] & pstats < pgr[i+1]))

  citot <- quantile(cis,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))

  return(list(ptot = ptot,citot = citot))
}

prnt.stats <- function(x,digits=2){
  cat("P-value distribution \n-----------------------\n")
  pvdf <- data.frame(Proportion = round(x$ptot,digits))
  rownames(pvdf) <- c("0.00 - 0.01","0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.25", "0.25 - 0.50", "0.50 - 1.00")
  print(pvdf)

  cat("\n\nCI width quantiles \n---------------------\n")
  cidf <- data.frame(Quantile = round(unname(x$citot),digits))
  rownames(cidf ) <- names(x$citot)
  print(cidf)
}

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

prnt.dfbeta <- function(x, digits = 2){
  cat("Robustness \n-----------------------\n")
  cat("Mean minimum number of observations required to: \n\n")
  dfb <- data.frame(x=c("Make estimate not significant","Change the sign of the estimate","Create wrong sign and significant estimate"),
                    Number = round(c(mean(x[[1]]),mean(x[[3]]),mean(x[[5]])),digits = digits),
                    Proportion = round(c(mean(x[[2]]),mean(x[[4]]),mean(x[[6]])),digits = digits))
  print(dfb)
  
}

