.summarize.errors <- function(out,
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
      m.err <- mean(ss.ests.sgn)/true
    }
  }
  
  # type S
  s.err <- NA
  if(length(ss.ests)>0){
    s.err <- mean(sign(ss.ests)!=sign(true))
  }
  
  setNames(c(pwr,m.err,s.err),c("p","m","s"))
  
}

.prnt.errors <- function(errors,digits=2){
  cat("Errors\n----------------------------------\n")
  
  errdf <- data.frame(Value = round(errors,digits))
  rownames(errdf) <- c("Type 2 (Power)","Type M (Exaggeration ratio)","Type S (Wrong sign)")
  print(errdf)
}

.summarize.stats <- function(out,
                             par){
  
  tstats <- lapply(out,function(i)i$est/i$SE)
  pstats <- lapply(tstats,function(i)(1-pnorm(abs(i)))*2)
  pstats <- Reduce(rbind,pstats)
  pstats <- pstats[,par]
  
  ses <- Reduce(rbind,lapply(out,function(i)i$SE))
  cis <- ses*thresh*2
  cis <- cis[,par]
  
  pgr <- c(0,0.01,0.05,0.1,0.25,0.5,1)
  ptot <- rep(NA,6)
  for(i in 1:6)ptot[i] <- mean((pstats >= pgr[i] & pstats < pgr[i+1]))
  
  citot <- quantile(cis,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
  
  return(list(ptot = ptot,citot = citot))
}

.prnt.stats <- function(x,digits=2){
  cat("P-value distribution \n-----------------------\n")
  pvdf <- data.frame(Proportion = round(x$ptot,digits))
  rownames(pvdf) <- c("0.00 - 0.01","0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.25", "0.25 - 0.50", "0.50 - 1.00")
  print(pvdf)
  
  cat("\n\nCI width quantiles \n---------------------\n")
  cidf <- data.frame(Quantile = round(unname(x$citot),digits))
  rownames(cidf ) <- names(x$citot)
  print(cidf)
}
