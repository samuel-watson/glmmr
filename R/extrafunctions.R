

# function to identify group membership
match_rows <- function(x,target,by){
  if(ncol(target)==1){
    tstr <- target[,by]
    xstr <- x[,by]
  } else {
    xstr <- Reduce(paste0,as.data.frame(apply(x[,by],2,function(i)paste0(i,".0000."))))
    tstr <- Reduce(paste0,as.data.frame(apply(target[,by],2,function(i)paste0(i,".0000."))))
  }
  Z <- matrix(0,nrow=length(xstr),ncol=length(tstr))
  mat <- lapply(tstr,function(i)which(xstr==i))
  for(i in 1:length(mat))Z[mat[[i]],i] <- 1
  return(Z)
}

## model non-linear functons below

fexp <- function(x){
  if(length(x$pars)!=2)stop("two parameters required for fexp")
  x$pars[1]*exp(-x$pars[2]*x$data)
}

pexp <- function(x){
  x$pars[1]^x$data
}

gr <- function(x){
  I(x$data==0)*x$pars[1]^2
}

## create block matrix

blockmat <- function(...){
  matlist <- list(...)
  n <- length(matlist)
  N <- 0:(n-1)
  rlist <- list()
  for(i in 1:n){
    N <- (N+1)%%n
    N[N==0] <- n
    rlist[[i]] <- Reduce(cbind,matlist[N])
  }
  Reduce(rbind,rlist)
}

log_mvnd <- function(s,logdet,invD){
  Matrix::drop(-(length(s)/2)*log(2*pi)-0.5*logdet-0.5*Matrix::t(s)%*%invD%*%s)
}