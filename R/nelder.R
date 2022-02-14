
#GENERATE CRISS CROSS
cross_df <- function(df1,df2){
  cnames <- c(colnames(df1),colnames(df2))
  df1 <- as.data.frame(df1[rep(1:nrow(df1),each=nrow(df2)),])
  df3 <- cbind(df1,df2[1:nrow(df2),])
  colnames(df3) <- cnames
  return(df3)
}

#GENERATE NESTING
nest_df <- function(df1,df2){
  df3 <- cbind(df1[rep(1:nrow(df1),each=nrow(df2)),],df2)
  colnames(df3)[1:ncol(df1)] <- colnames(df1)
  if(ncol(df1)>1)ids <- Reduce(paste0,df3[,1:ncol(df1)]) else ids <- df3[,1]
  df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))] <- apply(as.data.frame(df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))]),2,
                                                     function(i)as.numeric(as.factor(paste0(ids,i))))
  colnames(df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))]) <- colnames(df2)
  return(df3)
}

nelder <- function(formula){
  if(formula[[1]]=="~")formula <- formula[[2]]
  f1l <- formula[[2]]
  f1r <- formula[[3]]
  
  if(as.character(f1l[[1]])%in%c("*",">")){
    df1 <- Recall(f1l)
  } else if(as.character(f1l[[1]])%in%c("(")){
    df1 <- Recall(f1l[[2]])
  } else {
    df1 <- data.frame(a = seq(1,f1l[[2]]))
    colnames(df1) <- as.character(f1l[[1]])
  }
  
  if(as.character(f1r[[1]])%in%c("(")){
    df2 <- Recall(f1r[[2]])
  } else {
    df2 <- data.frame(a = seq(1,f1r[[2]]))
    colnames(df2) <- as.character(f1r[[1]])
  }
  
  if(formula[[1]] == "*"){
    df <- cross_df(df1,df2)
  } else if(formula[[1]] == ">"){
    df <- nest_df(df1,df2)
  }
  rownames(df) <- NULL
  return(df)
}

cycles <- function(a){
  fa <- a
  for(i in 1:(length(a)-1)){
    a <- c(a[2:length(a)],a[1])
    fa <- c(fa,a)
  }
  fa
}