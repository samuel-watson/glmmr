
staircase_crt <-  function(J,
                           M,
                           beta=c(rep(0,J+1),0),
                           icc,
                           cac = NULL,
                           iac = NULL,
                           var = 1,
                           family = gaussian()){
  if(missing(icc))stop("icc must be set as a minimum")
  
  ndesigns <- length(icc) * ifelse(!is.null(cac[1]),length(cac),1) *
    ifelse(!is.null(iac[1]),length(iac),1)
  
  if(!is.null(cac[1]) && !is.na(cac[1])){
    wp_var <- icc[1]*var*(1-cac[1])
    bp_var <- icc[1]*var[1]*cac[1]
  } else {
    bp_var <- icc[1]*var
  }
  if(!is.null(iac[1]) && !is.na(iac[1])){
    ind_var <- var*(1-icc[1])*iac[1]
    sigma <- var*(1-icc[1])*(1-iac[1])
  } else {
    sigma <- var*(1-icc[1])
  }
  t <- J+1
  if(!is.null(iac[1]) && !is.na(iac[1])){
    df <- nelder(formula(paste0("~ (J(",J,") > ind(",M,")) * t(",t,")")))
  } else {
    df <- nelder(formula(paste0("~ (J(",J,") * t(",t,")) > ind(",M,")")))
  }
  
  df <- df[order(df$J,df$t),]
  int <- c()
  for(i in 1:J){
    int <- c(int,rep(c(rep(0,t-(i)),rep(1,i)),1))
  }
  df$int <- rep(int,each=M)
  
  # remove the non-staircase observations
  df$J <- J+1 - df$J
  for(j in 1:J){
    df <- rbind(df[df$J == j&df$t %in%c(j,j+1),],
                df[df$J != j,])
  }
  df <- df[df$t %in% c(2:J),]
  df$t <- df$t - 1
  df <- df[order(df$J,df$t),]
  
  if(is.null(cac[1]) || is.na(cac[1])){
    if(is.null(iac[1]) || is.na(iac[1])){
      f1 <- "~(1|gr(J))"
      pars <- list(list(sqrt(bp_var)))
    } else {
      f1 <- "~(1|gr(J)) + (1|gr(ind))"
      pars <- list(list(sqrt(bp_var)),list(sqrt(ind_var)))
    }
  } else {
    if(is.null(iac[1]) || is.na(iac[1])){
      f1 <- "~ (1|gr(J)) + (1|gr(J*t))"
      pars <- list(list(sqrt(bp_var)),list(sqrt(wp_var)))
    } else {
      f1 <- "~ (1|gr(J)) + (1|gr(J*t)) + (1|gr(ind))"
      pars <- list(list(sqrt(bp_var)),list(sqrt(wp_var)),list(sqrt(ind_var)))
    }
  }
  
  d1 <- Design$new(
    covariance = list(
      data=df,
      formula = f1,
      parameters = pars
    ),
    mean.function = list(
      formula = "~ factor(t) + int - 1",
      data = df,
      family = family,
      parameters = beta
      
    ),
    var_par = sigma
  )
  
  if(ndesigns>1){
    ds1 <- DesignSpace$new(d1)
    if(is.null(cac))cac <- NA
    if(is.null(iac))iac <- NA
    dsvalues <- expand.grid(icc=icc,cac=cac,iac=iac)
    
    for(i in 1:(ndesigns-1)){

      if(!is.null(dsvalues$cac[i+1]) && !is.na(dsvalues$cac[i+1])){
        wp_var <- dsvalues$icc[i+1]*var*(1-dsvalues$cac[i+1])
        bp_var <- dsvalues$icc[i+1]*var*dsvalues$cac[i+1]
      } else {
        bp_var <- dsvalues$icc[i+1]*var
      }
      if(!is.null(dsvalues$iac[i+1]) && !is.na(dsvalues$iac[i+1])){
        ind_var <- var*(1-dsvalues$icc[i+1])*dsvalues$iac[i+1]
        sigma <- var*(1-dsvalues$icc[i+1])*(1-dsvalues$iac[i+1])
      } else {
        sigma <- var*(1-dsvalues$icc[i+1])
      }
      
      if(is.null(dsvalues$cac[i+1]) || is.na(dsvalues$cac[i+1])){
        if(is.null(dsvalues$iac[i+1]) || is.na(dsvalues$iac[i+1])){
          f1 <- "~(1|gr(J))"
          pars <- list(list(bp_var))
        } else {
          f1 <- "~(1|gr(J)) + (1|gr(ind))"
          pars <- list(list(bp_var),list(ind_var))
        }
      } else {
        if(is.null(dsvalues$iac[i+1]) || is.na(dsvalues$iac[i+1])){
          f1 <- "~ (1|gr(J)) + (1|gr(J*t))"
          pars <- list(list(bp_var),list(wp_var))
        } else {
          f1 <- "~ (1|gr(J)) + (1|gr(J*t)) + (1|gr(ind))"
          pars <- list(list(bp_var),list(wp_var),list(ind_var))
        }
      }
      
      
      ds1$add(
        Design$new(
          covariance = Covariance$new(
            data=df,
            formula = f1,
            parameters = pars
          ),
          mean.function = d1$mean_function$clone(),
          var_par = 1
        )
      )
      
      # $parameters <- pars
      # ds1$.__enclos_env__$private$designs[[i+1]]$covariance$formula <- f1
      
    }
    return(ds1)
  } else {
    return(d1)
  }
}
