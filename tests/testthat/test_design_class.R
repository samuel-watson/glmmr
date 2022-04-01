testthat::test_that("design class summary functions",{
  df <- nelder(~ ((int(2)*t(3)) > cl(3)) > ind(5))
  df$int <- df$int - 1
  mf1 <- MeanFunction$new(
    formula = ~ int + factor(t) - 1,
    data=df,
    parameters = rep(0,4),
    family = gaussian()
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)) + (1|gr(cl*t)),
    parameters = c(0.25,0.1)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  expect_equal(sum(des$fitted()),0)
  expect_equal(des$n(),90)
  ncl <- des$n_cluster()
  expect_equal(ncl[1,2],18)
  expect_s3_class(des$sim_data(),"numeric")
  expect_equal(length(des$sim_data()),90)
  
  des$subset_cols(1:3)
  expect_equal(ncol(des$mean_function$X),3)
  des$subset_rows(1:50)
  expect_equal(des$n(),50)
})

testthat::test_that("MCML tests",{
  df <- nelder(~(cl(10)*t(5)) > ind(10))
  df$int <- 0
  df[df$cl > 5, 'int'] <- 1
  mf1 <- MeanFunction$new(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family =gaussian()
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)),
    parameters = c(0.25)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  ysim <- des$sim_data()
  out <- des$MCML(y = ysim,
                  se.method = "none",
                  verbose=FALSE)
  expect_s3_class(out$coefficients$est,"numeric")
  expect_equal(length(out$coefficients$est),8)
  expect_true(!is.na(out$aic))
  expect_true(!is.na(out$Rsq$cond))
  expect_true(!is.na(out$Rsq$marg))
  
  out <- des$MCML(y = ysim,
                  se.method = "none",
                  method="mcnr",
                  verbose=FALSE)
  expect_s3_class(out$coefficients$est,"numeric")
  expect_equal(length(out$coefficients$est),8)
  
  out <- des$MCML(y = ysim,
                  se.method = "approx",
                  verbose=FALSE)
  expect_s3_class(out$coefficients$SE,"numeric")
  expect_true(out$coefficients$SE[1]>0 & out$coefficients$SE[1]<1)
  
  hess <- do.call(f_lik_hess,append(des$covariance$.__enclos_env__$private$D_data,
                                    list(as.matrix(des$covariance$Z),
                                         as.matrix(des$mean_function$X),
                                         ysim,
                                         out$re.samps,
                                         out$coefficients$est[8],
                                         family=des$mean_function$family[[1]],
                                         link=des$mean_function$family[[2]],
                                         start = out$coefficients$est[1:8],
                                         lower = c(rep(-Inf,6),rep(1e-6,2)),
                                         upper = c(rep(Inf,8)),
                                         tol=1e-4)))
  expect_s3_class(hess,"matrix")
  grad <- do.call(f_lik_grad,append(des$covariance$.__enclos_env__$private$D_data,
                                    list(as.matrix(des$covariance$Z),
                                         as.matrix(des$mean_function$X),
                                         ysim,
                                         out$re.samps,
                                         out$coefficients$est[8],
                                         family=des$mean_function$family[[1]],
                                         link=des$mean_function$family[[2]],
                                         start = out$coefficients$est[1:8],
                                         lower = c(rep(-Inf,6),rep(1e-6,2)),
                                         upper = c(rep(Inf,8)),
                                         tol=1e-4)))
  expect_s3_class(grad,"matrix")
  
  out <- des$MCML(y = des$sim_data(),
                  se.method = "perm",
                  permutation.par = 6,
                  verbose=FALSE,
                  options = list(
                    perm_type="unw",
                    perm_iter=10,
                    perm_parallel=FALSE,
                    perm_ci_steps=10
                  ))
  expect_true(out$coefficients$SE[6]>0 & out$coefficients$SE[6]<1)
  out <- des$MCML(y = des$sim_data(),
                  se.method = "none",
                  options = list(
                    sim_lik_step=TRUE
                  ))
  expect_s3_class(out$coefficients$est,"numeric")
  expect_equal(length(out$coefficients$est),8)
  
  des <- Design$new(
    covariance = cov1,
    mean.function = list(
      formula = ~ factor(t) + int - 1,
      data=df,
      parameters = c(rep(0,5),0.6),
      family =binomial()
    ),
    var_par = 1
  )
  ysim <- des$sim_data()
  out <- des$MCML(y = ysim,
                  se.method = "none",
                  verbose=FALSE)
  expect_s3_class(out$coefficients$est,"numeric")
  expect_equal(length(out$coefficients$est),8)
  
  des <- Design$new(
    covariance = cov1,
    mean.function = list(
      formula = ~ factor(t) + int - 1,
      data=df,
      parameters = c(rep(0,5),0.6),
      family =poisson()
    ),
    var_par = 1
  )
  ysim <- des$sim_data()
  out <- des$MCML(y = ysim,
                  se.method = "none",
                  verbose=FALSE)
  expect_s3_class(out$coefficients$est,"numeric")
  expect_equal(length(out$coefficients$est),8)
  
})

testthat::test_that("analysis tests",{
  df <- nelder(~(cl(10)*t(5)) > ind(10))
  df$int <- 0
  df[df$cl > 5, 'int'] <- 1
  mf1 <- MeanFunction$new(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family =gaussian()
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)),
    parameters = c(0.25)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  
  test <- des$analysis(type="sim",
                       iter=2,
                       par=6,
                       parallel = FALSE,
                       verbose = FALSE,
                       options = list(no_warnings=TRUE),
                       m = 10)
  
  expect_s3_class(test,"glmmr.sim")
  expect_s3_class(test$coefficients[[1]]$est,"numeric")
  expect_s3_class(test$coefficients[[1]]$SE,"numeric")
  expect_true(test$coefficients[[1]]$SE[1]>0 & test$coefficients[[1]]$SE[1]<1)
  expect_s3_class(test$dfbeta[[1]],"matrix")
  expect_true(all(!is.na(test$dfbeta[[1]])))
  
  test <- des$analysis(type="sim_approx",
                       iter=2,
                       par=6,
                       parallel = FALSE,
                       verbose = FALSE,
                       options = list(no_warnings=TRUE),
                       m = 10)
  
  expect_s3_class(test$coefficients[[1]]$est,"numeric")
  expect_s3_class(test$coefficients[[1]]$SE,"numeric")
  expect_true(test$coefficients[[1]]$SE[1]>0 & test$coefficients[[1]]$SE[1]<1)
  expect_s3_class(test$dfbeta[[1]],"matrix")
  expect_true(all(!is.na(test$dfbeta[[1]])))
  
  cov2 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl))+(1|gr(cl*t)),
    parameters = c(0.25,0.1)
  )
  des2 <- Design$new(
    covariance = cov2,
    mean.function = mf1,
    var_par = 1
  )
  
  test <- des$analysis(type="sim",
                       iter=2,
                       par=6,
                       sim_design = des2,
                       parallel = FALSE,
                       verbose = FALSE,
                       options = list(no_warnings=TRUE),
                       m = 10)
  
  expect_true(all(!is.na(test$sim_mean_formula)))
  expect_s3_class(test,"glmmr.sim")
  expect_s3_class(test$coefficients[[1]]$est,"numeric")
  expect_s3_class(test$coefficients[[1]]$SE,"numeric")
  expect_true(test$coefficients[[1]]$SE[1]>0 & test$coefficients[[1]]$SE[1]<1)
  expect_s3_class(test$dfbeta[[1]],"matrix")
  expect_true(all(!is.na(test$dfbeta[[1]])))
})

testthat::test_that("MCMC tests",{
  df <- nelder(~(cl(6)*t(5)) > ind(5))
  df$int <- 0
  df[df$cl > 3, 'int'] <- 1
  mf1 <- MeanFunction$new(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family =gaussian()
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)),
    parameters = c(0.25)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  
  ysim <- des$sim_data()
  capture.output(suppressWarnings(nfit1 <- des$MCMC(y=ysim,
                                                    warmup_iter = 10,
                                                    sampling_iter = 10,
                                                    parallel = FALSE,
                                                    chains=1)),file=tempfile())
  expect_s3_class(nfit1,"stanfit")
  
})

testthat::test_that("bayesian analysis tests",{
  df <- nelder(~(cl(6)*t(5)) > ind(5))
  df$int <- 0
  df[df$cl > 3, 'int'] <- 1
  mf1 <- MeanFunction$new(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family =gaussian()
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)*ar1(t)),
    parameters = c(0.25,0.8)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  
  cov2 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl))+(1|gr(cl*t)),
    parameters = c(0.25,0.1)
  )
  des2 <- Design$new(
    covariance = cov2,
    mean.function = mf1,
    var_par = 1
  )
  
  test <- des$analysis_bayesian(10,6,warmup_iter = 20,sampling_iter = 20,
                                 parallel = FALSE,verbose=FALSE)
  
  expect_s3_class(test,"glmmr.sim")
  expect_s3_class(test$posterior_var,"numeric")
  expect_equal(length(test$posterior_var),10)
  expect_s3_class(test$sbc_ranks,"numeric")
  expect_equal(length(test$sbc_ranks),10)
  expect_true(all(test$sbc_ranks >=0 & test$sbc_ranks <= 20))
  expect_s3_class(test$posterior_threshold,"numeric")
  expect_equal(length(test$posterior_threshold),10)
  expect_true(all(test$posterior_threshold>=0 & test$posterior_threshold <=1))
  
  test <- des$analysis_bayesian(10,6,warmup_iter = 20,sampling_iter = 20,
                                 sim_design = des2,
                                 parallel = FALSE,verbose=FALSE)
  
  expect_s3_class(test,"glmmr.sim")
  expect_s3_class(test$posterior_var,"numeric")
  expect_equal(length(test$posterior_var),10)
  expect_s3_class(test$sbc_ranks,"numeric")
  expect_equal(length(test$sbc_ranks),10)
  expect_true(all(test$sbc_ranks >=0 & test$sbc_ranks <= 20))
  expect_s3_class(test$posterior_threshold,"numeric")
  expect_equal(length(test$posterior_threshold),10)
  expect_true(all(test$posterior_threshold>=0 & test$posterior_threshold <=1))
  
})

testthat::test_that("permutation tests",{
  
  df <- nelder(~(cl(6)*t(5)) > ind(5))
  df$int <- 0
  df[df$cl > 3, 'int'] <- 1
  
  treatf <- function(){
    tr <- sample(rep(c(0,1),each=3),6,replace = FALSE)
    rep(tr,each=25)
  }
  
  mf1 <- MeanFunction$new(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family =gaussian(),
    treat_var = "int",
    random_function = treatf
  )
  cov1 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)),
    parameters = c(0.25)
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  
  perm1 <- des$permutation_test(ysim,6,0.61,type="unw",
                                iter = 10,nsteps = 10,verbose = FALSE,
                                start=0.8)
  expect_true(perm1$p>=0,perm1$p<=1)
  expect_s3_class(perm1$lower,"numeric")
  
  perm1 <- des$permutation_test(ysim,6,0.61,type="cov",
                                iter = 10,nsteps = 10,verbose = FALSE,
                                start=0.8)
  expect_true(perm1$p>=0,perm1$p<=1)
  expect_s3_class(perm1$lower,"numeric")
  
})