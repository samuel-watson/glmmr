testthat::test_that("nelder function produces data frames",{
  df <- nelder(~cl(3) > t(4))
  expect_s3_class(df,"data.frame")
  df2 <- nelder(~(cl(4) * t(3)) > ind(5))
  expect_s3_class(df2,"data.frame")
})

testthat::test_that("covariance class creation works",{
  cov1 <- Covariance$new(
    data = nelder(~cl(3) > t(4)),
    formula = ~ (1|gr(cl)),
    parameters = list(list(0.05))
  )
  expect_s4_class(cov1$D,"dgCMatrix")
  expect_s4_class(cov1$Z,"dgCMatrix")
  expect_warning(Covariance$new(
    data = nelder(~(cl(4) * t(3)) > ind(5)),
    formula = ~ (1|gr(cl)*pexp(t)) + (1|gr(ind)),
    parameters = list(list(0.05,0.8),list(0.01))
  ))
  cov2 <- Covariance$new(
    data = nelder(~(cl(4) > ind(5)) * t(3)),
    formula = ~ (1|gr(cl)*pexp(t)) + (1|gr(ind)),
    parameters = list(list(0.05,0.8),list(0.01))
  )
  expect_s4_class(cov2$D,"dgCMatrix")
  expect_s4_class(cov2$Z,"dgCMatrix")
})

testthat::test_that("mean function class creation works",{
  df <- nelder(~ ((int(2)*t(3)) > cl(3)) > ind(5))
  df$int <- df$int - 1
  mf1 <- MeanFunction$new(
    formula = ~ int + factor(t) - 1,
    data=df,
    parameters = rep(0,4),
    family = gaussian()
  )
  expect_s4_class(mf1$X,"dgCMatrix")
  expect_equal(mf1$n(),90)
})

testthat::test_that("design class creation works",{
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
    formula = ~ (1|gr(cl)) + (1|gr(cl)*gr(t)),
    parameters = list(list(0.05),list(1,0.01))
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  expect_s4_class(des$Sigma,"dgCMatrix")
  expect_equal(des$n(),90)
})

testthat::test_that("design class functions",{
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
    formula = ~ (1|gr(cl)) + (1|gr(cl)*gr(t)),
    parameters = list(list(0.05),list(1,0.01))
  )
  des <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  expect_equal(round(des$power(par=1,value = 1,alpha = 0.05),3),0.986)
  
})