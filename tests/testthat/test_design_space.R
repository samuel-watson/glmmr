testthat::test_that("set up design space options",{
  des <- stepped_wedge(5,5,icc=0.05)
  ds <- DesignSpace$new(des)
  expect_equal(sum(  ds$experimental_condition-seq_len(150)),0)
  expect_equal(ds$weights,1)
  des2 <- stepped_wedge(5,5,icc=0.05,cac=0.8)
  ds <- DesignSpace$new(des,des2)
  expect_equal(sum(  ds$experimental_condition-seq_len(150)),0)
  expect_true(all(ds$weights==0.5))
})

testthat::test_that("test optimal design",{
  df <- nelder(~(cl(6)*t(5)) > ind(5))
  df$int <- 0
  df[df$t >= df$cl, 'int'] <- 1
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
  ds <- DesignSpace$new(des)
  
  opt <- ds$optimal(30,C=c(rep(0,5),1),verbose=FALSE)
  
  expect_equal(length(opt),30)
  expect_true(all(opt[1:2]<=20))
  
  ds <- DesignSpace$new(des,experimental_condition = df$cl)
  opt <- ds$optimal(30,C=c(rep(0,5),1),verbose=FALSE)
  
  expect_equal(round(opt[1],4),0.2561)
  
  opt <- ds$optimal(4,C=c(rep(0,5),1),verbose=FALSE,force_hill = TRUE)
  
  expect_true(all(opt == c(1,3,5,6)))
  
  cov2 <- Covariance$new(
    data = df,
    formula = ~ (1|gr(cl)*ar1(t)),
    parameters = c(0.25,0.8)
  )
  des2 <- Design$new(
    covariance = cov1,
    mean.function = mf1,
    var_par = 1
  )
  
  ds <- DesignSpace$new(des,des2)
  opt <- ds$optimal(30,C=list(c(rep(0,5),1),c(rep(0,5),1)),
                    verbose=FALSE,robust_function = "weighted")
  expect_equal(length(opt),30)
  expect_true(all(opt[1:2]<=20))
  
  opt <- ds$optimal(30,C=list(c(rep(0,5),1),c(rep(0,5),1)),
                    verbose=FALSE,robust_function = "minimax")
  expect_equal(length(opt),30)
  expect_true(all(opt[1:2]<=20))
  
})

