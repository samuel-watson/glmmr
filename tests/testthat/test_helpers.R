testthat::test_that("stepped wedge helper function",{
  des <- stepped_wedge(6,10,icc=0.05)
  expect_s3_class(des,"Design")
  des <- stepped_wedge(6,10,icc=0.05,cac=0.8)
  expect_s3_class(des,"Design")
  des <- stepped_wedge(6,10,icc=0.05,cac=0.8,iac=0.1)
  expect_s3_class(des,"Design")
})

testthat::test_that("parallel helper function",{
  des <- parallel_crt(6,10,icc=0.05)
  expect_s3_class(des,"Design")
  des <- parallel_crt(6,10,icc=0.05,cac=0.8)
  expect_s3_class(des,"Design")
  des <- parallel_crt(6,10,icc=0.05,cac=0.8,iac=0.1)
  expect_s3_class(des,"Design")
})

testthat::test_that("staircase helper function",{
  des <- staircase_crt(6,10,icc=0.05)
  expect_s3_class(des,"Design")
  des <- staircase_crt(6,10,icc=0.05,cac=0.8)
  expect_s3_class(des,"Design")
  des <- staircase_crt(6,10,icc=0.05,cac=0.8,iac=0.1)
  expect_s3_class(des,"Design")
})