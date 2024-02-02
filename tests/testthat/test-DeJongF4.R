library(testthat)
library(xegaSelectGene)

test_that("Problem environment DeJongF4 produces the correct results",
  {expect_identical(DeJongF4Factory()$name(), "DeJongF4")
  expect_setequal(DeJongF4Factory()$bitlength(), rep(64, 30))
  expect_setequal(DeJongF4Factory()$genelength(), 1920)
  expect_setequal(DeJongF4Factory()$lb(), rep(-1.28, 30))
  expect_setequal(DeJongF4Factory()$ub(), rep(1.28, 30))
  expect_output(DeJongF4Factory()$describe(), regexp="De Jong")
  expect_setequal(DeJongF4Factory()$maxp(), c(1.28, -1.28))
  expect_equal(DeJongF4Factory()$solution()$minimum, 0.0)
  expect_equal(DeJongF4Factory()$solution()$maximum, 1248.225)
  expect_setequal(DeJongF4Factory()$solution()$minpoints[[1]], rep(0, 30))
  expect_setequal(DeJongF4Factory()$solution()$maxpoints[[1]], rep(1.28, 30))
  }
) 

