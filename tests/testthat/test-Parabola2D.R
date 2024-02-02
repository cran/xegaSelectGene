library(testthat)
library(xegaSelectGene)

test_that("Problem environment DelayedP produces the correct results",
  {expect_identical(DelayedPFactory()$name(), "DelayedP")
  expect_setequal(DelayedPFactory()$bitlength(), rep(20, 2))
  expect_setequal(DelayedPFactory()$genelength(), 40)
  expect_setequal(DelayedPFactory()$lb(), rep(-4.5, 2))
  expect_setequal(DelayedPFactory()$ub(), rep(4.5, 2))
  expect_equal(DelayedPFactory()$f(c(1,1)), 2)
  expect_output(DelayedPFactory()$describe(), regexp="parabola")
  expect_equal(DelayedPFactory()$solution()$minimum, 0.0)
  expect_equal(DelayedPFactory()$solution()$maximum, 40.5)
  expect_setequal(DelayedPFactory()$solution()$minpoints[[1]], rep(0, 2))
  expect_setequal(DelayedPFactory()$solution()$maxpoints[[1]], rep(4.5, 2))
  }
) 

test_that("Problem environment Parabola2D produces the correct results",
  {expect_identical(Parabola2DFactory()$name(), "Parabola2D")
  expect_setequal(Parabola2DFactory()$bitlength(), rep(20, 2))
  expect_setequal(Parabola2DFactory()$genelength(), 40)
  expect_setequal(Parabola2DFactory()$lb(), rep(-4.5, 2))
  expect_setequal(Parabola2DFactory()$ub(), rep(4.5, 2))
  expect_equal(Parabola2DFactory()$f(c(1,1)), 2)
  expect_output(Parabola2DFactory()$describe(), regexp="unimodal")
  expect_equal(Parabola2DFactory()$solution()$minimum, 0.0)
  expect_equal(Parabola2DFactory()$solution()$maximum, 40.5)
  expect_setequal(Parabola2DFactory()$solution()$minpoints[[1]], rep(0, 2))
  expect_setequal(Parabola2DFactory()$solution()$maxpoints[[1]], rep(4.5, 2))
  }
) 

test_that("Problem environment Parabola2DEarly produces the correct results",
  {expect_identical(Parabola2DEarlyFactory()$name(), "Parabola2DEarly")
  expect_setequal(Parabola2DEarlyFactory()$bitlength(), rep(20, 2))
  expect_setequal(Parabola2DEarlyFactory()$genelength(), 40)
  expect_setequal(Parabola2DEarlyFactory()$lb(), rep(-4.5, 2))
  expect_setequal(Parabola2DEarlyFactory()$ub(), rep(4.5, 2))
  expect_equal(Parabola2DEarlyFactory()$f(c(1,1)), 2)
  expect_output(Parabola2DEarlyFactory()$describe(), regexp="unimodal")
  expect_equal(Parabola2DEarlyFactory()$solution()$minimum, 0.0)
  expect_equal(Parabola2DEarlyFactory()$solution()$maximum, 40.5)
  expect_setequal(Parabola2DEarlyFactory()$solution()$minpoints[[1]], rep(0, 2))
  expect_setequal(Parabola2DEarlyFactory()$solution()$maxpoints[[1]], rep(4.5, 2))
  }
) 

test_that("Problem environment Parabola2DEarly (min) produces the correct results",
  {
 solution<-list()
 solution$phenotypeValue<-0.001
 lF<-list()
 lF$TerminationEps<-function() {0.01}
 lF$Max<-function() {FALSE}
 lF$penv<-Parabola2DEarlyFactory()
  expect_identical(Parabola2DEarlyFactory()$terminate(solution, lF), TRUE)
 solution$phenotypeValue<-0.1
  expect_identical(Parabola2DEarlyFactory()$terminate(solution, lF), FALSE)
  }
) 

test_that("Problem environment Parabola2DEarly (max) produces the correct results",
  {
 solution<-list()
 solution$phenotypeValue<-40.495
 lF<-list()
 lF$TerminationEps<-function() {0.01}
 lF$Max<-function() {TRUE}
 lF$penv<-Parabola2DEarlyFactory()
  expect_identical(Parabola2DEarlyFactory()$terminate(solution, lF), TRUE)
 solution$phenotypeValue<-40.0
  expect_identical(Parabola2DEarlyFactory()$terminate(solution, lF), FALSE)
  }
) 

test_that("Problem environment Parabola2DErr produces the correct results",
  {expect_identical(Parabola2DErrFactory()$name(), "Parabola2D")
  expect_setequal(Parabola2DErrFactory()$bitlength(), rep(20, 2))
  expect_setequal(Parabola2DErrFactory()$genelength(), 40)
  expect_setequal(Parabola2DErrFactory()$lb(), rep(-4.5, 2))
  expect_setequal(Parabola2DErrFactory()$ub(), rep(4.5, 2))
  a<-Parabola2DErrFactory()
  expect_error(while (TRUE) {a$f(c(1,1))})
  expect_error(while (TRUE) {a$f(c(1,1))})
  expect_error(while (TRUE) {a$f(c(1,1))})
  expect_output(Parabola2DErrFactory()$describe(), regexp="error")
  expect_equal(Parabola2DErrFactory()$solution()$minimum, 0.0)
  expect_equal(Parabola2DErrFactory()$solution()$maximum, 40.5)
  expect_setequal(Parabola2DErrFactory()$solution()$minpoints[[1]], rep(0, 2))
  expect_setequal(Parabola2DErrFactory()$solution()$maxpoints[[1]], rep(4.5, 2))
  }
) 
