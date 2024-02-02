library(testthat)
library(xegaSelectGene)

test_that("newCounter and counter OK",
          {a<-newCounter()
	   a(); a(); a()
           expect_equal(a("Show"), 3)
          }
)

test_that("newTimer and timer OK",
          {a<-newTimer()
	   a(); Sys.sleep(1.0); a()
	   expect_lt(a("TimeUsed"), 2.0)
	   expect_gt(a("TimeUsed"), 1.0)
           expect_equal(a("Count"), 1)
	   a(); Sys.sleep(1.0); a()
	   expect_lt(a("TimeUsed"), 3.0)
	   expect_gt(a("TimeUsed"), 2.0)
           expect_equal(a("Count"), 2)
          }
)

test_that("counted Function OK",
          {a<-function(s) {Sys.sleep(s)}
           t<-newCounter()
	   b<-Counted(a, t) 
           expect_equal(t("Show"), 0)
	   b(1); b(1.5)
           expect_equal(t("Show"), 2)
          }
)

test_that("timed Function OK",
          {a<-function(s) {Sys.sleep(s)}
           t<-newTimer()
	   b<-Timed(a, t) 
           expect_lt(t("TimeUsed"), 1.0)
           expect_equal(t("Count"), 0)
	   b(1); b(1.5)
           expect_lt(t("TimeUsed"), 3.0)
           expect_gt(t("TimeUsed"), 2.0)
           expect_equal(t("Count"), 2)
          }
)


