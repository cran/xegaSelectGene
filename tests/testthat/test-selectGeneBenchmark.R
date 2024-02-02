library(testthat)
library(xegaSelectGene)

test_that("runSelectBenchmark OK",
          {lF<-NewlFselectGenes()
	   df<-runSelectBenchmarks(lim=c(10, 100), both=TRUE) 
           expect_equal(dim(df), c(27, 3))
          }

)

test_that("predictSelectTime OK",
          {lF<-NewlFselectGenes()
	   popsizes<-as.integer(seq(from=100, to=1000, length.out=6))
	   a<-runSelectBenchmarks(popsizes, both=TRUE) 
	   b<-predictSelectTime(a, method="SUS", 155)
           expect_equal(length(b), 2)
	   c<-predictSelectTime(a, method="SUS", c(155,2000))
           expect_equal(length(c), 2)
          }
)

