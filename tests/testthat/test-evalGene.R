library(testthat)
library(xegaSelectGene)

parm<-function(x) {function() {x}}

test_that("NewlFevalGenes OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
           expect_equal(lF$Max(), 1)
	   g<-list(gene1=c(1, 1))
	   expect_equal(sum(lF$DecodeGene(g, lF)), 2)
           expect_identical(lF$penv$name(), "Parabola2D")
          }
)

test_that("EvalGeneU OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, 0, 0))
           expect_equal(EvalGeneU(g1, lF)$fit, 0)
           expect_identical(EvalGeneU(g1, lF)$evaluated, TRUE)
          }
)

test_that("EvalGeneU Fails OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  lF$CWorstFitness<-parm(0.0)
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, "a", 0))
         g2<-EvalGeneU(g1, lF)
          expect_equal(g2$fit, lF$CWorstFitness())
          expect_identical(g2$evalFail, TRUE)
          expect_identical(g2$evaluated, TRUE)
          }
)

test_that("EvalGeneDet OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, 0, 0))
           expect_equal(EvalGeneDet(g1, lF)$fit, 0)
	  g2<-EvalGeneDet(g1, lF)
           expect_identical(g2$evaluated, TRUE)
           expect_equal(EvalGeneDet(g2, lF)$fit, 0)
          }
)

test_that("EvalGeneDet OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  lF$CWorstFitness<-parm(0.0)
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, "a", 0))
          g2<-EvalGeneDet(g1, lF)
          expect_equal(g2$fit, lF$CWorstFitness())
          expect_identical(g2$evalFail, TRUE)
          expect_identical(g2$evaluated, TRUE)
	  }
	  )


test_that("EvalGeneStoch Fail First OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  lF$CWorstFitness<-parm(0.0)
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, "a", 0))
          g2<-EvalGeneStoch(g1, lF)
           expect_identical(g2$evaluated, TRUE)
           expect_identical(g2$evalFail, TRUE)
           expect_equal(g2$obs, 0)
           expect_equal(g2$var, 0)
           expect_equal(g2$fit, lF$CWorstFitness())
          }
)

test_that("EvalGeneStoch Fail Second OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  lF$CWorstFitness<-parm(0.0)
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, 0, 0))
           expect_equal(EvalGeneStoch(g1, lF)$fit, 0)
	  g2<-EvalGeneStoch(g1, lF)
           expect_identical(g2$evaluated, TRUE)
           expect_equal(g2$obs, 1)
           expect_equal(g2$var, 0)
           expect_equal(g2$fit, 0)
	  g2$gene1=c(0, 0, "a")
	  g3<-EvalGeneStoch(g2, lF)
           expect_equal(g3$obs, 1)
           expect_equal(g3$var, 0)
           expect_equal(g3$fit, 0)
          }
)

test_that("EvalGeneFactory OK",
          {lF<-NewlFevalGenes(Parabola2DFactory())
	  g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(0, 0, 0))
	  EvalFun<-EvalGeneFactory()
           expect_equal(EvalFun(g1, lF)$fit, 0)
	  EvalFun<-EvalGeneFactory("EvalGeneU")
           expect_equal(EvalFun(g1, lF)$fit, 0)
	  EvalFun<-EvalGeneFactory("Deterministic")
           expect_equal(EvalFun(g1, lF)$fit, 0)
	  EvalFun<-EvalGeneFactory("Stochastic")
           expect_equal(EvalFun(g1, lF)$fit, 0)
	   expect_error(EvalGeneFactory("Stchastic"))
          }
)

test_that("testEvalGeneStoch OK",
          {lF<-NewlFevalGenes(DeJongF4Factory())
          set.seed(1)
          g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(1.0, -1.5))
	  g10<-testEvalGeneStoch(g1, lF, 10)
	  expect_identical(g10$evaluated, TRUE)
	  expect_identical(g10$evalFail, FALSE)
	  expect_identical(g10$obs, 10)
	  expect_gt(g10$fit, 7)
	  expect_lt(g10$fit, 15)
          expect_lt(g10$sigma, 1.5)
	  expect_gt(g10$sigma, 0.5)
	  expect_equal(g10$sigma, sqrt(g10$var/g10$obs))
          }
)


