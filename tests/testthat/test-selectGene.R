library(testthat)
library(xegaSelectGene)

FitnessRange<-1000

test_that("NewlFselectGenes OK",
          {lF<-NewlFselectGenes()
	   expect_identical(lF$SelectionContinuation(), TRUE)
           expect_equal(lF$Max(), 1)
           expect_equal(lF$Offset(), 1)
           expect_equal(lF$Eps(), 0.01)
           expect_equal(lF$TournamentSize(), 2)
           expect_equal(lF$SelectionBias(), 1.5)
          }
)

test_that("SelectPropFitOnln OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFitOnln(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectPropFitOnln(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFitOnln((fit^2), lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectPropFit OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFit(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectPropFit(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFit(fit^2, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectPropFitM OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFitM(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectPropFitM(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFitM(fit^2, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectPropFitDiffOnln OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFitDiffOnln(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectPropFitDiffOnln(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFitDiffOnln((mean(fit)+fit)^2, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectPropFitDiff OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFitDiff(fit, lF)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   f2<-SelectPropFitDiff(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFitDiff((fit+mean(fit))^2, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectPropFitDiffM OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectPropFitDiffM(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectPropFitDiffM(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectPropFitDiffM((fit+mean(fit))^2, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectUniform OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectUniform(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectUniform(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectUniform(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit)-sqrt(var(fit)))
           expect_lt(mean(fit[f100]), mean(fit)+sqrt(var(fit)))
          }
)

test_that("SelectUniformP OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectUniformP(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectUniformP(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectUniformP(fit, lF, 100)
           expect_identical(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectDuel OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectDuel(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectDuel(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectDuel(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectTournament 2 OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectTournament(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectTournament(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectTournament(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectTournament 5 OK",
          {lF<-NewlFselectGenes()
	   lF$TournamentSize<-parm(5)
           expect_equal(lF$TournamentSize(), 5)
           fit<-sample(FitnessRange, 100, replace=TRUE)
           fit<-fit - mean(fit)
           f1<-SelectTournament(fit, lF)
           expect_identical(fit[f1]%in%fit, TRUE)
           fit<-sample(FitnessRange, 100, replace=TRUE)
           f2<-SelectTournament(fit, lF, 2)
           expect_identical(fit[f2]%in%fit, rep(TRUE,2))
           f100<-SelectTournament(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectSUS (negative fitness) OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectSUS(fit, lF)
	   expect_identical(fit[f1]%in%fit, TRUE)
	  })

test_that("SelectSUS OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   f2<-SelectSUS(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectSUS(fit^2, lF, 10)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)


test_that("SelectLRSelective OK",
          {lF<-NewlFselectGenes()
	   lF$SelectionBias<-parm(2)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectLRSelective(fit, lF)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   f2<-SelectLRSelective(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectLRSelective(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectLinearRankTSR OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   fit<-fit - mean(fit)
	   f1<-SelectLinearRankTSR(fit, lF)
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   expect_identical(fit[f1]%in%fit, TRUE)
	   f2<-SelectLinearRankTSR(fit, lF, 2)
	   expect_identical(fit[f2]%in%fit, rep(TRUE,2))
	   f100<-SelectLinearRankTSR(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("TransformSelect SelectSUS OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-TransformSelect(fit, lF, SelectSUS)
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory() OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory()
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory Uniform OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory(method="Uniform")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit)-sqrt(var(fit)))
           expect_lt(mean(fit[f100]), mean(fit)+sqrt(var(fit)))
          }
)

test_that("SelectGeneFactory UniformP OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory(method="UniformP")
	   f100<-cSelect(fit, lF, 100)
           expect_equal(mean(fit[f100]), mean(fit))
	   f110<-cSelect(fit, lF, 110)
           expect_identical((mean(fit[f110])==mean(fit)), FALSE)
          }
)

test_that("SelectGeneFactory ProportionalOnln OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory(method="ProportionalOnln")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory Proportional OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory(method="Proportional")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory ProportionalM OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("ProportionalM")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory PropFitDiffOnln OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("PropFitDiffOnln")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory PropFitDiff OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("PropFitDiff")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory PropFitDiffM OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("PropFitDiffM")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory Duel OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("Duel")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory Tournament OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("Tournament")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory STournament OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("STournament")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory LRSelective OK",
          {lF<-NewlFselectGenes()
	   lF$SelectionBias<-function() {1.99}
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("LRSelective")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory LRTSR OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("LRTSR")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory SUS OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   cSelect<-SelectGeneFactory("SUS")
	   f100<-cSelect(fit, lF, 100)
           expect_gt(mean(fit[f100]), mean(fit))
          }
)

test_that("SelectGeneFactory HUGO OK",
          {lF<-NewlFselectGenes()
	   fit<-sample(FitnessRange, 100, replace=TRUE)
	   expect_error(SelectGeneFactory("HUGO")) 
          }
)
