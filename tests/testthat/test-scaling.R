library(testthat)
library(xegaSelectGene)

parm<-function(x) {function() {x}}

test_that("ScaleFitness OK",
{lF<-list()
   lF$Offset<-function() {0.01}
   fit<-sample(10, 15, replace=TRUE)
   expect_setequal(ScaleFitness(fit, 0, lF), rep(1, 15))
   expect_setequal(ScaleFitness(fit, 1, lF), fit)
   expect_setequal(ScaleFitness(fit, 0.5, lF), (fit^0.5))
   fitneg<-fit-10
   expect_setequal(ScaleFitness(fitneg, 0, lF), rep(1, 15))
   expect_setequal(ScaleFitness(fitneg, 1, lF), (fitneg+0.01+abs(min(fitneg))))
   expect_setequal(ScaleFitness(fitneg, 2, lF), (fitneg+0.01+abs(min(fitneg)))^2)
}
)

test_that("ScalingFitness OK",
{lF<-list()
   lF$Offset<-function() {0.01}
   lF$ScalingExp<-function() {2}
   fit<-sample(10, 15, replace=TRUE)
   expect_setequal(ScalingFitness(fit, lF), (fit^2))
   fitneg<-fit-10
   expect_setequal(ScalingFitness(fitneg, lF), (fitneg+0.01+abs(min(fitneg)))^2)
}
)

test_that("ThresholdScaleFitness OK",
{
lF<-list()
lF$Offset<-parm(0.0001)
lF$ScalingThreshold<-parm(0.05)
lF$RDM<-parm(1.0)
lF$ScalingExp<-parm(0.5)
lF$ScalingExp2<-parm(2)
fit<-sample(10, 15, replace=TRUE)
expect_setequal(ThresholdScaleFitness(fit, lF), fit)
# increase pressure
lF$RDM<-parm(1.2)
expect_setequal(ThresholdScaleFitness(fit, lF), fit^0.5)
# decrease pressure
lF$RDM<-parm(0.5)
expect_setequal(ThresholdScaleFitness(fit, lF), fit^2)
}
)

test_that("ContinuousScaleFitness OK",
{
lF<-list()
lF$Offset<-parm(0.0001)
lF$ScalingThreshold<-parm(0.05)
lF$RDM<-parm(1.0)
lF$RDMWeight<-parm(1.0)
fit<-sample(10, 15, replace=TRUE)
expect_setequal(ContinuousScaleFitness(fit, lF), fit)
# increase pressure
lF$RDM<-parm(1.2)
lF$RDMWeight<-parm(1.2)
expect_setequal(ContinuousScaleFitness(fit, lF), fit^(1.2*1.2))
# decrease pressure
lF$RDM<-parm(0.5)
lF$RDMWeight<-parm(0.9)
expect_setequal(ContinuousScaleFitness(fit, lF), fit^(0.5*0.9))
}
)

test_that("ScalingFactory OK",
          {
lF<-list()
lF$Offset<-parm(0.0001)
lF$ScalingThreshold<-parm(0.05)
lF$ScalingExp<-parm(0.5)
lF$ScalingExp2<-parm(2)
lF$RDM<-parm(1.0)
lF$RDMWeight<-parm(1.0)
fit<-sample(10, 15, replace=TRUE)
ScalingFun<-ScalingFactory()
expect_equal(ScalingFun(fit, lF), fit)
ScalingFun<-ScalingFactory("NoScaling")
expect_equal(ScalingFun(fit, lF), fit)
ScalingFun<-ScalingFactory(method="ConstantScaling")
expect_equal(ScalingFun(fit, lF), ScalingFitness(fit, lF))
ScalingFun<-ScalingFactory(method="ThresholdScaling")
expect_equal(ScalingFun(fit, lF), ThresholdScaleFitness(fit, lF))
ScalingFun<-ScalingFactory(method="ContinuousScaling")
expect_equal(ScalingFun(fit, lF), ContinuousScaleFitness(fit, lF))
expect_error(ScalingFactory(method="CntinuousScaling"))
          }
)

test_that("DispersionRatio OK",
{ lF<-NewlFselectGenes()
p<-matrix(0, nrow=3, ncol=8)
p[1,]<-c(14.1,  0.283,  5.53, 14.0, 19.4, 38.1, 90.2, 6.54)
p[2,]<-c(20.7,  0.794, 14.63, 19.0, 26.5, 38.8, 71.4, 5.27)
p[3,]<-c(24.0,  6.007, 16.89, 24.1, 29.2, 38.8, 73.4, 6.50)
dm<-DispersionMeasureFactory("var")
expect_equal(DispersionRatio(p, dm, lF), 1.0280112)
lF$ScalingDelay<-parm(2)
expect_equal(DispersionRatio(p, dm, lF), 0.81374723)
lF$DRmin<-parm(0.9)
expect_equal(DispersionRatio(p, dm, lF), 0.9)
lF$DRmin<-parm(0.2)
lF$DRmax<-parm(0.7)
expect_equal(DispersionRatio(p, dm, lF), 0.7)
          }
)

test_that("DispersionMeasureFactory OK",
          {
fit<-sample(10, 15, replace=TRUE)
populationStats<-c(mean(fit), fivenum(fit), var(fit), mad(fit, constant=1))
DM<-DispersionMeasureFactory()
expect_equal(DM(populationStats), var(fit))
DM<-DispersionMeasureFactory(method="var")
expect_equal(DM(populationStats), var(fit))
DM<-DispersionMeasureFactory(method="std")
expect_equal(DM(populationStats), var(fit)^0.5)
DM<-DispersionMeasureFactory(method="mad")
expect_equal(DM(populationStats), mad(fit, constant=1))
DM<-DispersionMeasureFactory(method="cv")
expect_equal(DM(populationStats), (var(fit)^0.5)/mean(fit))
DM<-DispersionMeasureFactory(method="range")
expect_equal(DM(populationStats), (max(fit)-min(fit)))
DM<-DispersionMeasureFactory(method="iqr")
expect_equal(DM(populationStats), (fivenum(fit)[4]- fivenum(fit)[2]))
expect_error(DispersionMeasureFactory(method="CntinuousScaling"))
          }
)


