
#
# (c) 2023 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-Level Functions
#                 Independent of gene representation.
#          Package: selectGene
#
#  - Scaling functions: identity, dynamic, linear, power law, logarithmic, ... 
#

#' Scaling Fitness
#'
#' @description Fitness is transformed by a power function 
#'              \code{fit^k}.
#'              If \code{k} is 
#'              \itemize{
#'               \item less than 1: Selection pressure is decreased.
#'               \item 1:  Selection pressure remains constant.
#'               \item larger than 1: Selection pressure is increased.
#'               \item 0:  Fitness is constant. Random selection.
#'               \item smaller than 0: Fitness is 
#'              \code{1/(fit^k)}.
#'                      }
#'              
#' @details Power functions are used for contrast sharpening or softening
#'          in image analysis. 
#'          For fuzzy sets representing the value of a linguistic variable,
#'          the power function has been used as concentration or dilation
#'          transformations for modeling adverbs.
#'
#' @references Wenstop, Fred (1980) 
#'          Quantitative Analysis with Linguistic Variables.
#'          Fuzzy Sets and Systems, 4(2), pp. 99-115.
#'          <doi:10.1016/0165-0114(80)90031-7>
#' 
#' @param fit A fitness vector.
#' @param k   Scaling exponent.
#' @param lF  Local configuration.
#' 
#' @return A scaled fitness vector.
#'
#' @family Scaling
#'
#' @examples
#' lF<-list()
#' lF$Offset<-parm(0.0001)
#' fit<-sample(10, 20, replace=TRUE)
#' fit
#' ScaleFitness(fit, 0.5, lF)
#' @export
ScaleFitness<-function(fit, k, lF) 
 { minimum<-min(fit)
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
  return(fit^k)
 }

#' Abstract interface for ScaleFitness.
#'
#' @description The scaling constant \code{k} is  
#'              set by the function \code{lF$ScalingExp}.
#'
#' @param fit  A fitness vector.
#' @param lF   Local configuration.
#'
#' @return Scaled fitness vector
#'
#' @family Scaling
#'
#' @examples
#' lF<-list()
#' lF$Offset<-parm(0.0001)
#' lF$ScalingExp<-parm(2)
#' fit<-sample(10, 20, replace=TRUE)
#' fit
#' ScalingFitness(fit, lF)
#' @export
ScalingFitness<-function(fit, lF) 
{ ScaleFitness(fit, lF$ScalingExp(), lF)}

#' Dispersion Ratio Based Threshold Fitness Scaling.  
#'
#' @description Fitness is transformed by a power function 
#'              with a scaling exponent.
#'              The choice of the scaling exponent depends on 
#'              the ratio of the dispersion measures of the 
#'              current and the previous population fitness.
#'           
#' @details  The scaling exponent is selected by the following rule:
#'           \itemize{ 
#'           \item If \code{lF$RDM()>1+lF$ScalingThreshold()} 
#'                 then choose the scaling exponent \code{lF$ScalingExp()}.   
#'                 The scaling exponent should be larger than 1 to increase 
#'                 the selection pressure.
#'           \item If \code{lF$RDM()<1+lF$SCalingThreshold}
#'           and \code{lF$RDM()>1-lF$SCalingThreshold}, the fitness is not scaled.
#'           \item If \code{lF$RDM()<1-lF$SCalingThreshold}
#'                 then choose the scaling exponent \code{lF$ScalingExp2()}.   
#'                 The scaling exponent should be smaller than 1 to decrease 
#'                 the selection pressure.
#'           }              
#'
#' @param fit   Fitness vector.
#' @param lF    Local configuration.
#' 
#' @return Scaled fitness vector.
#'
#' @family Scaling
#' @family Adaptive Parameter
#'
#' @examples
#' lF<-list()
#' lF$Offset<-parm(0.0001)
#' lF$ScalingThreshold<-parm(0.05)
#' lF$RDM<-parm(1.0)
#' lF$ScalingExp<-parm(2.0)
#' lF$ScalingExp2<-parm(0.5)
#' fit<-sample(10, 20, replace=TRUE)
#' fit
#' ThresholdScaleFitness(fit, lF)
#' lF$RDM<-parm(1.2)
#' ThresholdScaleFitness(fit, lF)
#' lF$RDM<-parm(0.8)
#' ThresholdScaleFitness(fit, lF)
#' @export
ThresholdScaleFitness<-function(fit, lF) 
 { 
	 if (lF$RDM()>1+lF$ScalingThreshold())
      # increase pressure
      {return(ScaleFitness(fit, lF$ScalingExp(), lF))}
	 if (lF$RDM()<1-lF$ScalingThreshold())
      # decrease pressure
      {return(ScaleFitness(fit, lF$ScalingExp2(), lF))}
      return(fit)
}

#' Dispersion Ratio Based Continuous Fitness Scaling. 
#'
#' @description The scaling exponent is the product of 
#'              \code{lF$RDMWeight()} and \code{lF$RDM()}.
#'
#' @param fit A fitness vector.
#' @param lF  Local configuration.
#'
#' @return Scaled fitness vector.
#'
#' @family Scaling
#' @family Adaptive Parameter
#'
#' @examples
#' lF<-list()
#' lF$Offset<-parm(0.0001)
#' lF$RDMWeight<-parm(2)
#' lF$RDM<-parm(1.2)
#' fit<-sample(10, 20, replace=TRUE)
#' fit
#' ContinuousScaleFitness(fit, lF)
#' @export
ContinuousScaleFitness<-function(fit, lF) 
      {return(ScaleFitness(fit,(lF$RDMWeight()*lF$RDM()), lF))}

#' Scaling Factory
#' 
#' @param method   A scaling method. Available methods are:
#'        \itemize{ 
#'         \item "NoScaling": Identity (Default).       
#'         \item "ConstantScaling": \code{fit^k} with constant exponent.
#'               Function \code{ConstantScaling()}.
#'         \item "ThresholdScaling": 
#'         \itemize{
#'         \item 
#'         If the dispersion ratio is larger than \code{1+threshold}, 
#'         use a constant scaling exponent with a value below 1 
#'         (decrease of selection pressure).
#'               Function \code{ThresholdScaling()}.
#'         \item If the dispersion ratio is lower than \code{1-threshold}, 
#'         use a constant scaling exponent with a value above 1
#'         (increase of selection pressure).
#'         \item Else use a scaling exponent of 1. This means no scaling.
#'         }
#'         \item "ContinuousScaling": Use weighted dispersion ratio 
#'         as scaling exponent.
#'         Function \code{ContinuousScaling()}.
#'         }
#'
#' @return A scaling function. 
#'
#' @family Configuration
#'
#' @examples
#' fit<-sample(10, 20, replace=TRUE)
#' lF<-list()
#' lF$ScalingExp<-parm(2)
#' Scale<-ScalingFactory()
#' fit
#' Scale(fit, lF)
#' Scale<-ScalingFactory("ConstantScaling")
#' Scale(fit, lF)
#' @export
ScalingFactory<-function(method="NoScaling") {
if (method=="NoScaling") {f<-function(x, lF) {x}}
if (method=="ConstantScaling") {f<-ScalingFitness}
if (method=="ThresholdScaling") {f<-ThresholdScaleFitness}
if (method=="ContinuousScaling") {f<-ContinuousScaleFitness}
if (!exists("f", inherits=FALSE)) 
	{stop("Scaling label ", method, " does not exist")}
return(f)
}

#' Dispersion Ratio
#'
#' @description The dispersion ratio is computed as 
#'      the ratio \code{DM(t)/DM(k)}
#'      where \code{DM(t)} is the dispersion measure of period t and 
#'      \code{DM(k)} the dispersion measure of period \code{max(1, (t-k))}.
#'      \code{k} is specified by \code{lF$ScalingDelay()}.
#'     
#' @details The dispersion ratio may take unreasonably high and low values
#'          leading to numerical underflow or overflow 
#'          of fitness values. Therefore,
#'          we use hard thresholding to force 
#'          the dispersion ratio into the interval 
#'          \code{[lF$DRmin(), lF$DRmax()]}.
#'          The default interval is \code{[0.5, 2.0]}.
#'
#' @param popStat   Population statistics.
#' @param DM        Dispersion function. 
#' @param lF        Local configuration.
#'
#' @return Dispersion ratio. 
#'
#' @family Scaling
#'
#' @examples
#' p<-matrix(0, nrow=3, ncol=8)
#' p[1,]<-c(14.1,  0.283,  5.53, 14.0, 19.4, 38.1, 90.2, 6.54)
#' p[2,]<-c(20.7,  0.794, 14.63, 19.0, 26.5, 38.8, 71.4, 5.27)
#' p[3,]<-c(24.0,  6.007, 16.89, 24.1, 29.2, 38.8, 73.4, 6.50)
#' F<-list()
#' F$ScalingDelay<-function() {1}
#' F$DRmax<-function() {2.0}
#' F$DRmin<-function() {0.5}
#' dm<-DispersionMeasureFactory("var")
#' DispersionRatio(p, dm, F)
#' F$ScalingDelay<-function() {2}
#' DispersionRatio(p, dm, F)
#' @export
DispersionRatio<-function(popStat, DM, lF)
{
  t<-nrow(popStat)
  k<-max(1, (t-lF$ScalingDelay()))
#  cat("DM(popStat[t,]):",DM(popStat[t,]), 
#      "divided", "DM(popStat[k,]):", DM(popStat[k,]), 
#      "is", (DM(popStat[t,])/DM(popStat[k,])), "\n")
  ratio<-DM(popStat[t,])/DM(popStat[k,])
### In case of overflow/underflow, set ratio to 1.0
  if (is.na(ratio)) {ratio<-1.0} # nocov
  ##### lF$DRmin / lF$DRmax ??
  ### smoothing by tanh centered on 1.
  ### hard thresholding!
  if (ratio>lF$DRmax()) {ratio<-lF$DRmax()}
  if (ratio<lF$DRmin()) {ratio<-lF$DRmin()}
  return(ratio)
}

#' Configure dispersion measure.
#'
#' \code{DispersionMeasureFactory()} returns a function 
#' for the dispersion measure as specified by a label.
#'     If an invalid label is 
#'     specified, the configuration fails.
#'
#' @param method   A dispersion measure.
#'     \itemize{
#'     \item  "var": Variance (Default).
#'     \item "std": Standard deviation.
#'     \item "mad": Median absolute deviation (\code{mad(vec, constant=1)}).
#'     \item "cv": Coefficient of variation". 
#'     \item "range": Range.
#'     \item "iqr": Inter quartile range 
#'            (approximated by the lower and upper hinge of \code{fivenum}). 
#'     }
#'     If an invalid label is 
#'     specified, the configuration fails.
#' 
#' @return A function which computes the dispersion measure from the vector of 
#'         population statistics produced by \code{xegaObservePopulation}
#'         of package \code{xegaPopulation}.
#' 
#' @family Configuration
#'
#' @examples
#' require(stats)
#' fit<-sample(30, 20, replace=TRUE)
#' populationStats<-c(mean(fit), fivenum(fit), var(fit), mad(fit, constant=1))
#' dm<-DispersionMeasureFactory("var")
#' dm(populationStats)
#' dm<-DispersionMeasureFactory("range")
#' dm(populationStats)
#' @export 
DispersionMeasureFactory<-function(method="var")
{
if (method == "var") {f<-function(popstatvec) {popstatvec[7]}}
if (method == "std") {f<-function(popstatvec) {popstatvec[7]^0.5}}
if (method == "mad") {f<-function(popstatvec) {popstatvec[8]}}
if (method == "cv") {f<-function(popstatvec) {(popstatvec[7]^0.5)/popstatvec[1]}}
if (method == "range") {f<-function(popstatvec) {(popstatvec[6]-popstatvec[2])}}
if (method == "iqr") {f<-function(popstatvec) {(popstatvec[5]-popstatvec[3])}}
if (!exists("f", inherits=FALSE)) 
	{stop("Dispersion measure label ", method, " does not exist")}
return(f)
}

# end of file
