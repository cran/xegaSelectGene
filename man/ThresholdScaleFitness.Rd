% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{ThresholdScaleFitness}
\alias{ThresholdScaleFitness}
\title{Dispersion Ratio Based Fitness Scaling.}
\usage{
ThresholdScaleFitness(fit, lF)
}
\arguments{
\item{fit}{Fitness vector.}

\item{lF}{Local configuration.}
}
\value{
Scaled fitness vector.
}
\description{
Fitness is transformed by a power function 
             \code{fit^lF$ScalingExp}.
             If \code{lF$ScalingExp} is 
             \itemize{
              \item less than 1: Selection pressure is decreased.
              \item 1:  Selection pressure remains constant.
              \item larger than 1: Selection pressure is increased.
              \item 0:  Fitness is constant. Random selection.
              \item smaller than 1: Fitness is 
             \code{1/fit^lF$ScalingExp}.
                     }
}
\examples{
lF<-list()
lF$Offset<-parm(0.0001)
lF$ScalingThreshold<-parm(0.05)
lF$RDM<-parm(1.0)
lF$ScalingExp<-parm(0.5)
lF$ScalingExp2<-parm(2)
fit<-sample(10, 20, replace=TRUE)
fit
ThresholdScaleFitness(fit, lF)
lF$RDM<-parm(1.2)
ThresholdScaleFitness(fit, lF)
lF$RDM<-parm(0.8)
ThresholdScaleFitness(fit, lF)
}
\seealso{
Other Scaling: 
\code{\link{ContinuousScaleFitness}()},
\code{\link{DispersionRatio}()},
\code{\link{ScaleFitness}()},
\code{\link{ScalingFitness}()}

Other Adaptive Parameter: 
\code{\link{ContinuousScaleFitness}()}
}
\concept{Adaptive Parameter}
\concept{Scaling}
