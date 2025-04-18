% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timer.R
\name{Timed}
\alias{Timed}
\title{Transformation into a timed function}
\usage{
Timed(FUN, timer)
}
\arguments{
\item{FUN}{A function whose run time should be measured.}

\item{timer}{A timer generated by \code{newTimer()}.}
}
\value{
A timed function.
}
\description{
\code{Timed()} takes two functions as arguments, 
    namely the function whose time and call frequency 
    should be measured and a timer object created by \code{newTimer()}.
    It returns a timed function.
}
\examples{
    test<-function(seconds) {Sys.sleep(seconds)} 
    testTimer<-newTimer()
    testTimed<-Timed(test, testTimer)
    testTimer("Count"); testTimer("TimeUsed")
    testTimed(1); testTimed(2)
    testTimer("Count") 
    testTimer("TimeUsed")
}
\seealso{
Other Performance Measurement: 
\code{\link{Counted}()},
\code{\link{newCounter}()},
\code{\link{newTimer}()}
}
\concept{Performance Measurement}
