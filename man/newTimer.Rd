% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timer.R
\name{newTimer}
\alias{newTimer}
\title{Timer for R code chunks.}
\usage{
newTimer()
}
\value{
\code{newTimer()} returns a timer function.

\code{a_timer_function()} returns the used time in seconds
        (invisible).

\code{a_timer_function("TimeUsed")} returns the used time in seconds.

\code{a_timer_function("Count")} returns the number of executions
         of a timed function and/or a timed block of R-Code in seconds.
}
\description{
\code{newTimer()} sets up a timer object with two 
             internal state variables, namely \code{count} 
             to count the number of timer calls and
             \code{tUsed} to calculate the total time spent in a code block
             between two timer calls.
}
\details{
\itemize{
      \item 
            Generate a timer: 
            \code{a<-newTimer()} sets up the timer \code{a()}.
            The timer \code{a()} supports three methods:
      \enumerate{
       \item \code{a()} or 
             \code{a("Measure")} or 
             \code{a(method="Measure")} 
             \strong{starts} the timer when called 1st, 3rd, 5th, ... time
             and \strong{stops} the timer 
             when called the 2nd, 4th, 6th, ... time.
             The calls can be manually inserted 
             before and after a block of R-code for profiling.
       \item \code{a("TimeUsed")} or
             \code{a(method="TimeUsed")} returns the time used in seconds.
       \item \code{a("Count")} or
             \code{a(method="Count")} returns the number of times 
             the function/block or R-code has been executed.
          }
       \item 
       The second way of usage is with the \code{Timed()} function:  
      \enumerate{
      \item Generate a timer: 
            \code{a<-newTimer()} sets up the timer \code{a()}.
       \item You convert a function \code{b()} into a timed function
             \code{bTimed()} by 
             \code{bTimed<-Timed(a, b)}. 
       \item You use \code{bTimed()} instead of \code{b()}.
       \item At the end, you can query the aggregated time and 
             the aggregated number of executions by 
             \code{a("TimeUsed")} and 
             \code{a("Count")}, respectively. 
       }
        }
}
\examples{
   a<-newTimer() 
   a(); Sys.sleep(2); a()
   a("TimeUsed")
   a("Count")
}
\seealso{
Other Performance Measurement: 
\code{\link{Counted}()},
\code{\link{Timed}()},
\code{\link{newCounter}()}
}
\concept{Performance Measurement}
