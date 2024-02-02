#
# (c) 2023 Andreas Geyer-Schulz
#          Simple Genetic Algorithm xegaX packages in R. V 0.1
#          Layer: General Profiling functions.
#

#' Counter
#' 
#' @description \code{newCounter} sets up a counter object with one 
#'              internal state variable, namely \code{count} 
#'              to count the number of counter calls.
#' 
#' @details 
#'             Generate a counter: 
#'             \code{a<-newCounter()} sets up the counter \code{a}.
#'             The counter \code{a} supports three methods:
#'       \enumerate{
#'        \item \code{a()} or 
#'              \code{a("Measure")} or 
#'              \code{a(method="Measure")} 
#'              \strong{starts} the timer when called 1st, 3rd, 5th, ... time
#'              and \strong{stops} the timer 
#'              when called the 2nd, 4th, 6th, ... time.
#'              The calls can be manually inserted 
#'              before and after a block of R-code for profiling.
#'        \item \code{a("Count")} or
#'              \code{a(method="Count")} returns the number of times 
#'              the function/block or R-code has been executed.
#'           }
#'
#' @return \code{newCounter()} returns a counter function.
#' @return \code{a_counter_function()} returns the 
#'          number of times it has been called
#'         (invisible).
#' @return \code{a_counter_function("Show")} returns the number of executions
#'         of the \code{a_counter_function}.
#'
#' @family Performance Measurement
#'
#' @examples 
#'    a<-newCounter() 
#'    a(); a()
#'    a("Show")
#' @export
newCounter<-function()
{ count<-0
  Counter<-function(method="Count")
  { if (method=="Count")
	{   count<<-count+1
	    invisible<-count }
         if (method=="Show")
         { count }
  }
  return(Counter)
}

#' Timer for R code chunks.
#' 
#' @description \code{newTimer} sets up a timer object with two 
#'              internal state variables, namely \code{count} 
#'              to count the number of timer calls and
#'              \code{tUsed} to calculate the total time spent in a code block
#'              between two timer calls.
#' 
#' @details 
#'       \itemize{
#'       \item 
#'             Generate a timer: 
#'             \code{a<-newTimer()} sets up the timer \code{a}.
#'             The timer \code{a} supports three methods:
#'       \enumerate{
#'        \item \code{a()} or 
#'              \code{a("Measure")} or 
#'              \code{a(method="Measure")} 
#'              \strong{starts} the timer when called 1st, 3rd, 5th, ... time
#'              and \strong{stops} the timer 
#'              when called the 2nd, 4th, 6th, ... time.
#'              The calls can be manually inserted 
#'              before and after a block of R-code for profiling.
#'        \item \code{a("TimeUsed")} or
#'              \code{a(method="TimeUsed")} returns the time used in seconds.
#'        \item \code{a("Count")} or
#'              \code{a(method="Count")} returns the number of times 
#'              the function/block or R-code has been executed.
#'           }
#'        \item 
#'        The second way of usage is with the \code{Timed} function:  
#'       \enumerate{
#'       \item Generate a timer: 
#'             \code{a<-newTimer()} sets up the timer \code{a}.
#'        \item You convert a function \code{b} into a timed function
#'              \code{bTimed} by 
#'              \code{bTimed<-Timed(a, b)}. 
#'        \item You use \code{bTimed} instead of \code{b}.
#'        \item At the end, you can query the aggregated time and 
#'              the aggregated number of executions by 
#'              \code{a("TimeUsed")} and 
#'              \code{a("Count")}, respectively. 
#'        }
#'         }
#'
#' @return \code{newTimer()} returns a timer function.
#' @return \code{a_timer_function()} returns the used time in seconds
#'         (invisible).
#' @return \code{a_timer_function("TimeUsed")} returns the used time in seconds.
#' @return \code{a_timer_function("Count")} returns the number of executions
#'          of a timed function and/or a timed block of R-Code in seconds.
#'
#' @family Performance Measurement
#'
#' @examples 
#'    a<-newTimer() 
#'    a(); Sys.sleep(2); a()
#'    a("TimeUsed")
#'    a("Count")
#' @export
newTimer<-function()
{ tUsed<-0
  tStart<-0
  count<-0
  Timer<-function(method="Measure")
  {
	if (method=="Measure")
	{	
          if (count%%2)
          { tEnd<-Sys.time()
            tUsed<<-tUsed+as.numeric(tEnd-tStart, units="secs")
            count<<-count+1
            invisible(tUsed)
          }
          else
          { tStart<<-Sys.time()
            count<<-count+1
            invisible(tUsed)
          }
	 }
         if (method=="TimeUsed")
         { return(tUsed) }
         if (method=="Count")
         { return(count/2) }

  }
  return(Timer)
}

#' Transformation into a counted function
#'
#' @description
#'     \code{Counted} takes two functions as arguments: 
#'     The function whose call frequency 
#'     should be measured and a counter object created by \code{newCounter()}.
#'     It returns a counted function. 
#'
#' @param FUN    A function whose run time should be measured.
#' @param counter A counter generated by \code{newCounter()}.  
#' 
#' @return A counted function.
#'
#' @family Performance Measurement
#' 
#' @examples
#'     test<-function(v) {sum(v)} 
#'     testCounter<-newCounter()
#'     testCounted<-Counted(test, testCounter)
#'     testCounter("Show")
#'     testCounted(sample(10,10)); testCounted(sample(10,10))
#'     testCounter("Show")
#' @export
Counted<-function(FUN, counter)
{
  cFUN<-function(...){z<-FUN(...);counter();return(z)}
  return(cFUN)
}

#' Transformation into a timed function
#'
#' @description
#'     \code{Timed} takes two functions as arguments, 
#'     namely the function whose time and call frequency 
#'     should be measured and a timer object created by \code{newTimer()}.
#'     It returns a timed function. 
#'
#' @param FUN    A function whose run time should be measured.
#' @param timer  A timer generated by \code{newTimer()}.  
#' 
#' @return A timed function.
#'
#' @family Performance Measurement
#' 
#' @examples
#'     test<-function(seconds) {Sys.sleep(seconds)} 
#'     testTimer<-newTimer()
#'     testTimed<-Timed(test, testTimer)
#'     testTimer("Count"); testTimer("TimeUsed")
#'     testTimed(1); testTimed(2)
#'     testTimer("Count") 
#'     testTimer("TimeUsed")
#' @export
Timed<-function(FUN, timer)
{
  tFUN<-function(...){timer();z<-FUN(...);timer();return(z)}
  return(tFUN)
}
