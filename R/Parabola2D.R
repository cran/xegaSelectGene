
#  Problem Environments: Parabola2D
# (c) 2023 Andreas Geyer-Schulz

#' Factory for a 2-dimensional quadratic parabola with delayed execution.
#'
#' @description This list of functions sets up the problem environment
#'              for a 2-dimensional
#'              quadratic parabola with a delayed execution of 0.1s.
#'              This function aims to test strategies
#'              of distributed/parallel execution of functions.
#'
#' @details The factory contains examples of all functions 
#'          which form the interface of a problem environment to
#'          the genetic algorithm with binary-coded genes
#'          of package \code{xega}.
#'
#' @return A problem environment represented as a list of functions:
#'         \itemize{
#'         \item \code{$name()}: The name of the problem environment.
#'         \item \code{$bitlength()}: The vector of the number of bits of 
#'                                    each parameter of the function.
#'         \item \code{$genelength()}: The number of bits of the gene.
#'         \item \code{$lb()}: The vector of lower bounds of the parameters.
#'         \item \code{$ub()}: The vector of upper bounds of the parameters.
#'         \item \code{$f(parm, gene=0, lF=0)}): The fitness function. 
#'         }
#'         Additional elements:
#'         \itemize{
#'         \item \code{$describe()}: Print a description of the problem environment to the console.
#'         \item \code{$solution()}: The solution structure. A named list with \code{minimum}, \code{maximum} and
#'                                   2 lists of equivalent solutions: \code{minpoints}, \code{maxpoints}. 
#'         }
#'
#' @family Problem Environments
#' @seealso Parabola2D
#'
#' @examples
#' DelayedP<-DelayedPFactory()
#' DelayedP$f(c(2.2, 1.0))
#' @export
DelayedPFactory<-function(){
self<-list()
self<-c(self, name=function() {"DelayedP"})
self<-c(self, bitlength=function() {c(20, 20)})
self<-c(self, genelength=function() {sum(self$bitlength())})
self<-c(self, lb=function() {c(-4.5, -4.5)})
self<-c(self, ub=function() {c(4.5, 4.5)})
self<-c(self, f=function(parm, gene=0, lF=0) {
               Sys.sleep(0.1)
               sum(parm^{2}) })
self<-c(self, describe=function() {
cat("DelayedP is a 2-dimensional parabola with spherical constant-cost contours.", "\n")
cat("It is a continuous, convex, unimodal, 2-dimensional quadratic function.", "\n")
cat("The function sleeps 0.1s before it is executed.", "\n")
})
self<-c(self, solution=function() {
                 s<-list()
                 s[["minimum"]]<-0
                 s[["minpoints"]]<-list(rep(0, 2))
                 s[["maximum"]]<-40.5
                 s[["maxpoints"]]<-list( c(4.5, 4.5),
                                      c(4.5, -4.5),
                                      c(-4.5, 4.5),
                                      c(-4.5, -4.5))
                 return(s)
                   })
return(self)
}

#' Factory for a 2-dimensional quadratic parabola.
#'
#' @description This list of functions sets up the problem environment
#'              for a 2-dimensional
#'              quadratic parabola.
#'
#' @details The factory contains examples of all functions 
#'          which form the interface of a problem environment to
#'          the simple genetic algorithm with binary-coded genes
#'          of package \code{xega}.
#'
#' @inherit DelayedPFactory return
#'
#' @family Problem Environments
#' @seealso DelayedP, Parabola2DErr
#'
#' @examples
#' Parabola2D<-Parabola2DFactory()
#' Parabola2D$f(c(2.2, 1.0))
#' @export
Parabola2DFactory<-function(){
self<-list()
self<-c(self, name=function() {"Parabola2D"})
self<-c(self, bitlength=function() {c(20, 20)})
self<-c(self, genelength=function() {sum(self$bitlength())})
self<-c(self, lb=function() {c(-4.5, -4.5)})
self<-c(self, ub=function() {c(4.5, 4.5)})
self<-c(self, f=function(parm, gene=0, lF=0) {
               sum(parm^{2}) })
self<-c(self, describe=function() {
cat("Parabola2D is a 2-dimensional parabola with spherical constant-cost contours.", "\n")
cat("It is a continuous, convex, unimodal, 2-dimensional quadratic function.", "\n")
})
self<-c(self, solution=function() {
                 s<-list()
                 s[["minimum"]]<-0
                 s[["minpoints"]]<-list(c(0, 0))
                 s[["maximum"]]<-40.5
                 s[["maxpoints"]]<-list( c(4.5, 4.5),
                                      c(4.5, -4.5),
                                      c(-4.5, 4.5),
                                      c(-4.5, -4.5))
                 return(s)
                   })
return(self)
}

#' Factory for a 2-dimensional quadratic parabola with early termination check.
#'
#' @description This list of functions sets up the problem environment
#'              for a 2-dimensional
#'              quadratic parabola.
#'
#' @details The factory contains examples of all functions 
#'          which form the interface of a problem environment to
#'          the simple genetic algorithm with binary-coded genes
#'          of package \code{xega}.
#'          This factory provides examples of a termination condition,
#'          a description function, and a solution function:
#'         \itemize{
#'         \item \code{terminate(solution)}
#'               checks for an early termination condition.
#'         \item \code{describe()}
#'               shows a description of the function.
#'         \item \code{solution()}
#'               returns a list with the 
#'               \code{minimum} and the
#'               \code{maximum} values as well as the lists
#'               \code{minpoints} and \code{maxpoints} of 
#'               the minimal and the maximal points.
#'         }
#'
#' @inherit DelayedPFactory return
#'
#' @family Problem Environments
#' @seealso DelayedP, Parabola2DErr
#'
#' @examples
#' Parabola2D<-Parabola2DEarlyFactory()
#' Parabola2D$f(c(2.2, 1.0))
#' @export
Parabola2DEarlyFactory<-function(){
self<-list()
self<-c(self, name=function() {"Parabola2DEarly"})
self<-c(self, bitlength=function() {c(20, 20)})
self<-c(self, genelength=function() {sum(self$bitlength())})
self<-c(self, lb=function() {c(-4.5, -4.5)})
self<-c(self, ub=function() {c(4.5, 4.5)})
self<-c(self, f=function(parm, gene=0, lF=0) {
               sum(parm^{2}) })
### TODO. 
### The solution is the result of bestInPopulation in xegaPopulation.R
### We stop as soon as we are within 5 percent of the optimum.
self<-c(self, terminate=function(solution, lF) {
	     if (lF$Max()==TRUE) {opt<-lF$penv$solution()$maximum}  
	     else {opt<-lF$penv$solution()$minimum}  
	     eps<-max((lF$TerminationEps()*opt), lF$TerminationEps())
	    # cat("Solution:", solution$phenotypeValue)
	    # cat(" [", (opt-eps), (opt+eps), "]\n")
	     if ((solution$phenotypeValue>(opt-eps)) &
	        (solution$phenotypeValue<(opt+eps)))
	     {return(TRUE)} else {return(FALSE)}	     
	       })
### TODO.
self<-c(self, describe=function() {
cat("Parabola2D is a 2-dimensional parabola with spherical constant-cost contours.", "\n")
cat("It is a continuous, convex, unimodal, 2-dimensional quadratic function.", "\n")
})
self<-c(self, solution=function() {
                 s<-list()
                 s[["minimum"]]<-0
                 s[["minpoints"]]<-list(c(0, 0))
                 s[["maximum"]]<-40.5
                 s[["maxpoints"]]<-list( c(4.5, 4.5),
                                      c(4.5, -4.5),
                                      c(-4.5, 4.5),
                                      c(-4.5, -4.5))
                 return(s)
                   })
return(self)
}


#' Factory for a randomly failing 2-dimensional quadratic parabola.
#'
#' @description This list of functions sets up the problem environment
#'              for a 2-dimensional
#'              quadratic parabola which produces an error 
#'              with a probability of 0.5.
#'
#' @details The factory contains examples of all functions 
#'          which form the interface of a problem environment to
#'          the simple genetic algorithm with binary-coded genes
#'          of package \code{xega}.
#'
#' @inherit DelayedPFactory return
#'
#' @family Problem Environments
#' @seealso DelayedP
#'
#' @examples
#' Parabola2DErr<-Parabola2DErrFactory()
#' @importFrom stats runif
#' @export
Parabola2DErrFactory<-function(){
self<-list()
self<-c(self, name=function() {"Parabola2D"})
self<-c(self, bitlength=function() {c(20, 20)})
self<-c(self, genelength=function() {sum(self$bitlength())})
self<-c(self, lb=function() {c(-4.5, -4.5)})
self<-c(self, ub=function() {c(4.5, 4.5)})
self<-c(self, f=function(parm, gene=0, lF=0) {
	       if (0.5>runif(1)) {"a"+3}	
               sum(parm^{2}) 
		   })
self<-c(self, describe=function() {
cat("Parabola2D is a 2-dimensional parabola with spherical constant-cost contours.", "\n")
cat("It is a continuous, convex, unimodal, 2-dimensional quadratic function.", "\n")
cat("It produces an error 50 percent of the time.", "\n")
})
self<-c(self, solution=function() {
                 s<-list()
                 s[["minimum"]]<-0
                 s[["minpoints"]]<-list(c(0, 0))
                 s[["maximum"]]<-40.5
                 s[["maxpoints"]]<-list( c(4.5, 4.5),
                                      c(4.5, -4.5),
                                      c(-4.5, 4.5),
                                      c(-4.5, -4.5))
                 return(s)
                   })
return(self)
}

