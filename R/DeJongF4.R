
#' Factory for function F4 (30-dimensional quartic with noise)
#'
#' @description This function factory sets up the problem environment
#'              for De Jong's function F4.
#' F4 is a 30-dimensional quartic function with Gaussian noise.
#' It is a continuous, convex, unimodal, high-dimensional quartic function
#' with Gaussian noise. For validation,  \eqn{\epsilon = 3*\sigma} 
#' will work most of the time.
#' Note: There exist \eqn{2^{30}} maxima (without noise)!
#'
#' @references De Jong, Kenneth A. (1975):
#' \emph{An Analysis of the Behavior of a Class of Genetic Adaptive Systems.}
#' PhD thesis, Michigan, Ann Arbor, pp. 203-206.
#' <https://deepblue.lib.umich.edu/handle/2027.42/4507>
#'
#' @inherit DelayedPFactory return
#' 
#' @family Problem Environments
#'
#' @examples
#' DeJongF4<-DeJongF4Factory()
#' DeJongF4$name()
#' DeJongF4$bitlength()
#' DeJongF4$genelength()
#' DeJongF4$lb()
#' DeJongF4$ub()
#' DeJongF4$f(c(2.01, -1.05, 4.95, -4.3, -3.0))
#' DeJongF4$f(c(2.01, -1.05, 4.95, -4.3, -3.0))
#' DeJongF4$describe()
#' DeJongF4$solution()
#' @importFrom stats rnorm
#' @export
DeJongF4Factory<-function() {
self<-list()
self<-c(self, name=function() {"DeJongF4"})
self<-c(self, bitlength=function() {rep(64,30)})
self<-c(self, genelength=function() {sum(self$bitlength())})
self<-c(self, lb=function() {rep(-1.28,30)})
self<-c(self, ub=function() {rep(1.28,30)})
self<-c(self, f=function(parm, gene=0, lF=0) 
	    {  sum(seq(1:length(parm))*parm^{4})+stats::rnorm(1) })
#	    {  sum(seq(1:length(parm))*parm^{4})})
self<-c(self, describe=function() {
cat("See", "\n")
cat("De Jong (1975) ")
cat("An Analysis of the Behavior of a Class of Genetic Adaptive Systems.", "\n")
cat("PhD thesis, Michigan, Ann Arbor, pp. 203-206", "\n\n")
cat("F4 is a 30-dimensional quartic function with Gaussian noise.", "\n")
cat("It is a continuous, convex, unimodal, high-dimensional quartic function", "\n")
cat("with Gaussian noise. For validation:  eps=3*sigma will work most of the time.", "\n")
cat("Note: There exist 2^30 maxima (without noise)!", "\n")
                   })
self<-c(self, maxp=function()
{
	c(1.28, -1.28)[sample((1:2), 30, replace=TRUE)]
}
)
self<-c(self, solution=function() {
                 s<-list()
                 s[["minimum"]]<-0
                 s[["minpoints"]]<-list(rep(0, 30))
                 s[["maximum"]]<-1248.225
                 s[["maxpoints"]]<-list(rep(1.28,30),
				      self$maxp(),	
                                      rep(-1.28,30))
                 return(s)
                   })
return(self)
}

