
#
# (c) 2023 Andreas Geyer-Schulz
#          selectGeneBenchmark.R:
#          Benchmark of 
#          Package: selectGene

#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-Level Functions
#                 Independent of gene representation.
#          Package: selectGene
#

#' Test a gene selection function 
#'
#' @description \code{testSelectGene()} implements 
#'              testing a selection function.
#'              It collects the results of the repeated execution 
#'              of the selection function
#'              given a fitness function.
#' 
#' @param fit             Fitness vector.
#' @param method          String specifying the selection function. 
#'                        See \code{SelectGeneFactory()}.
#' @param howOften        Integer. 
#' @param lF              Local configuration. 
#' @param continuation    Convert to index function? 
#' @param verbose         Boolean. Default: \code{FALSE}. 
#'                        If \code{TRUE}, the exection time of the transformation of the selection function
#'                        into a quasi-continuation function is printed on the console. 
#'
#' @return 
#' \itemize{
#'         \item \code{$fit} fitness vector.
#'         \item \code{$newPop} indices of survivors in fitness vector.
#'         \item \code{$time} time in seconds. 
#'         \item \code{$size} population size.
#'         \item \code{$method} selection method used.
#'         }
#'
#' @family Benchmark Selection Functions
#'
#' @examples
#' fit1<-rep(10,10)
#' fit2<-fit1+runif(rep(10,1))
#' fit3<-sample(100, 10, replace=TRUE)
#' testSelectGene(fit2, method="Tournament", howOften=100) 
#' testSelectGene(fit3, method="Tournament", howOften=10) 
#' @export
testSelectGene<-function(fit, 
			 method="Uniform", 
			 howOften=100,
			 lF=NewlFselectGenes(), continuation=TRUE, verbose=FALSE) 
{
SelectGene<-SelectGeneFactory(method=method)
v<-rep(0,howOften)

gc()

selectTimer<-newTimer()

if (continuation)
{ 

selectTimer()
	SelectGene<-TransformSelect(fit, lF, SelectGene)
selectTimer()
if (verbose) cat("TransformSelect:", selectTimer("TimeUsed"), "\n")
}

selectTimer()
for (i in 1: howOften)
{
    v[i]<-SelectGene(fit, lF)
}
selectTimer()

return(list(fit=fit, 
	    newPop=v, 
	    time=selectTimer("TimeUsed"), 
	    size=howOften, 
	    method=method, 
	    continuation=continuation))

#count<-rowSums(outer((1:length(fit)), v, FUN="=="))
#freq<-count/size
#
#t<-t(matrix(c(fit, freq, count), nrow=3, byrow=TRUE))
#
#return(t)
}

#' Benchmark and stress test of selection functions. 
#'
#' @description Times a selection function 
#'              for populations of size 10 to \eqn{10^{limit}}.
#'
#' @param method          Selection function. Default: Uniform.
#' @param continuation    Convert to index function? Default: TRUE.
#' @param limit           Vector of population sizes.
#' @param verbose         Boolean. Default: \code{FALSE}. 
#'                        If \code{TRUE}, the function benchmarked and the population size 
#'                        are printed to the console.
#'
#' @return Vector of execution times in seconds.
#'
#' @family Benchmark Selection Functions
#' @examples
#' selectBenchmark(method="Uniform", continuation=TRUE, limit=c(10, 100, 1000))
#' selectBenchmark(method="SUS", continuation=TRUE, limit=c(5000, 10000, 15000))
#' selectBenchmark(method="SUS", continuation=FALSE, limit=seq(from=100, to=1000, length.out=5))
#' @export
selectBenchmark<-function(method="Uniform", continuation=TRUE, limit=c(10, 100, 1000), verbose=FALSE)
{
	result<-vector("double", length(limit))
	for (i in (1:length(limit)))
	{     
		size<-limit[i]
		if (verbose) cat(method,"Selecting n=",size, "candidates\n")
		fit<-sample(1000, size, replace=TRUE)  
		a<-try(testSelectGene(fit, 
				      method=method, 
				      howOften=size, 
				      continuation=continuation)) 
		result[i]<-a$time
	}
	return(result)
	#unlist(result[sapply(result, function(x) !inherits(x, "try-error"))])
}

#' Script for testing a single selection functions
#'
#' @param name  Name is one of the following strings.
#'
#'              \enumerate{
#'              \item "Uniform" benchmarks \code{SelectUniform}.
#'              \item "ProportionalOnln" benchmarks \code{SelectPropFitOnln}.
#'              \item "Proportional" benchmarks \code{SelectPropFit}.
#'              \item "ProportionalM" benchmarks \code{SelectPropFitM}.
#'              \item "PropFitDiffOnln" benchmarks \code{SelectPropFitDiffOnln}.
#'              \item "PropFitDiff" benchmarks \code{SelectPropFitDiff}.
#'              \item "PropFitDiffM" benchmarks \code{SelectPropFitDiffM}.
#'              \item "Tournament" benchmarks \code{SelectTournament}.
#'              \item "Duel" benchmarks \code{SelectDuel}.
#'              \item "LinearRank" benchmarks \code{SelectLinearRank}.
#'              \item "SUS" benchmarks \code{SelectSUS}.
#'              }
#'
#' @param limit      Vector of population sizes.
#' @param both       For \code{both=TRUE} the selection function 
#'                   is benchmarked with and without transformation. 
#'                   For \code{both=FALSE}, only the transformed selection functions
#'                   are benchmarked.
#' @param verbose         Boolean. Default: \code{FALSE}. 
#'                        If \code{TRUE}, the function benchmarked and the population size 
#'                        are printed to the console.
#' 
#' @section Warning:
#'    The time to run the function for \code{lim>6} explodes 
#'    for all benchmark functions with higher than linear complexity.
#'    (e.g. \code{PropFit()}, \code{PropFitdiff()}, and \code{Tournament()}).
#' 
#' @return A data frame sorted in ascending order of time of last column. 
#'  
#' @family Benchmark Selection Functions
#'
#' @examples
#' runOneBenchmark("Duel", 5, both=FALSE)
#' runOneBenchmark("PropFitDiffOnln")
#' @export
runOneBenchmark<-function(name, limit=c(10, 100, 1000), both=TRUE, verbose=FALSE)
{
d<-data.frame()
n<-data.frame()
n<-rbind(n, paste(name," C"))
d<-rbind(d, selectBenchmark(method=name, continuation=TRUE, limit=limit, verbose=verbose))
if (both)
{ n<-rbind(n, name)
d<-rbind(d, selectBenchmark(method=name, continuation=FALSE, limit=limit, verbose=verbose))}
df<-cbind(n,d)
names(df)<-c("Benchmark", unlist(lapply(limit, toString)))
return(df)
}

#' Script for testing all selection functions
#'
#'
#'
#' @param lim      Vector of population sizes. 
#' @param both     For \code{both=TRUE} the selection function 
#'                 is benchmarked with and without transformation. 
#'                 For \code{both=FALSE}, only the transformed selection functions
#'                 are benchmarked.
#' @param verbose         Boolean. Default: \code{FALSE}. 
#'                        If \code{TRUE}, the function benchmarked and the population size 
#'                        are printed to the console.
#' 
#'
#' @return A data frame sorted in ascending order of the time of
#'         the last column. 
#'         The fastest selection methods come first.
#'         The first row contains the population sizes with which 
#'         the benchmark has been performed.
#'         The data frame has \code{1+length(lim)} columns:
#'         \itemize{
#'         \item "Benchmark": The name of the benchmarked selection 
#'                            function. A "C" after the name indicates
#'                            that the selection function has been 
#'                            transformed into a lookup function.
#'         \item \code{length(lim)} columns with the execution times in seconds.
#'         }
#'
#' @family Benchmark Selection Functions
#'
#' @examples
#' runSelectBenchmarks(lim=c(10, 100), both=TRUE, verbose=TRUE)
#' runSelectBenchmarks(lim=c(10, 100), both=FALSE)
#' @export
runSelectBenchmarks<-function(lim=c(10, 100), both=TRUE, verbose=FALSE)
{
df<-data.frame()
try(df<-rbind(df,runOneBenchmark("Uniform", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("Proportional", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("ProportionalOnln", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("ProportionalM", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("PropFitDiffOnln", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("PropFitDiff", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("PropFitDiffM", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("SUS", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("LRTSR", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("LRSelective", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("Duel", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("Tournament", lim, both, verbose)))
try(df<-rbind(df,runOneBenchmark("STournament", lim, both, verbose)))
# record with vector of  population sizes (independent variable)
n<-data.frame()
n<-rbind(n, "popsize")
d<-data.frame()
d<-rbind(d, lim)
pf<-cbind(n,d)
names(pf)<-c("Benchmark", unlist(lapply(lim, toString)))
# bind to matrix of timings
ndf<-df[order(df[,length(lim)+1]),]
ndf<-rbind(pf, ndf)
return(ndf)
}

#### TODO!

#' Predict the time use of a selection method for a popsize.
#'
#' @section Warning:
#'
#' Uses a quadratic regression model.
#' But the complexities of the functions are of orders
#'     O(1), O(n), O(n.ln(n)) and O(n^2). 
#'
#' @param df       Data frame. 
#' @param method   Selection method.
#' @param popsize  Population size.
#'
#' @return List with
#' \itemize{
#' \item
#'              \code{$model}:  The result of \code{stats::lm}.
#' \item
#'              \code{$predict}:  The result of \code{stats::predict}.
#' }
#'
#' @family Benchmark Selection Functions
#'
#'@examples
#' popsizes<-as.integer(seq(from=100, to=200, length.out=5))
#' a<-runSelectBenchmarks(popsizes, both=TRUE)
#' b<-predictSelectTime(a, method="SUS", 155)
#' summary(b$model)
#' b$predicted
#' c<- predictSelectTime(a, method="SUS  C", c(155, 500))
#' summary(c$model)
#' c$predicted
#'@importFrom stats lm predict
#'@export 
predictSelectTime<-function(df, method="Uniform", popsize=100000)
{
 l<-ncol(df)
  ix<-unlist(lapply(df[,1], function(x) {is.element(x, "popsize")}))
  iy<-unlist(lapply(df[,1], function(x) {is.element(x, method)}))
 y<-as.vector(unlist(df[iy,2:l]))
 x<-as.vector(unlist(df[ix,2:l]))
 dt<-data.frame(matrix(c(y, x,(x * log(x)) , (x^2)), ncol=4))
 names(dt)<-c("Y", "X1", "X2", "X3")
 model<-stats::lm(dt$Y~dt$X1+dt$X2+dt$X3)
 dt<-data.frame(X1=popsize, X2=popsize*log(popsize), X3=popsize^2)
 pred<-stats::predict(object=model, newdata=dt, interval="prediction")
 return(list(model=model, predicted=pred))
}	

# end of file
