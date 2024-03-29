% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectGeneBenchmark.R
\name{runSelectBenchmarks}
\alias{runSelectBenchmarks}
\title{Script for testing all selection functions}
\usage{
runSelectBenchmarks(lim = c(10, 100), both = TRUE, verbose = FALSE)
}
\arguments{
\item{lim}{Vector of population sizes.}

\item{both}{For \code{both=TRUE} the selection function 
is benchmarked with and without transformation. 
For \code{both=FALSE}, only the transformed selection functions
are benchmarked.}

\item{verbose}{Boolean. Default: \code{FALSE}. 
If \code{TRUE}, the function benchmarked and the population size 
are printed to the console.}
}
\value{
A data frame sorted in ascending order of the time of
        the last column. 
        The fastest selection methods come first.
        The first row contains the population sizes with which 
        the benchmark has been performed.
        The data frame has \code{1+length(lim)} columns:
        \itemize{
        \item "Benchmark": The name of the benchmarked selection 
                           function. A "C" after the name indicates
                           that the selection function has been 
                           transformed into a lookup function.
        \item \code{length(lim)} columns with the execution times in seconds.
        }
}
\description{
Script for testing all selection functions
}
\examples{
runSelectBenchmarks(lim=c(10, 100), both=TRUE, verbose=TRUE)
runSelectBenchmarks(lim=c(10, 100), both=FALSE)
}
\seealso{
Other Benchmark Selection Functions: 
\code{\link{predictSelectTime}()},
\code{\link{runOneBenchmark}()},
\code{\link{selectBenchmark}()},
\code{\link{testSelectGene}()}
}
\concept{Benchmark Selection Functions}
