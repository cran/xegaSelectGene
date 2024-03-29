% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeJongF4.R
\name{DeJongF4Factory}
\alias{DeJongF4Factory}
\title{Factory for function F4 (30-dimensional quartic with noise)}
\usage{
DeJongF4Factory()
}
\value{
A problem environment represented as a list of functions:
        \itemize{
        \item \code{$name()}: The name of the problem environment.
        \item \code{$bitlength()}: The vector of the number of bits of 
                                   each parameter of the function.
        \item \code{$genelength()}: The number of bits of the gene.
        \item \code{$lb()}: The vector of lower bounds of the parameters.
        \item \code{$ub()}: The vector of upper bounds of the parameters.
        \item \code{$f(parm, gene=0, lF=0)}): The fitness function. 
        }
        Additional elements:
        \itemize{
        \item \code{$describe()}: Print a description of the problem environment to the console.
        \item \code{$solution()}: The solution structure. A named list with \code{minimum}, \code{maximum} and
                                  2 lists of equivalent solutions: \code{minpoints}, \code{maxpoints}. 
        }
}
\description{
This function factory sets up the problem environment
             for De Jong's function F4.
F4 is a 30-dimensional quartic function with Gaussian noise.
It is a continuous, convex, unimodal, high-dimensional quartic function
with Gaussian noise. For validation,  \eqn{\epsilon = 3*\sigma} 
will work most of the time.
Note: There exist \eqn{2^{30}} maxima (without noise)!
}
\examples{
DeJongF4<-DeJongF4Factory()
DeJongF4$name()
DeJongF4$bitlength()
DeJongF4$genelength()
DeJongF4$lb()
DeJongF4$ub()
DeJongF4$f(c(2.01, -1.05, 4.95, -4.3, -3.0))
DeJongF4$f(c(2.01, -1.05, 4.95, -4.3, -3.0))
DeJongF4$describe()
DeJongF4$solution()
}
\references{
De Jong, Kenneth A. (1975):
\emph{An Analysis of the Behavior of a Class of Genetic Adaptive Systems.}
PhD thesis, Michigan, Ann Arbor, pp. 203-206.
<https://deepblue.lib.umich.edu/handle/2027.42/4507>
}
\seealso{
Other Problem Environments: 
\code{\link{DelayedPFactory}()},
\code{\link{Parabola2DEarlyFactory}()},
\code{\link{Parabola2DErrFactory}()},
\code{\link{Parabola2DFactory}()},
\code{\link{envXOR}},
\code{\link{lau15}},
\code{\link{newEnvXOR}()},
\code{\link{newTSP}()}
}
\concept{Problem Environments}
