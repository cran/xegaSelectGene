% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectGene.R
\name{parm}
\alias{parm}
\title{Factory for constants}
\usage{
parm(x)
}
\arguments{
\item{x}{A constant.}
}
\value{
The constant function.
}
\description{
\code{parm()} builds a constant function for \code{x}. 
            See Wickham (2019).
}
\examples{
TournamentSize<-parm(2)
TournamentSize()
Eps<-parm(0.01)
Eps()
}
\references{
Wickham, Hadley (2019): Advanced R, CRC Press, Boca Raton.
}
