% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectGene.R
\name{SelectDuel}
\alias{SelectDuel}
\title{Deterministic duel.}
\usage{
SelectDuel(fit, lF, size = 1)
}
\arguments{
\item{fit}{Fitness vector.}

\item{lF}{Local configuration.}

\item{size}{Number of selected genes. Default: 1.}
}
\value{
The index vector of the selected genes.
}
\description{
\code{SelectDuel()} implements selection
             by a tournament between 2
             randomly selected genes. The best gene always wins.
             This is the version of tournament selection 
             with the least selection pressure.
}
\details{
This is an O(n) implementation of tournament selection with 
         a tournament size of 2.

A special case of tournament selection.
}
\examples{
fit<-sample(10, 15, replace=TRUE)
SelectDuel(fit, NewlFselectGenes()) 
SelectDuel(fit, NewlFselectGenes(), length(fit)) 
}
\seealso{
Other Selection Functions: 
\code{\link{SelectLRSelective}()},
\code{\link{SelectLinearRankTSR}()},
\code{\link{SelectPropFit}()},
\code{\link{SelectPropFitDiff}()},
\code{\link{SelectPropFitDiffM}()},
\code{\link{SelectPropFitDiffOnln}()},
\code{\link{SelectPropFitM}()},
\code{\link{SelectPropFitOnln}()},
\code{\link{SelectSTournament}()},
\code{\link{SelectSUS}()},
\code{\link{SelectTournament}()},
\code{\link{SelectUniform}()},
\code{\link{SelectUniformP}()}
}
\concept{Selection Functions}
