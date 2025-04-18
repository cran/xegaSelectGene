% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectGene.R
\name{SelectGeneFactory}
\alias{SelectGeneFactory}
\title{Configure the selection function of a genetic algorithm.}
\usage{
SelectGeneFactory(method = "PropFitDiffOnln")
}
\arguments{
\item{method}{A string specifying the selection function.}
}
\value{
A selection function for genes.
}
\description{
\code{SelectGeneFactory()} implements selection
             of one of the gene selection functions in this 
             package by specifying a text string.
             The selection fails ungracefully (produces
             a runtime error), if the label does not match.
             The functions are specified locally.             

             Current support:

             \enumerate{
             \item "Uniform" returns \code{SelectUniform()}.
             \item "UniformP" returns \code{SelectUniformP()}.
             \item "ProportionalOnln" returns \code{SelectPropFitOnln()}.
             \item "Proportional" returns \code{SelectPropFit()}.
             \item "ProportionalM" returns \code{SelectPropFitM()}.
             \item "PropFitDiffOnln" returns \code{SelectPropFitDiffOnln()}.
             \item "PropFitDiff" returns \code{SelectPropFitDiff()}.
             \item "PropFitDiffM" returns \code{SelectPropFitDiffM()}.
             \item "Tournament" returns \code{SelectTournament()}.
             \item "STournament" returns \code{SelectSTournament()}.
             \item "Duel" returns \code{SelectDuel()}.
             \item "LRSelective" returns \code{SelectLRSelective()}.
             \item "LRTSR" returns \code{SelectLinearRankTSR()}.
             \item "SUS" returns \code{SelectSUS()}.
             }
}
\details{
If \code{SelectionContinuation()==TRUE} then:
             \enumerate{
               \item In package xegaPopulation in
                     function \code{NextPopulation()},
                     first the functions 
                       \code{SelectGene()} and \code{SelectMate()} 
                       are transformed by \code{TransformSelect()} to
                       a continuation function with
                     embedded index vector and counter.
               \item For each call in \code{ReplicateGene()},
                     \code{SelectGene} and \code{SelectMate()} 
                     return the index of the selected gene.
                     }
}
\examples{
SelectGene<-SelectGeneFactory("Uniform")
fit<-sample(10, 15, replace=TRUE)
SelectGene(fit, lFselectGenes) 
sel<-"Proportional"
SelectGene<-SelectGeneFactory(method=sel)
fit<-sample(10, 15, replace=TRUE)
SelectGene(fit, lFselectGenes) 
}
\seealso{
Other Configuration: 
\code{\link{DispersionMeasureFactory}()},
\code{\link{EvalGeneFactory}()},
\code{\link{ScalingFactory}()}
}
\concept{Configuration}
