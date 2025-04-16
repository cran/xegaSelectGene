
#
# (c) 2023 Andreas Geyer-Schulz
#          Simple Genetic Algorithm in R. V 0.1
#          Layer: Gene-Level Functions
#                 Independent of gene representation.
#          Package: selectGene
#

#
# Part I.1. Function Factories for Constants.
# Part I.2. Object with local constant functions.
#

#' Factory for constants
#' 
#'@description \code{parm()} builds a constant function for \code{x}. 
#'             See Wickham (2019).
#' 
#'@references  Wickham, Hadley (2019): Advanced R, CRC Press, Boca Raton.
#'
#'@param x     A constant.
#' 
#'@return The constant function. 
#' 
#'@examples
#' TournamentSize<-parm(2)
#' TournamentSize()
#' Eps<-parm(0.01)
#' Eps()
#'@export
parm<-function(x){function() {return(x)}}

#' Generate local functions and objects.
#' 
#'@description \code{NewlFselectGenes()} returns 
#'              the list of functions which contains 
#'              a definition of all local objects required for the use
#'              of selection functions. We reference this object 
#'              as local configuration. When adding additional 
#'              selection functions, this must be extended
#'              by the constant (functions) needed to configure them.
#'
#'@return Local configuration \code{lF}.
#'
#'@examples
#' lF<-NewlFselectGenes()
#' lF$Max()
#'@export
NewlFselectGenes<-function() 
{
	 list(
	 SelectionContinuation = parm(TRUE),
	 Max=parm(1),
         Offset = parm(1),
         Eps = parm(0.01), 
	 TournamentSize = parm(2),
	 SelectionBias = parm(1.5),
	 MaxTSR = parm(1.9),
	 DRmax = parm(2.0),
	 DRmin = parm(0.5),
         ScalingDelay = parm(1)
	 ##
         )
}

#
# Part 2. Selection functions.  
#

#' Selection proportional to fitness O(n ln(n)).
#'
#' @description \code{SelectPropFitOnln()} implements selection
#'              proportional to fitness. Negative fitness
#'              vectors are shifted to \eqn{R^+}.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              Holland's schema theorem uses this selection function.
#'              See John Holland (1975) for further information.
#'
#'  @details This is a fast implementation with equivalent results 
#'           to the functions \code{SelectPropFit()} and 
#'           \code{SelectPropFitM()}.
#'           Its runtime is \eqn{O(n . ln(n))}.
#'
#' @section Credits:
#'          The code of this function has been written by 
#'          Fabian Aisenbrey.
#'
#' @references  Holland, John (1975): 
#'              \emph{Adaptation in Natural and Artificial Systems},  
#'              The University of Michigan Press, Ann Arbor.
#'              (ISBN:0-472-08460-7)
#'
#' @param fit     Fitness vector.
#' @param lF      Local configuration.
#' @param size    Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFitOnln(fit, NewlFselectGenes()) 
#' SelectPropFitOnln(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFitOnln<-
function(fit, lF, size=1)
{ minimum<-min(fit)
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
  weights<-fit/sum(fit)
  randoms<-sort(stats::runif(size))
  index<-rep(0, size)
  outputIndex<-1; weightIndex<-1; cumWeight<-0.0;
  while (outputIndex<=size && weightIndex <= length(weights)) {
	  cumWeight <- cumWeight + weights[weightIndex]
  while (outputIndex<= size && randoms[outputIndex]<cumWeight) {
        index[outputIndex]<-weightIndex
        outputIndex=outputIndex+1}
  weightIndex<-weightIndex+1
}
#  return(index[sample((1:size), size)])
  return(index)
}

#' Selection proportional to fitness \eqn{O(n^2)}.
#'
#' @description \code{SelectPropFit()} implements selection
#'              proportional to fitness. Negative fitness
#'              vectors are shifted to \eqn{R^+}.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              Holland's schema theorem uses this selection function.
#'              See John Holland (1975) for further information.
#'
#' @section Warning:
#'              There is a potential slow for-loop in the code.
#' 
#' @references  Holland, John (1975): 
#'              \emph{Adaptation in Natural and Artificial Systems},  
#'              The University of Michigan Press, Ann Arbor.
#'              (ISBN:0-472-08460-7)
#'
#' @param fit   Fitness vector.
#' @param lF    Local configuration.
#' @param size  Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFit(fit, NewlFselectGenes()) 
#' SelectPropFit(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFit<-
function(fit, lF, size=1)
{ minimum<-min(fit)
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
  slots<-cumsum(fit)/sum(fit)
  index<-rep(0, size)
  for (i in (1:size))
  { index[i]<-1+sum(1*slots<stats::runif(1)) }
  return(index)
}

#' Selection proportional to fitness (vector/matrix).
#'
#' @description \code{SelectPropFitM()} implements selection
#'              proportional to fitness. Negative fitness
#'              vectors are shifted to \eqn{R^+}.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              Holland's schema theorem uses this selection function.
#'              See John Holland (1975) for further information.
#'
#' @section Warning:
#'              The code is completely written in vector/matrix 
#'              operations. 
#'              \code{outer} uses \eqn{O(n^2)} memory cells. 
#' 
#' @references  Holland, John (1975): 
#'              \emph{Adaptation in Natural and Artificial Systems},  
#'              The University of Michigan Press, Ann Arbor.
#'              (ISBN:0-472-08460-7)
#'
#' @param fit   Fitness vector.
#' @param lF    Local configuration.
#' @param size  Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFitM(fit, NewlFselectGenes()) 
#' SelectPropFitM(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFitM<-
function(fit, lF, size=1)
{ if (min(fit)<= 0) {fit<-lF$Offset()+abs(min(fit))+fit}
  1+colSums(1*outer((cumsum(fit)/sum(fit)), stats::runif(rep(1,size)), FUN="<"))
}

#' Selection proportional to fitness differences O(n ln(n)).
#'
#' @description \code{SelectPropFitDiffOnln()} implements selection
#'              proportional to fitness differences. Negative fitness
#'              vectors are shifted to \eqn{R^+}.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              Holland's schema theorem uses this selection function.
#'              See John Holland (1975) for further information.
#'
#' @details This is a fast implementation which gives exactly the same 
#'           results as the functions \code{SelectPropFitDiff()} 
#'           and \code{SelectPropDiffFitM()}.
#'           Its runtime is \eqn{O(n . ln(n))}.
#'
#'          An epsilon (\code{lF$Eps()}) is added to the fitness 
#'          difference vector. This guarantees numerical stability,
#'          even if all genes in the population have the same fitness. 
#'
#' @section Credits:
#'          The code of this function has been adapted by 
#'          Fabian Aisenbrey.
#'
#' @section Warning:
#'              There is a potential slow for-loop in the code.
#' 
#' @references  Holland, John (1975): 
#'              \emph{Adaptation in Natural and Artificial Systems},  
#'              The University of Michigan Press, Ann Arbor.
#'              (ISBN:0-472-08460-7)
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFitDiffOnln(fit, NewlFselectGenes()) 
#' SelectPropFitOnln(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFitDiffOnln<-
function(fit, lF, size=1)
{ minimum<-min(fit)
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
  fitdiff<-lF$Eps()+fit-minimum
  weights<-fitdiff/sum(fitdiff)
  randoms<-sort(stats::runif(size))
  index<-rep(0, size)
  outputIndex<-1; weightIndex<-1; cumWeight<-0.0;
  while (outputIndex<=size && weightIndex <= length(weights)) {
	  cumWeight <- cumWeight + weights[weightIndex]
  while (outputIndex<= size && randoms[outputIndex]<cumWeight) {
        index[outputIndex]<-weightIndex
        outputIndex=outputIndex+1}
  weightIndex<-weightIndex+1
}
#  return(index[sample((1:size), size)])
  return(index)
}

#' Selection proportional to fitness differences. 
#'
#' @description \code{SelectPropFitDiff()} implements selection
#'              proportional to fitness differences.
#'              It selects a gene out of the population
#'              with a probability proportional to the fitness
#'              difference to the gene with minimal fitness.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              The fitness of survival of the gene with 
#'              minimal fitness is set by \code{lF$Eps()} 
#'              to \code{0.01} per default.
#'              See equation (7.45) Andreas Geyer-Schulz (1997), p. 205.
#'
#' @references   Geyer-Schulz, Andreas (1997):
#'      \emph{Fuzzy Rule-Based Expert Systems and Genetic Machine Learning},
#'      Physica, Heidelberg.
#'      (ISBN:978-3-7908-0830-X)
#'
#' @section Note: 
#'        \code{SelectPropFitDiff()} is a dynamic scaling function.
#'        Complexity: \eqn{O(n^2)}.
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFitDiff(fit, NewlFselectGenes()) 
#' SelectPropFitDiff(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFitDiff<-
function(fit, lF, size=1)
{ minimum<-min(fit)
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
  fitdiff<-lF$Eps()+fit-minimum
  slots<-cumsum(fitdiff)/sum(fitdiff)
  index<-rep(0, size)
  for (i in (1:size))
  { index[i]<-1+sum(1*slots<stats::runif(1)) }
  return(index)
}

#' Selection proportional to fitness differences. 
#'
#' @description \code{SelectPropFitDiffM()} implements selection
#'              proportional to fitness differences.
#'              It selects a gene from the population
#'              with a probability proportional to the fitness
#'              difference to the gene with minimal fitness.
#'              The default of the function \code{lF$Offset()} is \code{1}. 
#'              The fitness of survival of the gene with 
#'              minimal fitness is set by \code{lF$Eps()} 
#'              to \code{0.01} per default.
#'              See equation (7.45) Andreas Geyer-Schulz (1997), p. 205.
#'
#' @section Warning:
#'              \code{outer} uses \eqn{O(n^2)} memory cells. 
#'
#' @references   Andreas Geyer-Schulz (1997):
#'          \emph{Fuzzy Rule-Based Expert Systems and Genetic Machine Learning},
#'                Physica, Heidelberg.
#'           <978-3-7908-0830-X>
#'
#' @section Note: 
#'        \code{SelectPopFitDiff()} is a dynamic scaling function.
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectPropFitDiffM(fit, NewlFselectGenes()) 
#' SelectPropFitDiffM(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectPropFitDiffM<- 
function(fit, lF, size=1)
{1+colSums(1*
outer((cumsum(lF$Eps()+(fit-min(fit)))/(lF$Eps()+sum(fit-min(fit)))),
       stats::runif(rep(1,size)),FUN="<"))}

#' Selection with uniform probability.
#'
#' @description \code{SelectUniform()} implements selection
#'              by choosing a gene with equal probability. 
#'
#' @details  This selection function is useful:
#'              \enumerate{
#'              \item
#'              To specify mating behavior in crossover operators.
#'              \item
#'              For computer experiments without selection pressure.
#'              \item 
#'              For computing random search solutions as a benchmark.
#'              }
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectUniform(fit, NewlFselectGenes()) 
#' SelectUniform(fit, NewlFselectGenes(), length(fit)) 
#' @export
SelectUniform<-
function(fit, lF, size=1)
{ sample(length(fit), size=size, replace=TRUE)}

#' Selection with uniform probability without replacement.
#'
#' @description \code{SelectUniformP()} implements selection
#'              by choosing a gene with equal probability without 
#'              replacement.
#'              Usage:
#'              \enumerate{
#'              \item
#'              To specify mating behavior in crossover operators.
#'              \item
#'              For computer experiments without selection pressure.
#'              }
#'
#' @details Selection without replacement guarantees that 
#'          vectors of different indices are selected.
#'          A vector of the size of the population is a permutation 
#'          of indices. This property is needed for the classic 
#'          variant of differential evolution.
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @references 
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134. 
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectUniformP(fit, NewlFselectGenes()) 
#' SelectUniformP(fit, NewlFselectGenes(), length(fit)) 
#' @export
SelectUniformP<-
function(fit, lF, size=1)
{ if (size <=length(fit)) 
	{replace=FALSE} 
else 
	{replace=TRUE} 
	sample(length(fit), size=size, replace=replace)}

#' Deterministic duel. 
#' 
#' @description \code{SelectDuel()} implements selection
#'              by a tournament between 2
#'              randomly selected genes. The best gene always wins.
#'              This is the version of tournament selection 
#'              with the least selection pressure.
#'
#' @details This is an O(n) implementation of tournament selection with 
#'          a tournament size of 2. 
#'
#' @details A special case of tournament selection.
#'              
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectDuel(fit, NewlFselectGenes()) 
#' SelectDuel(fit, NewlFselectGenes(), length(fit)) 
#' @export
SelectDuel<-function(fit, lF, size=1)
{ l<-length(fit)
cand1<-sample(1:l, size, replace=TRUE)
cand2<-sample(1:l, size, replace=TRUE)
cand1wins<-fit[cand1]>fit[cand2]
i1<-cand1[cand1wins]
i2<-cand2[!cand1wins]
return(c(i1, i2))
}

#' Deterministic tournament of size \code{k}.
#' 
#' @description \code{Tournament()} is implemented in two steps:
#' \enumerate{
#' \item
#' A subset of size k of the population is selected with uniform probability.
#' \item
#' A gene is selected with probability proportional to fitness.
#' }
#'
#' @details In each generation, the worst \code{k-1} genes 
#'          in a population do not survive.
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#'
#' @return Index of the best candidate.
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' Tournament(fit, NewlFselectGenes()) 
#' @export
Tournament<-function(fit, lF)
{ cand<-sample(length(fit), size=lF$TournamentSize())
 (cand[max(fit[cand])==fit[cand]])[1] }

#' Stochastic tournament of size \code{k}.
#' 
#' @description \code{STournament()} is implemented in two steps:
#' \enumerate{
#' \item
#' A subset of size k of the population is selected with uniform probability.
#' \item
#' A gene is selected with probability proportional to fitness.
#' }
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#'
#' @return Index of candidate.
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' STournament(fit, NewlFselectGenes()) 
#' @export
STournament<-function(fit, lF)
{ cand<-sample(length(fit), size=lF$TournamentSize())
  cand[SelectPropFit(fit[cand], lF)]}

#' Tournament selection. 
#'
#' @description \code{SelectTournament()} implements selection
#'              by doing a tournament between \code{lF$TournamentSize()} 
#'              randomly selected genes. The best gene always wins.
#'              The default of \code{lF$TournamentSize()} is \code{2}. This 
#'              is the version with the least selection pressure.
#'              
#'              \code{lF$TournamentSize()} must be less 
#'                    than the population size.
#'
#' @param fit Fitness vector.
#' @param lF  Local configuration.
#' @param size Number of selected genes. Default: \code{1}. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectTournament(fit, NewlFselectGenes()) 
#' SelectTournament(fit, NewlFselectGenes(), length(fit)) 
#' @export
SelectTournament<-function(fit, lF, size=1)
{ 
index<-rep(0, size)
for (i in (1:size)) 
{ index[i]<-Tournament(fit, lF) }
return(index)
}

#' Stochastic tournament selection. 
#'
#' @description \code{SelectSTournament()} implements selection
#'              through a stochastic tournament between 
#'               \code{lF$TournamentSize()} 
#'              randomly selected genes. A gene wins a tournament 
#'              with a probability proportional to its fitness.
#'              The default of \code{lF$TournamentSize()} is \code{2}. 
#'              A tournament
#'              with 2 participants has the least selection pressure.
#'              
#'              \code{lF$TournamentSize()} must be less 
#'                    than the population size.
#'
#' @param fit     Fitness vector.
#' @param lF      Local configuration.
#' @param size    Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectSTournament(fit, NewlFselectGenes()) 
#' SelectSTournament(fit, NewlFselectGenes(), length(fit)) 
#' @export
SelectSTournament<-function(fit, lF, size=1)
{ 
index<-rep(0, size)
for (i in (1:size)) 
{ index[i]<-STournament(fit, lF) }
return(index)
}

#' Stochastic universal sampling.
#'
#' @description \code{SelectSUS()} implements selection
#'      by Baker's stochastic universal sampling method. 
#'      SUS is a strictly sequential algorithm 
#'      which has zero bias and minimal spread.
#'      SUS uses a single random number for each generation.
#'      See Baker, James E. (1987), p. 16.
#'               
#' @references Baker, James E. (1987):
#'      Reducing Bias and Inefficiency in the Selection Algorithm.
#'      In Grefenstette, John J.(Ed.) 
#'      \emph{Proceedings of the Second International 
#'      Conference on Genetic Algorithms on Genetic Algorithms}, pp. 14-21.
#'      (ISBN:978-08058-0158-8)
#'
#' @param fit   Fitness vector.
#' @param lF    Local configuration.
#' @param size  Number of selected genes. Default: 1. 
#'
#' @return The index vector of the selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectSUS(fit, NewlFselectGenes()) 
#' SelectSUS(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectSUS<-function(fit, lF, size=1)
{ minimum<-min(fit)	
  if (minimum<= 0) {fit<-lF$Offset()+abs(minimum)+fit}
sum<-0; pos<-1
ptr<-stats::runif(1)
expVal<-fit/mean(fit)
index<-rep(0, size)
for (j in (1:length(fit)))
{ sum<-sum+expVal[j] 
  while (sum>ptr) {index[pos]<-j; pos<-pos+1; ptr<-ptr+1} }
return(index[sample(length(index), size=size, replace=(size>length(fit)))]) }

#' Linear rank selection with selective pressure.
#'
#' @description \code{SelectLRSelective()} implements selection
#'      by Whitley's linear rank selection with selective pressure 
#'      for the GENITOR algorithm. 
#'      See Whitley, Darrell (1989), p. 121.
#'               
#' @details The selection pressure is configured by the constant function 
#'          \code{lF$SelectionBias()}. Its values should be strictly larger 
#'          than 1 and preferably below 2. The default is set to 1.5.
#'          A value of 1.0 means uniform random selection. 
#'               
#' @references Whitley, Darrell (1989):
#' The GENITOR Algorithm and Selection Pressure.
#' Why Rank-Based Allocation of Reproductive Trials is Best.
#'      In Schaffer, J. David (Ed.) 
#'      \emph{Proceedings of the Third International 
#'      Conference on Genetic Algorithms on Genetic Algorithms}, pp. 116-121.
#'      (ISBN:1-55860-066-3)
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Size of return vector (default: 1).
#'
#' @return The index vector of selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectLRSelective(fit, NewlFselectGenes()) 
#' SelectLRSelective(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectLRSelective<- function(fit, lF, size=1) {
sB<-lF$SelectionBias()
i<-1.0+length(fit)*(sB-sqrt(sB*sB-4.0*(sB-1.0)*stats::runif(rep(1, size))))/2.0/(sB-1.0)
f<-sort(fit, decreasing=TRUE, index.return=TRUE)
return(f$ix[as.integer(i)])
}

#' Linear rank selection with interpolated target sampling rates.
#'
#' @description \code{SelectLinearRankTSR()} implements selection
#'      with interpolated target sampling rates.
#'               
#' @details The target sampling rate is a linear interpolation 
#'          between \code{lF$MaxTSR()} and \code{Min<-2-lF$MaxTSR()},
#'          because the sum of the target sampling rates is $n$. 
#'          The target sampling rates are computed and used as a fitness
#'          vector for stochastic universal sampling algorithm
#'          implemented by \code{SelectSUS()}.
#'          \code{lF$MaxTSR()} should be in [1.0, 2.0].
#'               
#'          TODO: More efficient implementation. We use two sorts!
#'
#' @references Grefenstette, John J. and Baker, James E. (1989):
#'      How Genetic Algorithms Work: A Critical Look at Implicit Parallelism
#'      In Schaffer, J. David (Ed.) 
#'      \emph{Proceedings of the Third International 
#'      Conference on Genetic Algorithms on Genetic Algorithms}, pp. 20-27.
#'      (ISBN:1-55860-066-3)
#'
#' @param fit    Fitness vector.
#' @param lF     Local configuration.
#' @param size   Size of return vector (default: 1).
#'
#' @return The index vector of selected genes.
#'
#' @family Selection Functions
#'
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' SelectLinearRankTSR(fit, NewlFselectGenes()) 
#' SelectLinearRankTSR(fit, NewlFselectGenes(), length(fit)) 
#' @importFrom stats runif
#' @export
SelectLinearRankTSR <- function(fit, lF, size=1) {
Max<-lF$MaxTSR()
Min<-2-Max
n<-length(fit)
tsr<-Min+(Max-Min)*((n:1)-1)/n 
f<-sort(fit, decreasing=TRUE, index.return=TRUE)
fix<-sort(f$ix, decreasing=FALSE, index.return=TRUE)
ftsr<-tsr[fix$ix]
return(SelectSUS(ftsr, lF, size))
}

#
# Part 3. Configuration of selection functions.
#

#' Convert a selection function into a continuation.
#' 
#' @description \code{TransformSelect()} precomputes 
#'              the indices of genes to be selected and
#'              converts the selection function into an access function to the 
#'              next index. 
#'              The access function provides a periodic random 
#'              index stream with a period length of the population size.
#'              In a genetic algorithm with a fixed size population,
#'              this avoids recomputation of the selection functions 
#'              for each gene 
#'              and its mate. 
#'               
#' @details  The motivation for this transformation is: 
#'           \enumerate{
#'           \item We avoid the recomputation of potentially 
#'                 expensive selection functions.
#'                 E.g. In population-based genetic algorithms, 
#'                 the selection function
#'                 is computed twice per generation instead of 
#'                 more than generation times the population size. 
#'           \item No additional control flow is needed.
#'           \item Dynamic reconfiguration is possible.
#'           \item All selection functions have a 
#'                 common abstract interface and, 
#'                 therefore, can be overloaded by 
#'                 specialized concrete implementations.
#'                 (Polymorphism).
#'           } 
#' 
#'           The implementation idea is adapted from the
#'           continuation passing style 
#'           in functional programming. See Reynolds, J. C. (1993).
#' 
#' @section Parallelization/Distribution:
#'           \enumerate{         
#'           \item We use this tranformation if only the evaluation of genes 
#'                 should be parallelized/distributed.
#'           \item If the complete replication of genes is parallelized, 
#'                 this transformation cannot be used in its current form.
#'                 The current implementations of the selection functions
#'                 can not easily be parallelized.
#'                 }
#'
#' @references Reynolds, J. C. (1993):
#'             The discoveries of continuations.
#'             \emph{LISP and Symbolic Computation} 6, 233-247. 
#'             <doi:10.1007/BF01019459>
#'
#' @param fit          Fitness vector.
#' @param lF           Local configuration.
#' @param SelectFUN    Selection function.
#'
#' @return A function with a  state which consists of 
#'         the precomputed gene index 
#'         vector, its \code{length}, and a \code{counter}. 
#'         The function increments the counter in the state 
#'         of its environment and
#'         returns the precomputed gene index at position 
#'         \code{modulo((counter+1),length)} in   
#'         the precomputed index vector in its environment.
#'         The function supports the same interface as a selection function.
#'         
#' @family Performance Optimization
#' 
#' @examples
#' fit<-sample(10, 15, replace=TRUE)
#' newselect<-TransformSelect(fit, NewlFselectGenes(), SelectSUS) 
#' newselect(fit, NewlFselectGenes())
#' newselect(fit, NewlFselectGenes(), 5)
#' newselect(fit, NewlFselectGenes(), 10)
#' newselect(fit, NewlFselectGenes(), 10)
#' @export
TransformSelect<-function(fit, lF, SelectFUN)
{ c<-0
  index<-SelectFUN(fit, lF, length(fit))
  ilen<-length(index)
  return(function(fit, lF, size=1) 
	 { i<-index[1+(c+seq(1:size))%%ilen]; c<<-c+size; i}) }

#' Configure the selection function of a genetic algorithm.
#'
#' @description \code{SelectGeneFactory()} implements selection
#'              of one of the gene selection functions in this 
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.             
#'
#'              Current support:
#' 
#'              \enumerate{
#'              \item "Uniform" returns \code{SelectUniform()}.
#'              \item "UniformP" returns \code{SelectUniformP()}.
#'              \item "ProportionalOnln" returns \code{SelectPropFitOnln()}.
#'              \item "Proportional" returns \code{SelectPropFit()}.
#'              \item "ProportionalM" returns \code{SelectPropFitM()}.
#'              \item "PropFitDiffOnln" returns \code{SelectPropFitDiffOnln()}.
#'              \item "PropFitDiff" returns \code{SelectPropFitDiff()}.
#'              \item "PropFitDiffM" returns \code{SelectPropFitDiffM()}.
#'              \item "Tournament" returns \code{SelectTournament()}.
#'              \item "STournament" returns \code{SelectSTournament()}.
#'              \item "Duel" returns \code{SelectDuel()}.
#'              \item "LRSelective" returns \code{SelectLRSelective()}.
#'              \item "LRTSR" returns \code{SelectLinearRankTSR()}.
#'              \item "SUS" returns \code{SelectSUS()}.
#'              }
#'
#' @details
#'             If \code{SelectionContinuation()==TRUE} then:
#'              \enumerate{
#'                \item In package xegaPopulation in
#'                      function \code{NextPopulation()},
#'                      first the functions 
#'                        \code{SelectGene()} and \code{SelectMate()} 
#'                        are transformed by \code{TransformSelect()} to
#'                        a continuation function with
#'                      embedded index vector and counter.
#'                \item For each call in \code{ReplicateGene()},
#'                      \code{SelectGene} and \code{SelectMate()} 
#'                      return the index of the selected gene.
#'                      }
#'
#' @param method A string specifying the selection function. 
#'
#' @return A selection function for genes.
#'
#' @family Configuration 
#'
#' @examples
#' SelectGene<-SelectGeneFactory("Uniform")
#' fit<-sample(10, 15, replace=TRUE)
#' SelectGene(fit, lFselectGenes) 
#' sel<-"Proportional"
#' SelectGene<-SelectGeneFactory(method=sel)
#' fit<-sample(10, 15, replace=TRUE)
#' SelectGene(fit, lFselectGenes) 
#' @export
SelectGeneFactory<-function(method="PropFitDiffOnln") {
if (method=="Uniform") {f<- SelectUniform}
if (method=="UniformP") {f<- SelectUniformP}
if (method=="ProportionalOnln") {f<-SelectPropFitOnln}
if (method=="Proportional") {f<-SelectPropFit}
if (method=="ProportionalM") {f<-SelectPropFitM}
if (method=="PropFitDiffOnln") {f<-SelectPropFitDiffOnln}
if (method=="PropFitDiff") {f<-SelectPropFitDiff}
if (method=="PropFitDiffM") {f<-SelectPropFitDiffM}
if (method=="Duel") {f<- SelectDuel}
if (method=="Tournament") {f<- SelectTournament}
if (method=="STournament") {f<- SelectSTournament}
if (method=="LRSelective") {f<- SelectLRSelective}
if (method=="LRTSR") {f<- SelectLinearRankTSR}
if (method=="SUS") {f<- SelectSUS}
if (!exists("f", inherits=FALSE)) 
	{stop("Selection label ", method, " does not exist")}
return(f)
}

# end of file
