#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions
#            Representation independent.
#     Package: xegaEvalGene
#

#' Generate local functions and objects
#' 
#'@description \code{NewlFevalGenes()} returns 
#'              the list of functions containing
#'              a definition of all local objects required for the use
#'              of evaluation functions. We reference this object 
#'              as local configuration. When adding additional 
#'              evaluation functions, this must be extended
#'              by the constant (functions) needed to configure them.
#'
#'@param  penv A problem environment.
#'@return The local configuration. A list of functions.
#'
#'@examples
#' Parabola2D<-Parabola2DFactory()
#' lF<-NewlFevalGenes(Parabola2D)
#' lF$Max()
#'@export
NewlFevalGenes<-function(penv)
{
         list(
         Max=function() {1},
         ReportEvalErrors=function() {TRUE},
	 DecodeGene=function(gene, lF) {gene$gene1},
         penv=penv
         )
}

#' Evaluates a gene in a problem environment
#'
#' @description \code{EvalGeneU()} evaluates a gene in
#'              a problem environment. 
#'
#' @details If the evaluation of the fitness function of the 
#'          problem environment fails, the following 
#'          strategy is used: 
#'          We catch the error
#'          and print it, ignore it: 
#'
#'          The error handler returns \code{NA}.
#'          
#'          We check for the error and update the gene:
#'          \itemize{
#'          \item \code{$evaluated}  TRUE.
#'          \item \code{$evalFail}   TRUE.
#'          \item \code{$fit} is set to the minimum fitness in the 
#'                 population.
#'          }
#'
#'          The boolean function \code{lF$ReportEvalErrors()} controls
#'          the output of error messages for evaluation failures. 
#'          Rationale: In grammatical evolution, the standard approach 
#'          ignores attempts the evaluate incomplete programs.
#'
#' @section Future improvement:
#'             Provide configurable error handlers.
#'             Rationale: Make debugging for new problem environments easier.
#'             Catch communication problems in distributed/parallel 
#'             environments.
#'          
#' @param gene  A gene.
#' @param lF    The local configuration of the genetic algorithm.
#'
#' @return A gene (with \code{$evaluated==TRUE}).
#'
#' @family Evaluation Functions
#'            
#' @examples
#' Parabola2D<-Parabola2DFactory()
#' lF<-NewlFevalGenes(Parabola2D)
#' g1<-list(evaluated=FALSE, fit=0, gene1=c(1.0, -1.5, 3.37))
#' g2<-list(evaluated=FALSE, fit=0, gene1=c(0.0, 0.0, 0.0))
#' EvalGeneU(g1, lF)
#' EvalGeneU(g2, lF)
#' @export
EvalGeneU<-function(gene, lF)
{
ng<-gene
ReportEvalErrors<-lF$ReportEvalErrors
ng$fit<-tryCatch(
	lF$Max()*lF$penv$f(lF$DecodeGene(gene, lF), gene, lF),
	error = function(e) 
	         {if (ReportEvalErrors())
                  {message("EvalGeneU:")
                   message(conditionMessage(e))}
	          NA})
if (is.na(ng$fit)) # ignore and report!
{  ng$fit<-lF$CWorstFitness()  
   ng$evalFail<-TRUE} else 
{ng$evalFail<-FALSE}
ng$evaluated<-TRUE
return(ng)
}

#' Evaluates a repaired gene in a problem environment.
#'
#' @description \code{EvalGeneR()} evaluates a repaired gene in
#'              a problem environment.
#'
#' @details If the decoder repairs a gene, the repaired gene 
#'          must replace the original gene.
#'
#' @param gene  A gene.
#' @param lF    The local configuration of the genetic algorithm.
#'
#' @return A gene (with \code{$evaluated==TRUE}).
#'
#' @family Evaluation Functions
#' @seealso   EvalGeneU
#'            
#' @examples
#' Parabola2D<-Parabola2DFactory()
#' lF<-NewlFevalGenes(Parabola2D)
#' g1<-list(evaluated=FALSE, fit=0, gene1=c(1.0, -1.5, 3.37))
#' g2<-list(evaluated=FALSE, fit=0, gene1=c(0.0, 0.0, 0.0))
#' EvalGeneR(g1, lF)
#' EvalGeneR(g2, lF)
#' @export
EvalGeneR<-function(gene, lF)
{
ng<-gene
dg<-lF$DecodeGene(gene, lF)
if (!is.null(names(dg)))
{   p<-dg$parm; ng<-dg$gene }
else { p<-dg }
ReportEvalErrors<-lF$ReportEvalErrors
ng$fit<-tryCatch(
	lF$Max()*lF$penv$f(p, ng, lF),
	error = function(e) 
	         {if (ReportEvalErrors())
                  {message("EvalGeneR:")
                   message(conditionMessage(e))}
	          NA})
if (is.na(ng$fit)) # ignore and report!
{  ng$fit<-lF$CWorstFitness()  
   ng$evalFail<-TRUE} else 
{ng$evalFail<-FALSE}
ng$evaluated<-TRUE
return(ng)
}

#' Evaluates a gene in a deterministic problem environment. 
#' 
#' @description \code{EvalGeneDet()} evaluates a gene in
#'              a problem environment if it has not been evaluated yet.
#'              The repeated evaluations of a gene are omitted.
#'
#' @details If the evaluation of the fitness function of the 
#'          problem environment fails, we catch the error and 
#'          return \code{NA}.
#'
#' @param gene A gene.
#' @param lF   The local configuration of the genetic algorithm.
#'
#' @return A gene (with \code{$evaluated==TRUE}).
#'
#' @family Evaluation Functions
#'            
#' @examples
#' Parabola2D<-Parabola2DFactory()
#' lF<-NewlFevalGenes(Parabola2D)
#' g1<-list(evaluated=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g2<-list(evaluated=FALSE, fit=0, gene1=c(0.0, 0.0))
#' g1a<-EvalGeneDet(g1, lF)
#' EvalGeneDet(g1a, lF)
#' g2a<-EvalGeneDet(g2, lF)
#' EvalGeneDet(g2a, lF)
#' @export
EvalGeneDet<-function(gene, lF)
{
if (gene$evaluated) gene else EvalGeneU(gene, lF)
}

#' Evaluates a gene in a stochastic problem environment. 
#' 
#' @description \code{EvalGeneStoch()} evaluates a gene in
#'              a stochastic problem environment.
#'
#' @details In a stochastic problem environment, the expected fitness
#'          is maximized. The computation of the expectation is 
#'          done by incrementally updating the mean.
#'          For this, need the number of evaluations of the gene 
#'          (\code{$obs} of the gene).
#'          In addition, we compute the incremental variance 
#'          of the expected fitness 
#'          stored in \code{$var}.  
#'          The standard deviation is then \code{gene$var/gene$obs}.
#'
#'          If the evaluation of the fitness function of the 
#'          problem environment fails, we catch the error and 
#'          return \code{NA} for the first evaluation of the gene.
#'          If the gene has been evaluated, we return the old gene.
#'
#' @param gene A gene.
#' @param lF   The local configuration of the genetic algorithm.
#'
#' @return A gene with the elements
#'   \itemize{
#'   \item  \code{$evaluated}: Boolean.
#'   \item  \code{$evalFail}:  Boolean. 
#'   \item  \code{$fit}:       Mean fitness of gene. 
#'   \item  \code{$gene1}:     Gene.  
#'   \item  \code{$obs}:       Number of evaluations of gene.
#'   \item  \code{$var}:       Variance of fitness.  
#'   \item  \code{$sigma}:     Standard deviation of fitness.
#'   }
#'  
#' @family Evaluation Functions 
#'
#' @examples
#' DeJongF4<-DeJongF4Factory()
#' lF<-NewlFevalGenes(DeJongF4)
#' g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g1
#' g2<-EvalGeneStoch(g1, lF)
#' g2
#' g3<-EvalGeneStoch(g2, lF)
#' g3
#' g4<-EvalGeneStoch(g3, lF)
#' g4
#' g5<-EvalGeneStoch(g4, lF)
#' g5
#' @export
EvalGeneStoch<-function(gene, lF)
{
if (!gene$evaluated) 
{
 ng<-EvalGeneU(gene, lF)	
if (ng$evalFail==TRUE) {ng$obs<-0; ng$var<-0; return(ng)}
 ng$obs<-1
 ng$var<-0
 ng$sigma<-0
return(ng)}
else
{
ng<-EvalGeneU(gene, lF)
if (ng$evalFail==TRUE) {return(ng)}
ng$obs<-gene$obs+1
fit<-ng$fit
ng$fit<- (gene$obs)/(ng$obs)*gene$fit + (fit/(ng$obs))
ng$var<-gene$var+(fit-gene$fit)*(fit-ng$fit)
ng$sigma<-sqrt(ng$var/ng$obs)
return(ng)}
}

#' Test of incremental mean, variance, and standard deviation.
#'
#' @param gene A gene.
#' @param lF   A local function list with a problem environment.
#' @param rep  Number of repeated evaluations.
#'
#' @return A gene.
#'
#' @family Tests
#'
#' @examples
#' DeJongF4<-DeJongF4Factory()
#' lF<-NewlFevalGenes(DeJongF4)
#' g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g10<-testEvalGeneStoch(g1, lF, rep=10)
#' g10
#' @export
testEvalGeneStoch<-function(gene, lF, rep)
{ ng<-gene
	for (i in (1:rep)) {ng1<-EvalGeneStoch(ng, lF); ng<-ng1}
  return(ng)}

#' Configure the evaluation function of a genetic algorithm.
#' 
#' @description \code{EvalGeneFactory()} implements the selection 
#'              of one of the evaluation functions for a gene 
#'              in this package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#' @param method   Available methods are: 
#'
#' \itemize{
#' \item "EvalGeneU": Evaluate gene (Default).
#'                    Function \code{EvalGeneU}.
#' \item "EvalGeneR": If the gene has been repaired by a decoder, 
#'                    the gene is replaced by the repaired gene.
#'                    Function \code{EvalGeneR}.
#' \item "Deterministic": A gene which has been evaluated is 
#'                        not reevaluated. 
#'                    Function \code{EvalGeneDet}.
#' \item "Stochastic": The fitness mean and 
#'                     the fitness variance
#'                     are incrementally updated.
#'                     Genes remaining in the population over several
#'                     generations, the fitness mean converges to the 
#'                     expected mean.
#'                    Function \code{EvalGeneStoch}.
#' }            
#'
#' @return An evaluation function.
#'
#' @family Configuration
#'
#' @examples
#' set.seed(5)
#' DeJongF4<-DeJongF4Factory()
#' lF<-NewlFevalGenes(DeJongF4)
#' EvalGene<-EvalGeneFactory("EvalGeneU")
#' g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g1
#' g2<-EvalGene(g1, lF)
#' g2
#' EvalGene<-EvalGeneFactory("Deterministic")
#' g3<-EvalGene(g2, lF)
#' g3
#' set.seed(5)
#' EvalGene<-EvalGeneFactory("Stochastic")
#' g1<-list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g1
#' g2<-EvalGene(g1, lF)
#' g2
#' g3<-EvalGene(g2, lF)
#' g3
#' @export
EvalGeneFactory<-function(method="EvalGeneU") {
if (method=="EvalGeneU") {f<- EvalGeneU}
if (method=="EvalGeneR") {f<- EvalGeneR}
if (method=="Deterministic") {f<- EvalGeneDet}
if (method=="Stochastic") {f<- EvalGeneStoch}
if (!exists("f", inherits=FALSE))
        {stop("Evaluation method label ", method, " does not exist")}
return(f)
}

#' Evaluate a gene
#'
#' @description \code{EvalGene()} is the abstract function which evaluates 
#'              a gene. 
#'      
#' @details For minimization problems, the fitness value 
#'          must be multiplied by -1. The constant function 
#'          \code{lF$Max()} returns 1 for a maximization and 
#'          -1 for a minimization problem.
#' 
#' @param gene       A gene (representation independent).
#' @param lF         The local configuration (a function factory provided by 
#'                          the \code{xegaX} package).
#'
#' @return A gene.
#'
#' @family Evaluation Functions
#'
#' @examples
#' DeJongF4<-DeJongF4Factory()
#' lF<-NewlFevalGenes(DeJongF4)
#' EvalGene<-EvalGeneFactory()
#' g1<-list(evaluated=FALSE, fit=0, gene1=c(1.0, -1.5))
#' g1
#' g2<-EvalGene(g1, lF)
#' g2
#' @export
EvalGene<-EvalGeneFactory(method="EvalGeneU")

