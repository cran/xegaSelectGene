
#' Selection functions for genetic algorithms.
#' 
#' The \code{selectGene} package provides selection and scaling functions
#' for genetic algorithms. 
#' All functions of this package are independent of the gene 
#' representation. 
#'
#' \itemize{
#'    \item Scaling functions and dispersion measures are in scaling.R
#'    \item Selection functions are in selectGene.R.
#'          For selection functions, a transformation to index access 
#'          functions is provided (a limited form of function continuation).
#'    \item Benchmark functions for selection functions are in 
#'          selectGeneBenchmark.R.
#'          Except for uniform selection, the continuation form of
#'          selection functions should be used.
#'    \item Evaluation functions are in evalGene.R.
#'    \item Counting and timing of function executions  
#'          are provided by transformation functions in timer.R
#'    \item Problem environments for examples and unit tests for
#'      \itemize{
#'         \item function optimization:
#'          DeJongF4.R (stochastic functions) and Parabola2D.R 
#'          (delayed execution for benchmarking of parallelism,
#'           deterministic function,
#'           deterministic function with early termination check, 
#'           function with random failures)
#'         \item combinatorial optimization: 
#'             newTSP.R (for the traveling salesman problem).
#'         \item boolean function learning:
#'             newXOR.R (for the XOR problem).
#'           }
#'         }
#'
#' @section Interface Scaling Functions:
#'
#'              All scaling functions must implement
#'              the following abstract interface:
#'
#'              \code{function name}(\code{fit}, \code{lF})
#'
#'              \strong{Parameters}
#'
#'              \itemize{
#'              \item \code{fit}  A fitness vector. 
#'              \item \code{lF}   Local configuration.
#'              }
#'
#'              \strong{Return Value}
#'
#'              Scaled fitness vector.
#'
#' @section Interface Dispersion Measures:
#'
#'              All dispersion measure functions must implement
#'              the following abstract interface:
#'
#'              \code{function name}(\code{popstatvec})
#'
#'              \strong{Parameters}
#'
#'              \itemize{
#'              \item \code{popstatvec}  Vector of population statistics. 
#'
#' The internal state of the genetic algorithm is described 
#' by a matrix of the history of population statistics.
#' Each row consists of 8 population statistics
#' (mean, min, Q1, median, Q3, max, var, mad). 
#' A row is a vector of population statistics.  
#'
#'              }
#'
#'              \strong{Return Value}
#'
#'              Dispersion measure (real).
#'
#' @section Interface Selection Functions:
#'
#'              All selection functions must implement
#'              the following abstract interface:
#'
#'              \code{function name}(\code{fit}, \code{lF}, \code{size})
#'
#'              \strong{Parameters}
#'
#'              \itemize{
#'              \item \code{fit}  a vector of fitness values.
#'              \item \code{lF}   a local function list.
#'              \item \code{size} the number of indices returned.
#'              }
#'
#'              \strong{Return Value}
#'
#'              A vector of indices of length \code{size}.
#'
#'              All selection functions are implemented
#'              WITHOUT a default assignment to \code{lF}.
#'
#'              A missing configuration should raise an error!
#'
#'              The default value of \code{size} is \code{1}.
#'
#'
#'@section Constants:
#'
#'              Some scaling and selection functions use constants which should 
#'              be configured. 
#'              We handle these constants by 
#'              constant functions created by \code{parm(constant)}.
#'              We store all of these functions in the list of 
#'              local functions \code{lF}.
#'              The rationale is to reduce the number of parameters
#'              of selection functions and to provide a uniform 
#'              interface for selection functions.
#'
#'@section Table of Scaling Constants:
#'
#' \tabular{rcl}{ 
#' \strong{Constant}   \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$Offset           \tab 1                \tab ScaleFitness \cr
#' lF$ScalingExp       \tab 1                \tab ScalingFitness,  \cr 
#'                     \tab                  \tab ThresholdScaleFitness  \cr 
#' lF$ScalingExp2      \tab 1                \tab ThresholdScaleFitness  \cr 
#' lF$ScalingThreshold \tab 1                \tab ThresholdScaleFitness  \cr 
#' lF$RDMWeight        \tab 1.0              \tab ContinuousScaleFitness  \cr 
#' lF$DRMin            \tab 0.5              \tab DispersionRatio  \cr 
#' lF$DRMax            \tab 2.0              \tab DispersionRatio  \cr 
#' lF$ScalingDelay     \tab 1                \tab DispersionRatio  \cr 
#' }
#'
#' \tabular{rcl}{ 
#' \strong{State Variable} \tab \strong{Start Value} \tab \strong{Used in} \cr 
#' lF$RDM           \tab 1.0                \tab ThresholdScaleFitness \cr
#'                  \tab                    \tab ContinuousScaleFitness \cr
#'                  \tab                    \tab xega::RunGA \cr
#' }
#'
#'@section Table of Selection Constants:
#' 
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' lF$SelectionContinuation \tab TRUE \tab xegaPopulation::xegaNextPopulation \cr 
#' lF$Offset           \tab 1                \tab SelectPropFitOnLn \cr
#'                     \tab                  \tab SelectPropFit \cr
#'                     \tab                  \tab SelectPropFitM \cr
#'                     \tab                  \tab SelectPropFitDiffOnLn \cr
#'                     \tab                  \tab SelectPropFitDiff \cr
#'                     \tab                  \tab SUS \cr
#'                     \tab                  \tab SelectLinearRankTSR \cr
#' lF$eps              \tab 0.01             \tab SelectPropFitDiffM \cr
#' lF$TournamentSize   \tab 2                \tab Tournament \cr
#'                     \tab                  \tab SelectTournament \cr
#'                     \tab                  \tab STournament \cr
#'                     \tab                  \tab SelectSTournament \cr
#' lF$SelectionBias    \tab 1.5              \tab SelectLRSelective \cr
#' lF$MaxTSR           \tab 1.5              \tab SelectLinearRankTSR 
#' }
#'
#'@section Parallel/Distributed Execution:
#'
#'     All selection functions in this package return
#'     \enumerate{
#'     \item the index of a selected gene.
#'       The configured selection function is executed each time
#'       a gene must be selected in the gene replication process.
#'       This allows a parallelization/distribution of the
#'       complete gene replication process and the fitness evaluation.
#'       However, the price to pay is a recomputation of the selection
#'       algorithms for each gene and each mate (which may be costly).
#'       The execution time of Baker's SUS function explodes
#'       when used in this way.
#'     \item a vector of indices of the selected genes.
#'        We compute
#'        a vector of indices for genes and their mates,
#'        and we replace the selection function with a
#'        quasi-continuation function
#'        with precomputed indices
#'        which when called, returns the next index.
#'        The selection computation is executed once for each generation
#'        without costly recomputation.
#'        The cost of selecting a gene and its mate is the cost of indexing 
#'        an integer in a vector.
#'        This version is faster for almost all selection functions
#'        (Sequential computation).
#'
#'        The parallelization of quasi-continuation function is not yet implemented. 
#'         }
#'
#' @section Constant Functions for Configuration:
#'
#' The following constant functions are expected to be in the local function list lF.
#'           \itemize{
#'           \item \code{Offset()} in \code{SelectPropFit}:
#'                 Since all fitness values must be larger than 0,
#'                 in case of negative fitness values, \code{Offset()} 
#'                 is the value of the minimum fitness value (default: 1).
#'           \item \code{Eps()} in \code{SelectPropFitDiff}:
#'                 \code{Eps()} is a very small value to eliminate 
#'                 differences of 0.
#'           \item \code{TournamentSize()} in \code{SelectTournament}:
#'                 Specifies the size of the tournament. Per default: 2.
#'           \item \code{SelectionBias()} in \code{SelectLinearRank}.
#'                  This constant must be larger than 1.0 and usually
#'                  should be set at most to 2.0.
#'                  Increasing \code{SelectionBias()}
#'                  increases selection pressure.
#'           Beyond 2.0, there is the danger of premature convergence.
#'           }
#'
#' @section Performance Measurement:
#'
#' The file \code{Timer.R}:  Functions for timing and counting.
#'
#' The file \code{selectGeneBenchmark.R}: A benchmark of selection functions. 
#'
#' @section Interface Function Evaluation and Methods:
#'
#'              All evaluation functions must implement
#'              the following abstract interface:
#'
#'              \code{function name}(\code{gene}, \code{lF})
#'
#'              \strong{Parameters}
#'
#'              \itemize{
#'              \item \code{gene}  a gene. 
#'              \item \code{lF}   a local function list.
#'              }
#'
#'              \strong{Return Value}
#'
#'              A gene.
#' 
#' The file \code{evalGene.R} contains different function evaluation methods.
#' \enumerate{
#' \item \code{EvalGeneU} evaluates a gene unconditionally. (Default.)
#' \item \code{EvalGeneR} evaluates a gene unconditionally
#'        and allows the repair of the gene by the decoder.
#'
#' \item \code{EvalGeneDet} memoizes the evaluation of a gene in the 
#'              in the gene. Genes are evaluated only once.
#'              This leads to a performance improvement for
#'              deterministic functions.
#'
#' \item \code{EvalGeneStoch} computes an incremental average of 
#'            the value of a gene.
#'            The average converges to the true value as the number 
#'            of repeated evaluations of a gene increases.
#' } 
#'
#' @section Gene Representation:
#'
#' A gene is a named list:
#'   \itemize{
#'    \item $gene1      the gene.
#'    \item $fit        the fitness value of the gene
#'                      (for EvalGeneDet and EvalGeneU) or
#'                      the mean fitness (for stochastic functions
#'                      evaluated with EvalGeneStoch).
#'    \item $evaluated  has the gene been evaluated?
#'    \item $evalFail   has the evaluation of the gene failed?
#'    \item $var        the cumulative variance of the fitness 
#'                      of all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item $sigma      the standard deviation of the fitness of 
#'                      all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item $obs        the number of evaluations of a gene.
#'                      (For stochastic functions)
#'   }
#'
#' @section The Architecture of the xegaX-Packages:
#' 
#' The xegaX-packages are a family of R-packages which implement 
#' eXtended Evolutionary and Genetic Algorithms (xega).  
#' The architecture has 3 layers, 
#' namely the user interface layer,
#' the population layer, and the gene layer: 
#' 
#' \itemize{
#' \item
#' The user interface layer (package \code{xega}) 
#' provides a function call interface and configuration support
#' for several algorithms: genetic algorithms (sga), 
#' permutation-based genetic algorithms (sgPerm), 
#' derivation-free algorithms as e.g. differential evolution (sgde), 
#' grammar-based genetic programming (sgp) and grammatical evolution
#' (sge). 
#' \item
#' The population layer (package \code{xegaPopulation}) contains
#' population-related functionality as well as support for 
#' population statistics dependent adaptive mechanisms and parallelization.
#' \item 
#' The gene layer is split into a representation-independent and 
#' a representation-dependent part:
#' \enumerate{
#' \item 
#'  The representation-indendent part (package \code{xegaSelectGene})
#'  is responsible for variants of selection operators, evaluation 
#'  strategies for genes, and profiling and timing capabilities.        
#' \item 
#'  The representation-dependent part consists of the following packages: 
#' \itemize{
#' \item \code{xegaGaGene} for binary coded genetic algorithms.
#' \item \code{xegaPermGene} for permutation-based genetic algorithms.
#' \item \code{xegaDfGene} for derivation-free algorithms as e.g. 
#'                         differential evolution.
#' \item \code{xegaGpGene} for grammar-based genetic algorithms.
#' \item \code{xegaGeGene} for grammatical evolution algorithms.
#' }
#' The packages \code{xegaDerivationTrees} and \code{xegaBNF} support
#' the last two packages:
#' \code{xegaBNF} essentially provides a grammar compiler, and 
#' \code{xegaDerivationTrees} is an abstract data type for derivation trees.
#' }} 
#'
#' @family Package Description
#'
#' @name xegaSelectGene
#' @aliases xegaSelectGene
#' @docType package
#' @title Package xegaSelectGene.
#' @author Andreas Geyer-Schulz
#' @section Copyright: (c) 2023 Andreas Geyer-Schulz
#' @section License: MIT
#' @section URL: <https://github.com/ageyerschulz/xegaSelectGene>
#' @section Installation: From CRAN by \code{install.packages('xegaSelectGene')}
NULL
