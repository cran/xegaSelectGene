
#' Generates a problem environment for the XOR problem.
#'
#' @return The problem environment for the XOR problem with 
#'         \itemize{
#' \item \code{$name}: 
#'         \code{"envXOR"}, the name of the problem environment.
#' \item \code{$buildtest(expr)}: 
#'                               The function which builds the environment 
#'                               for evaluating the expression 
#'                               by binding the variables to the parameters.
#' \item \code{$TestCases}: The truthtable of the XOR function.
#' \item \code{$f(expr, gene=NULL, lF=NULL)}: The fitness function.
#'                              \code{expr} is the string with the 
#'                              logical expression to be evaluated. 
#'         }
#'
#' @family Problem Environments
#'
#' @return The problem environment.
#' @examples
#' envXOR<-newEnvXOR()
#' envXOR$name()
#' a2<-"OR(OR(D1, D2), (AND(NOT(D1), NOT(D2))))"
#' a3<-"OR(OR(D1, D2), AND(D1, D2))"
#' a4<-"AND(OR(D1,D2),NOT(AND(D1,D2)))"
#' gp4<-"(AND(AND(OR(D2,D1),NOT(AND(D1,D2))),(OR(D2,D1))))"
#' envXOR$f(a2)
#' envXOR$f(a3)
#' envXOR$f(a4)
#' envXOR$f(gp4)
#'
## @importFrom xegaDerivationTrees treeLeaves
#' @export
newEnvXOR<-function()
{
self<-list()
self$name<-function() {"envXOR"}
self$buildTest<-function(expr) {
        f<-paste("function(v) {
        AND<-function(x,y){return(x & y)}
        NAND<-function(x,y){return(!(x & y))}
        OR<-function(x,y){return(x|y)}
        NOT<-function(x){return(!x)}
        D1<-v[1]
        D2<-v[2]
        return(", expr, ")}", sep="")
        return(eval(parse(text=f)))
}
self$TestCases<-matrix(c(0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0),
                nrow=4, ncol=3, byrow=TRUE)
self$f<-function(expr, gene=NULL, lF=NULL)
{ test<-self$buildTest(expr)
s<-0
for (i in 1:nrow(self$TestCases))
{ s<-s+
(self$TestCases[i,ncol(self$TestCases)]==test(self$TestCases[i,]))
}
if (identical(gene, NULL)) {return(s)}
  # b<-xegaDerivationTrees::treeLeaves(gene$gene1, lF$Grammar$ST)
  # s<-(s+(1/(b^2)))
return(s)
}
a<-self$name()
a<-self$buildTest("D1")
a<-self$TestCases
return(self)
}

#' The problem environment \code{envXOR} for programming the XOR function
#' either by grammar-based genetic programming or grammatical evolution.
#' 
#' @description 
#' The problem environment \code{envXOR} is a list with the following elements:
#' \itemize{
#' \item \code{envXOR$name}: 
#'         \code{"envXOR"}, the name of the problem environment.
#' \item \code{envXOR$buildtest(expr)}: 
#'                               The function which builds the environment 
#'                               for evaluating the expression with 
#'                               a binding of the variables to the parameters.
#' \item \code{envXOR$TestCases}: The truth table of the XOR function.
#' \item \code{envXOR$f(expr, gene=NULL, lF=NULL)}: The fitness function.
#'                              \code{expr} is the string with the 
#'                              logical expression to be evaluated. 

#' }
#'
#' @family Problem Environments
#'
#' @export
envXOR<-newEnvXOR()

