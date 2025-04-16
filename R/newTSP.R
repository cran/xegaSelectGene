#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaPermGene
#

#' Generate a TSP problem environment
#'
#' @description \code{newTSP()} generates the problem environment 
#'              for a traveling salesman problem (TSP).
#'
#' @details \code{newTSP()} provides several local permutation 
#'        improvement heuristics: 
#'        a greedy path of length k starting
#'        from city i, 
#'        the best greedy path of length k, 
#'        a random 2-Opt-move,
#'        and a sequence of random 2-Opt moves.
#'        They help to find bounds for the TSP 
#'        or to implement special purpose mutation operators. 
#'
#' @param D        A \code{n} x \code{n} distance matrix. 
#' @param Cities   The names of the cities.
#' @param Name     The name of the problem environment.
#' @param Solution Solution of problem (if known). 
#'        Default: \code{NA}
#' @param Path  Optimal permutation of cities (if known). As integer vector.
#'        Default: \code{NA}.
#'
#' @family Problem Environments
#'
#' @return A problem environment for the TSP.
#'  \enumerate{
#'  \item 
#'  \code{$name()}             a string with the name of the environment
#'  \item 
#'  \code{$cities()} a vector length \code{n} of city names.           
#'  \item 
#'  \code{$dist()}  the \code{n} x \code{n} distance matrix 
#'                            between \code{n} cities.
#'  \item 
#'  \code{$genelength()} the size of the permutation \code{n}.
#'                            E.g. for a TSP: the number of cities. 
#'  \item 
#'  \code{$f(permutation, gene, lF, tour=TRUE)} 
#'       the fitness function of the TSP. If \code{tour==FALSE},
#'       the path length is computed 
#'       (without the cost from city n to city 1).
#' 
#'       With a permutation of size \code{n} as argument. 
#'  \item 
#'  \code{$show(p)}    shows tour through the cities in 
#'                     path \code{p} with its cost. 
#'  \item 
#'  \code{$greedy(startPosition, k)} 
#'                   computes a \code{k}-step greedy minimal cost path 
#'                             beginning at the city \code{start}.
#'                             For  \code{k+1=n} the greedy solution gives
#'                             an upper bound for the TSP.
#'  \item 
#'  \code{$kBestGreedy(k, tour=TRUE)} computes the best greedy 
#'                             subtour with \code{k+1} cities.
#'                     For \code{tour=FALSE}, the best greedy
#'                             subpath with \code{k+1} cities is 
#'                             computed.
#'
#'  \item 
#'  \code{$rnd2Opt(permutation, maxTries=5)} returns a permutation 
#'                             improved by a single random 2-Opt-move  
#'                             after at most \code{maxTries=5} attempts.
#'
#'  \item 
#'  \code{$LinKernighan(permutation, maxTries=5, show=FALSE)} returns 
#'                         the best permutation found after
#'                         several random 2-Opt-moves
#'                         with at most \code{maxTries=5} attempts.
#'                         The loop stops after the first 2-Opt-move
#'                         which does not improve the solution.
#'
#'  \item 
#'  \code{$solution()}         known optimal solution.
#'  \item 
#'  \code{$path()}             known optimal round trip.
#'  \item 
#'  \code{$max()}            \code{FALSE}.
#'  \item {$globalOptimum()}  a named list with the following elements:
#'     \itemize{
#'     \item \code{$param}    the optimal permutation.
#'     \item \code{$value}    the known optimal solution.
#'     \item \code{$is.minimum}  \code{TRUE}. 
#'     }    
#'  }
#'
#' @examples
## a<-read.table("lau15distAGS.txt")
#' a<-matrix(0, nrow=15, ncol=15)
#' a[1,]<- c(0, 29, 82, 46, 68, 52, 72, 42, 51,  55,  29,  74,  23,  72,  46)
#' a[2,]<- c(29,  0, 55, 46, 42, 43, 43, 23, 23,  31,  41,  51,  11,  52,  21)
#' a[3,]<- c(82, 55,  0, 68, 46, 55, 23, 43, 41,  29,  79,  21,  64,  31,  51)
#' a[4,]<-c(46, 46, 68,  0, 82, 15, 72, 31, 62,  42,  21,  51,  51,  43,  64)
#' a[5,]<-c(68, 42, 46, 82,  0, 74, 23, 52, 21,  46,  82,  58,  46,  65,  23)
#' a[6,]<-c(52, 43, 55, 15, 74,  0, 61, 23, 55,  31,  33,  37,  51,  29,  59)
#' a[7,]<-c(72, 43, 23, 72, 23, 61,  0, 42, 23,  31,  77,  37,  51,  46,  33)
#' a[8,]<-c(42, 23, 43, 31, 52, 23, 42,  0, 33,  15,  37,  33,  33,  31,  37)
#' a[9,]<-c(51, 23, 41, 62, 21, 55, 23, 33,  0,  29,  62,  46,  29,  51,  11)
#' a[10,]<-c(55, 31, 29, 42, 46, 31, 31, 15, 29,  0,  51,  21,  41,  23,  37)
#' a[11,]<-c(29, 41, 79, 21, 82, 33, 77, 37, 62,  51,   0,  65,  42,  59,  61)
#' a[12,]<-c(74, 51, 21, 51, 58, 37, 37, 33, 46,  21,  65,   0,  61,  11,  55)
#' a[13,]<-c(23, 11, 64, 51, 46, 51, 51, 33, 29,  41,  42,  61,   0,  62,  23)
#' a[14,]<-c(72, 52, 31, 43, 65, 29, 46, 31, 51,  23,  59,  11,  62,   0,  59)
#' a[15,]<-c(46, 21, 51, 64, 23, 59, 33, 37, 11,  37,  61,  55,  23,  59,   0)
#' lau15<-newTSP(a, Name="lau15")
#' lau15$name()
#' lau15$genelength()
#' b<-sample(1:15, 15, FALSE)
#' lau15$f(b)
#' lau15$f(b, tour=TRUE)
#' lau15$show(b)
#' lau15$greedy(1, 14)
#' lau15$greedy(1, 1)
#'
#' @export
newTSP<-function(D, Name, Cities=NA, Solution=NA, Path=NA)
{ d<-dim(D)
if (!length(d)==2) stop("n times n matrix expected")  
if (!d[1]==d[2]) stop("n times n matrix expected")  
self<-list()
# constant functions
self$name<-parm(Name)
self$genelength<-parm(d[1])
self$dist<-parm(D)
#### 
if (all(is.na(Cities))) {cit<-1:d[1]} 
else {if (!length(Cities)==d[1]) {stop("List of n cities expected")} 
      else {cit<-Cities}}
self$cities<-parm(cit)
if (!all(is.na(Path)))  
   {if (!length(Path)==d[1]) {stop("Path of n cities expected")}}
if (!all(is.na(Path)))  
   { if (!all(Path %in% (1:d[1]))) {stop("Permutation of n integers expected")}}
self$solution<-parm(Solution)
self$path<-parm(Path)
self$max<-function() {return(FALSE)}
self$globalOptimum<-function()
  {l<-list()
  l$param<-Path; l$value<-Solution;l$is.minimum<-TRUE;return(l)}
# f
self$f<-function(permutation, gene=0, lF=0, tour=TRUE)
{ cost<-0
  l<-length(permutation)-1
  for (i in 1:l) 
     { cost<- cost+self$dist()[permutation[i], permutation[i+1]]}
  if (tour==TRUE) {cost<-cost+self$dist()[permutation[l+1], permutation[1]]}
  return(cost)}
# show. p is a path.
self$show<-function(p)
{ l<-length(p)-1
  pl<-0
  for (i in 1:l)
  {d<-self$dist()[p[i], p[i+1]] 
   pl<-pl+d
   cat(i, "From:", self$cities()[p[i]], " to ", self$cities()[p[i+1]], 
       " Distance: ", d, " ", pl, "\n")}
  d<-self$dist()[p[l+1], p[1]] 
   pl<-pl+d
   cat((l+1), "From:", self$cities()[p[l+1]], " to ", self$cities()[p[1]], 
       " Distance: ", d, " ", pl, "\n") }
# greedy
self$greedy<-function(startPosition, k)
{ # local functions
   without<-function(set, element) {set[!set==element]}
  # v a vector
   findMinIndex<-function(indexSet, v)
    {v<-indexSet[v==min(v)]
     return(v[sample(1:length(v),1)]) }

   nextPosition<-startPosition
   path<-as.vector(startPosition)
   indexSet<-without(1:self$genelength(), nextPosition)
   for (i in 1:k)
   {nextPosition<-findMinIndex(indexSet, self$dist()[nextPosition, indexSet])
    path<-c(path, nextPosition)
    indexSet<-without(indexSet, nextPosition)}
  return(path)}

self$kBestGreedy<-function(k, tour=TRUE)
{ l<-self$genelength()
  best<-self$greedy(1, k)
  costBest<-self$f(best, tour)
  for (i in 2:l)
  {
  new<-self$greedy(i, k)
  costNew<-self$f(new, tour)
  if (costNew<costBest) {best<-new; costBest<-costNew}
  }

  return(best)
}

self$rnd2Opt<-function(permutation, maxTries=5)
{

randomSplit<-function(l)
{ kpos<-sample(1:l, 2, replace=FALSE)
        if (kpos[1]>kpos[2]) {kpos[c(2, 1)]<-kpos}
        if ((kpos[2]==l) & (kpos[1]==1)) {kpos[2]<-l-1}
        if ((kpos[2]-kpos[1])==1)
           {if (kpos[2]==l) {kpos[1]<-kpos[1] -1} else {kpos[2]<-kpos[2]+1}}
        x1<-1:kpos[1]
        if (kpos[2]<l) {x2<-(kpos[2]+1):l} else {x2<-rep(0,0)}
        y<-(kpos[1]+1):kpos[2]
        return(c(x2, x1, y[length(y):1])) }

	cost1<-self$f(permutation)
	l<-length(permutation)
	for (i in 1:maxTries)
	{ newpermutation<-permutation[randomSplit(l)]
	cost2<-self$f(newpermutation)
	if (cost1>cost2) {return(newpermutation)}
	if (i==maxTries) {break} }
	return(permutation)
}

self$LinKernighan<-function(permutation, maxTries=5, show=FALSE)
{
	epsilon<-parm(0.000001)
	newpermutation<-permutation; i<-1
	repeat
	{ c1<-self$f(newpermutation); i<-i+1
	  newpermutation<-self$rnd2Opt(newpermutation, maxTries)
	  c2<-self$f(newpermutation)
	  if (show) {
	    cat(i,"p: ", newpermutation, 
		"c1:", c1, "c2:", c2, "diff:", (c1-c2), "\n")}
	  if (abs(c1-c2)<epsilon()) {break}
	}
	return(newpermutation)
}

# The next 8 statements are needed to evaluate the promises
# and to force a static binding despite lazy evaluation of
# arguments.
a<-self$name()
a<-self$genelength()
a<-self$dist()
a<-self$cities()
a<-self$solution()
a<-self$path()
a<-self$max()
a<-self$globalOptimum()
return(self) }

#
# Data for lau15
#
a<-matrix(0, nrow=15, ncol=15)

a[1,]<- c(0, 29, 82, 46, 68, 52, 72, 42, 51,  55,  29,  74,  23,  72,  46)
a[2,]<- c(29,  0, 55, 46, 42, 43, 43, 23, 23,  31,  41,  51,  11,  52,  21)
a[3,]<- c(82, 55,  0, 68, 46, 55, 23, 43, 41,  29,  79,  21,  64,  31,  51)
a[4,]<-c(46, 46, 68,  0, 82, 15, 72, 31, 62,  42,  21,  51,  51,  43,  64)
a[5,]<-c(68, 42, 46, 82,  0, 74, 23, 52, 21,  46,  82,  58,  46,  65,  23)
a[6,]<-c(52, 43, 55, 15, 74,  0, 61, 23, 55,  31,  33,  37,  51,  29,  59)
a[7,]<-c(72, 43, 23, 72, 23, 61,  0, 42, 23,  31,  77,  37,  51,  46,  33)
a[8,]<-c(42, 23, 43, 31, 52, 23, 42,  0, 33,  15,  37,  33,  33,  31,  37)
a[9,]<-c(51, 23, 41, 62, 21, 55, 23, 33,  0,  29,  62,  46,  29,  51,  11)
a[10,]<-c(55, 31, 29, 42, 46, 31, 31, 15, 29,  0,  51,  21,  41,  23,  37)
a[11,]<-c(29, 41, 79, 21, 82, 33, 77, 37, 62,  51,   0,  65,  42,  59,  61)
a[12,]<-c(74, 51, 21, 51, 58, 37, 37, 33, 46,  21,  65,   0,  61,  11,  55)
a[13,]<-c(23, 11, 64, 51, 46, 51, 51, 33, 29,  41,  42,  61,   0,  62,  23)
a[14,]<-c(72, 52, 31, 43, 65, 29, 46, 31, 51,  23,  59,  11,  62,   0,  59)
a[15,]<-c(46, 21, 51, 64, 23, 59, 33, 37, 11,  37,  61,  55,  23,  59,   0)

path<-c(1, 13, 2, 15, 9, 5, 7, 3, 12, 14, 10, 8, 6, 4, 11)

#' The problem environment \code{lau15} for a traveling salesman problem.
#'
#' @description 
#' 15 abstract cities for which a traveling salesman solution is sought.
#' Solution: A path with a length of 291.
#' 
#' The problem environment \code{lau15} is a list with the following functions: 
#'
#' \enumerate{
#' \item \code{lau15$name()}:  \code{"lau15"}, the name of the TSP  
#'                             problem environment.
#' \item \code{lau15$genelength()}: 15, the number of cities on the round trip.
#' \item \code{lau15$dist()}: The distance matrix of the problem.
#' \item \code{lau15$cities()}: A list of city names or the vector 
#'                         \code{1:lau15$genelength()}.
#' \item \code{lau15$f (permutation, gene = 0, lF = 0, tour = TRUE)}: 
#'                     The fitness function. Computes the roundtrip 
#'                     for permutation of cities.
#' \item \code{lau15$solution()}: 291, the known optimal solution of lau15.
#' \item \code{lau15$path()}: The permutation for the optimal roundtrip.
#' \item \code{lau15$show(p)}: Prints the roundtrip \code{p}.
#' \item \code{lau15$greedy(startposition, k)}: Computes a path of length
#'                          \code{k} starting at \code{startposition} 
#'                          by choosing the nearest city.
#' \item \code{lau15$kBestGreedy(k, tour=TRUE)}:
#'                          Computes the best greedy path/tour with 
#'                          k cities. 
#' \item \code{lau15$rnd2Opt(permutation, maxtries=5)}:
#'                          Tries to find a better permutation by 
#'                          at most 5 random 2-opt heuristics.
#' \item \code{lau15$LinKernighan(permutation, maxtries=5)}:
#'                A randomized Lin-Kernigan heuristic implemented 
#'                as a sequence of randomized 2-opt moves.
#' }
#'
#' @references 
#' Lau, H. T. (1986):
#' \emph{Combinatorial Heuristic Algorithms in FORTRAN}.
#' Springer, 1986. p. 61. <doi:10.1007/978-3-642-61649-5>
#'
#' @family Problem Environments
#'
#' @export
lau15<-newTSP(a, Name="lau15", Solution= 291, Path=path)

