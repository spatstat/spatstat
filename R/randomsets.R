#'
#'  randomsets.R
#'
#'  Generation of random sets
#'
#'  $Revision: 1.2 $ $Date: 2019/08/16 07:53:05 $


rthinclumps <- function(W, p, ...) {
  check.1.real(p)
  if(badprobability(p, TRUE))
    stop("p must be a valid probability between 0 and 1", call.=FALSE)
  if(!(is.im(W) || is.owin(W)))
    stop("W should be a window or pixel image", call.=FALSE)
  clumps <- connected(W, ...)
  keep <- (runif(length(levels(clumps))) < p)
  retained <- eval.im(keep[clumps])
  return(solutionset(retained))
}


  
