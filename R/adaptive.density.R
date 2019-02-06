#'
#'   adaptive.density.R
#'
#'   $Revision: 1.1 $  $Date: 2019/02/06 03:22:51 $
#'

adaptive.density <- function(X, ..., method=c("voronoi", "kernel")) {
  method <- match.arg(method)
  result <- switch(method,
                   voronoi = densityVoronoi(X, ...),
                   kernel  = densityAdaptiveKernel(X, ...))
  return(result)
}
