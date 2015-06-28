#'
#'    triangulate.R
#'
#'   Decompose a polygon into triangles
#'
#'    $Revision: 1.3 $  $Date: 2015/06/28 11:00:09 $
#'

triangulate.owin <- local({

  tricentre <- function(triangle) {
    as.numeric(lapply(vertices(triangle), mean))
  }
  
  triangulate.owin <- function(W) {
    stopifnot(is.owin(W))
    W <- as.polygonal(W)
    P <- as.ppp(vertices(W), W=Frame(W), check=FALSE)
    D <- delaunay(P)
    TD <- tiles(D)
    TC <- sapply(TD, tricentre)
    ok <- inside.owin(TC[1,], TC[2,], W)
    result <- tess(tiles=TD[ok], window=W, check=FALSE)
    return(result)
  }

  triangulate.owin
})
