#'
#'    triangulate.R
#'
#'   Decompose a polygon into triangles
#'
#'    $Revision: 1.4 $  $Date: 2015/11/21 11:13:00 $
#'

triangulate.owin <- local({

  is.triangle <- function(p) {
    return((length(p$bdry) == 1) && (length(p$bdry[[1]]$x) == 3))
  }

  triangulate.owin <- function(W) {
    stopifnot(is.owin(W))
    W <- as.polygonal(W, repair=TRUE)
    P <- as.ppp(vertices(W), W=Frame(W), check=FALSE)
    D <- delaunay(P)
    V <- intersect.tess(W, D)
    Candidates <- tiles(V)
    istri <- sapply(Candidates, is.triangle)
    Accepted <- Candidates[istri]
    if(any(!istri)) {
      # recurse
      Worries <- unname(Candidates[!istri])
      Fixed <- lapply(Worries, triangulate.owin)
      Fixed <- do.call(c, lapply(Fixed, tiles))
      Accepted <- append(Accepted, Fixed)
    }
    result <- tess(tiles=Accepted, window=W)
    return(result)
  }

  triangulate.owin
})
