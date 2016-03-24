#'
#'       minkowski.R
#' 
#'  Minkowski Sum
#'
#'  $Revision: 1.2 $ $Date: 2016/03/24 03:08:09 $


"%+%" <- MinkowskiSum <- local({

  poly2owin <- function(z) owin(poly=z, check=FALSE)

  sumconnected <- function(a, b, eps) {
    ## a and b are list(x,y) simply-connected polygons
    out <- polyclip::polyminkowski(a, b, x0=0, y0=0, eps=eps)
    if(length(out) == 1) return(out)
    ispos <- (sapply(out, Area.xypolygon) >= 0)
    if(sum(ispos) > 1) {
      browser()
      stop("Internal error: result of sumconnected is not simply connected",
           call.=FALSE)
    }
    return(out[ispos])
  }

  MinkowskiSum <- function(A, B) {
    A <- as.polygonal(A)
    B <- as.polygonal(B)
    ## determine common resolution for polyclip operations
    eps <- mean(c(sidelengths(Frame(A)), sidelengths(Frame(B))))/2^30
    ## separate into simply-connected pieces
    AA <- break.holes(A)$bdry
    BB <- break.holes(B)$bdry
    ## compute Minkowski sums of pieces
    pieces <- NULL
    for(b in BB) 
      pieces <- append(pieces, lapply(AA, sumconnected, b=b, eps=eps))
    ## form union in one step, to avoid artefacts
    result <- union.owin(solapply(pieces, poly2owin))
    return(result)
  }

  MinkowskiSum
})

