#'
#'       minkowski.R
#' 
#'  Minkowski Sum and related operations
#'
#'  $Revision: 1.6 $ $Date: 2016/08/27 09:22:06 $


"%(+)%" <- MinkowskiSum <- local({

  MinkowskiSum <- function(A, B) {
    if(is.ppp(A)) return(UnionOfShifts(B, A))
    if(is.ppp(B)) return(UnionOfShifts(A, B))
    ## extract lists of simply-connected polygons
    AA <- simplepolygons(A)
    BB <- simplepolygons(B)
    ## determine common resolution for polyclip operations
    eps <- mean(c(sidelengths(Frame(A)), sidelengths(Frame(B))))/2^30
    ## compute Minkowski sums of pieces
    pieces <- NULL
    for(b in BB) 
      pieces <- append(pieces, lapply(AA, MinkSumConnected, b=b, eps=eps))
    ## form union in one step, to avoid artefacts
    result <- union.owin(solapply(pieces, poly2owin))
    return(result)
  }

  poly2owin <- function(z) owin(poly=z, check=FALSE)

  MinkSumConnected <- function(a, b, eps) {
    ## a and b are list(x,y) simply-connected polygons
    out <- polyclip::polyminkowski(a, b, x0=0, y0=0, eps=eps)
    if(length(out) == 1) return(out)
    ispos <- (sapply(out, Area.xypolygon) >= 0)
    if(sum(ispos) > 1) {
      stop("Internal error: result of sumconnected is not simply connected",
           call.=FALSE)
    }
    return(out[ispos])
  }

  simplepolygons <- function(A) {
    if(is.psp(A)) return(psp2poly(A))
    ## convert to owin, then polygonal
    A <- as.polygonal(A)
    ## separate into simply-connected pieces
    AA <- break.holes(A)$bdry
    return(AA)
  }
  
  ## handle segment patterns as well 
  psp2poly <- function(X) apply(as.matrix(X$ends), 1, seg2poly)

  seg2poly <- function(z) with(as.list(z), list(x=c(x0, x1, x0), y=c(y0,y1,y0)))

  ##
  UnionOfShifts <- function(X, V) {
    #' compute the union or superposition of copies of X by vectors in V
    v <- as.matrix(coords(V))
    n <- nrow(v)
    Y <- vector(mode="list", length=n)
    for(i in seq_len(n)) 
      Y[[i]] <- shift(X, v[i,])
    Y <- as.solist(Y)
    if(is.owin(X)) {
      Z <- union.owin(Y)
    } else {
      #' X is a pattern of objects in a window
      W <- MinkowskiSum(Window(X), Window(V))
      Z <- superimpose(Y, W=W)
    }
    return(Z)
  }

  MinkowskiSum
})

dilationAny <- function(A, B) { MinkowskiSum(A, reflect(B)) }

"%(-)%" <- erosionAny <- function(A, B) {
  D <- Frame(A)
  Dplus <- grow.rectangle(D, 0.1 * shortside(D))
  Ac <- complement.owin(A, Dplus)
  AcB <- MinkowskiSum(Ac, reflect(B))
  if(is.subset.owin(D, AcB))
    return(emptywindow(D))
  C <- complement.owin(AcB[Dplus], Dplus)[D]
  return(C)
}
