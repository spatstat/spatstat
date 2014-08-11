#
#  is.subset.owin.R
#
#  $Revision: 1.10 $   $Date: 2013/11/01 06:49:39 $
#
#  Determine whether a window is a subset of another window
#
#  is.subset.owin()
#

is.subset.owin <- local({
  
  is.subset.owin <- function(A, B) {
    A <- as.owin(A)
    B <- as.owin(B)

    if(identical(A, B))
      return(TRUE)

    A <- rescue.rectangle(A)
    B <- rescue.rectangle(B)
  
    if(is.rectangle(B)) {
      # Some cases can be resolved using convexity of B
    
      # (1) A is also a rectangle
      if(is.rectangle(A)) {
        xx <- A$xrange[c(1,2,2,1)]
        yy <- A$yrange[c(1,1,2,2)]
        ok <- inside.owin(xx, yy, B)
        return(all(ok))
      } 
      # (2) A is polygonal
      # Then A is a subset of B iff,
      # for every constituent polygon of A with positive sign,
      # the vertices are all in B
      if(is.polygonal(A)) {
        ok <- unlist(lapply(A$bdry, okpolygon, B=B))
        return(all(ok))
      }
      # (3) Feeling lucky
      # Test whether the bounding box of A is a subset of B
      # Then a fortiori, A is a subset of B
      AA <- bounding.box(A)
      if(is.subset.owin(AA, B))
        return(TRUE)
    }

    if(!is.mask(A) && !is.mask(B)) {
      # rectangles or polygonal domains
      if(!all(inside.owin(vertices(A), , B)))
        return(FALSE)
      # all vertices of A are inside B.
      if(is.convex(B))
        return(TRUE)
      A <- as.polygonal(A)
      B <- as.polygonal(B)
      if(length(B$bdry) == 1 && length(A$bdry) == 1) {
        # two simply-connected sets 
        # check for boundary crossings
        bx <- crossing.psp(as.psp(A), as.psp(B))
        return(npoints(bx) == 0)
      } else {
        # compare area of intersection with area of A
        return(overlap.owin(A,B) >= area.owin(A))
      }
    }
  
   # Discretise
    a <- as.mask(A)
    b <- as.mask(B)
    xx <- as.vector(raster.x(a)[a$m])
    yy <- as.vector(raster.y(a)[a$m])
    ok <- inside.owin(xx, yy, b)
    return(all(ok))
    
  }

  okpolygon <- function(a, B) {
    if(area.xypolygon(a) < 0) return(TRUE)
    ok <- inside.owin(a$x, a$y, B)
    return(all(ok))
  }

  is.subset.owin
})
