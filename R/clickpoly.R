#
# clickpoly.R
#
#
# $Revision: 1.2 $  $Date: 2007/11/02 18:03:05 $
#
#
clickpoly <- function(add=FALSE, nv=NULL, np=1, ...) {
  if((!add) | dev.cur() == 1) {
    plot(0,0,type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), asp=1.0,
         axes=FALSE)
    rect(0,0,1,1)
  }
  gon <- list()
  stopifnot(np >= 1)
  for(i in 1:np) {
    if(np > 1)
      cat(paste(".... Polygon number", i, ".....\n"))
    if(!is.null(nv)) 
      cat(paste("click", nv, "times in window\n"))
    else
      cat(paste("to add points: click left mouse button in window\n",
                "      to exit: click middle mouse button\n",
                "[The last point should NOT repeat the first point]\n"))
    xy <- do.call("locator",
                  resolve.defaults(if(!is.null(nv)) list(n=nv) else list(),
                                   list(...),
                                   list(type="o")))
    if(area.xypolygon(xy) < 0)
      xy <- lapply(xy, rev)
    gon[[i]] <- xy
    plot(owin(poly=xy), add=TRUE)
  }
  result <- owin(poly=gon)
  plot(result, add=TRUE)
  return(result)
}

  
