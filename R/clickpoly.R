#
# clickpoly.R
#
#
# $Revision: 1.5 $  $Date: 2014/09/27 09:45:44 $
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
    if(Area.xypolygon(xy) < 0)
      xy <- lapply(xy, rev)
    gon[[i]] <- xy
    plot(owin(poly=xy), add=TRUE)
  }
  result <- owin(poly=gon)
  plot(result, add=TRUE)
  return(result)
}
  
clickbox <- function(add=TRUE) {
  cat("Click two corners of a box\n")
  if(!add) plot(owin(), main="Click two corners of a box") 
  a <- try(locator(1), silent=TRUE)
  if(inherits(a, "try-error")) {
    ## add=TRUE but there is no current plot
    plot.new()
    a <- locator(1)
  }
  abline(v=a$x)
  abline(h=a$y)
  b <- locator(1)
  abline(v=b$x)
  abline(h=b$y)
  ab <- concatxy(a, b)
  return(owin(range(ab$x), range(ab$y)))
}
