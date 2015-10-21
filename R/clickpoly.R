#
# clickpoly.R
#
#
# $Revision: 1.10 $  $Date: 2015/10/21 09:06:57 $
#
#

clickpoly <- function(add=FALSE, nv=NULL, np=1, ...) {
  if((!add) | dev.cur() == 1) {
    plot(0,0,type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), asp=1.0,
         axes=FALSE)
    rect(0,0,1,1)
  }
  spatstatLocator(0) ## check locator is enabled
  gon <- list()
  stopifnot(np >= 1)
  #
  for(i in 1:np) {
    if(np > 1)
      cat(paste(".... Polygon number", i, ".....\n"))
    if(!is.null(nv)) 
      cat(paste("click", nv, "times in window\n"))
    else
      cat(paste("to add points: click left mouse button in window\n",
                "      to exit: press ESC or click middle mouse button\n",
                "[The last point should NOT repeat the first point]\n"))
    xy <- do.call(spatstatLocator,
                  resolve.defaults(if(!is.null(nv)) list(n=nv) else list(),
                                   list(...),
                                   list(type="o")))
    if(Area.xypolygon(xy) < 0)
      xy <- lapply(xy, rev)
    gon[[i]] <- xy
    plotPolygonBdry(owin(poly=xy), ...)
  }
  result <- owin(poly=gon)
  plotPolygonBdry(result, ...)
  return(result)
}

clickbox <- function(add=TRUE, ...) {
  spatstatLocator(0) # check locator enabled
  cat("Click two corners of a box\n")
  if(!add) plot(owin(), main="Click two corners of a box") 
  a <- try(spatstatLocator(1), silent=TRUE)
  if(inherits(a, "try-error")) {
    ## add=TRUE but there is no current plot
    plot.new()
    a <- spatstatLocator(1, ...)
  }
  abline(v=a$x)
  abline(h=a$y)
  b <- spatstatLocator(1, ...)
  abline(v=b$x)
  abline(h=b$y)
  ab <- concatxy(a, b)
  result <- owin(range(ab$x), range(ab$y))
  plotPolygonBdry(result, ...)
  return(result)
}

plotPolygonBdry <- function(x, ...) {
  # filter appropriate arguments
  argh <- list(...)
  polyPars <- union(graphicsPars("lines"), graphicsPars("owin"))
  polyargs <- argh[names(argh) %in% polyPars]
  # change 'col' to 'border'
  nama <- names(polyargs)
  if(any(nama == "col") && !any(nama == "border"))
    names(polyargs)[nama == "col"] <- "border"
  # plot
  do.call(plot.owin,
          append(list(x=x, add=TRUE), polyargs))
}
