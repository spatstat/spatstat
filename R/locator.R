#'
#'     locator.R
#'
#'    $Revision: 1.1 $  $Date: 2017/01/07 09:23:51 $

spatstatLocator <- function(n, type=c("p","l","o","n"), ...) {
  #' remedy for failure of locator(type="p") in RStudio
  if(!identical(TRUE, dev.capabilities()$locator))
    stop("Sorry, this graphics device does not support the locator() function")
  # validate
  type <- match.arg(type)
  do.points <- type %in% c("p","o")
  do.lines <- type %in% c("l","o")
  argh <- list(...)
  pointsArgs <- c("cex", "col", "pch", "fg", "bg")
  segmentArgs <- graphicsPars("lines")
  # go
  res <- list(x=numeric(0), y = numeric(0))
  i <- 1L
  if(missing(n)) n <- Inf
  while(i<=n){
    tmp <- locator(n=1)
    if(is.null(tmp)) return(res)
    if(do.points)
      do.call.matched(points.default, append(tmp, argh), extrargs=pointsArgs)
    res$x <- c(res$x,tmp$x)
    res$y <- c(res$y,tmp$y)
    if(do.lines && i > 1L) {
      xy <- with(res, list(x0=x[i-1L], y0=y[i-1L], x1=x[i], y1=y[i]))
      do.call.matched(segments, append(xy, argh), extrargs=segmentArgs)
    }
    i <- i+1L
  }
  return(res)
}
  
