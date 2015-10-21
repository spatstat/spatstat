#
#   pairs.im.R
#
#   $Revision: 1.8 $   $Date: 2015/10/21 09:06:57 $
#

pairs.listof <- pairs.solist <- 
  pairs.im <- function(..., plot=TRUE) {
  argh <- list(...)
  cl <- match.call()
  ## unpack single argument which is a list of images
  if(length(argh) == 1) {
    arg1 <- argh[[1]]
    if(is.list(arg1) && all(unlist(lapply(arg1, is.im))))
      argh <- arg1
  }
  ## identify which arguments are images
  isim <- unlist(lapply(argh, is.im))
  nim <- sum(isim)
  if(nim == 0) 
    stop("No images provided")
  if(nim == 1) {
    ## one image: plot histogram
    h <- hist(..., plot=plot)
    return(invisible(h))
  }
  ## separate image arguments from others
  imlist <- argh[isim]
  rest   <- argh[!isim]
  ## determine image names for plotting
  imnames <- names(imlist)
  backupnames <- paste(cl)[c(FALSE, isim, FALSE)]
  if(length(backupnames) != nim)
    backupnames <- paste("V", seq_len(nim), sep="")
  if(length(imnames) != nim)
    imnames <- backupnames
  else if(any(needname <- !nzchar(imnames)))
    imnames[needname] <- backupnames[needname]
  ## extract pixel rasters and reconcile them
  imwins <- lapply(imlist, as.owin)
  names(imwins) <- NULL
  rasta    <- do.call("intersect.owin", imwins)
  ## extract image pixel values on common raster
  pixvals <- lapply(imlist, "[.im", i=rasta, raster=rasta, drop=TRUE)
  ## combine into data frame
  pixdf <- do.call("data.frame", pixvals)
  ## plot
  if(plot)
    do.call("pairs", resolve.defaults(list(x=pixdf),
                                      rest,
                                      list(labels=imnames, pch=".")))
  labels <- resolve.defaults(rest, list(labels=imnames))$labels
  colnames(pixdf) <- labels
  class(pixdf) <- c("plotpairsim", class(pixdf))
  return(invisible(pixdf))
}

plot.plotpairsim <- function(x, ...) {
  do.call("pairs.default",
          resolve.defaults(list(x=as.data.frame(x)),
                           list(...),
                           list(pch=".")))
  return(invisible(NULL))
}

print.plotpairsim <- function(x, ...) {
  cat("Object of class plotpairsim\n")
  cat(paste("contains pixel data for", commasep(sQuote(colnames(x))), "\n"))
  return(invisible(NULL))
}

panel.image <- function(x, y, ..., sigma=NULL) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  xx <- scaletointerval(x)
  yy <- scaletointerval(y)
  p <- ppp(xx, yy, window=square(1), check=FALSE)
  plot(density(p, sigma=sigma), add=TRUE, ...)
}

panel.contour <- function(x, y, ..., sigma=NULL) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  xx <- scaletointerval(x)
  yy <- scaletointerval(y)
  p <- ppp(xx, yy, window=square(1), check=FALSE)
  Z <- density(p, sigma=sigma)
  do.call("contour",
          resolve.defaults(list(x=Z, add=TRUE),
                           list(...),
                           list(drawlabels=FALSE)))
}

panel.histogram <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  do.call("rect",
          resolve.defaults(list(xleft   = breaks[-nB],
                                ybottom = 0,
                                xright  = breaks[-1],
                                ytop    = y),
                           list(...),
                           list(col="grey")))
}

  
  
  
