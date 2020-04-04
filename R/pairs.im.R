#
#   pairs.im.R
#
#   $Revision: 1.17 $   $Date: 2020/04/04 04:22:50 $
#

pairs.listof <- pairs.solist <- function(..., plot=TRUE) {
  argh <- expandSpecialLists(list(...), special=c("solist", "listof"))
  haslines <- any(sapply(argh, inherits, what="linim"))
  names(argh) <- good.names(names(argh), "V", seq_along(argh))
  if(haslines) {
    do.call(pairs.linim, append(argh, list(plot=plot)))
  } else {
    do.call(pairs.im, append(argh, list(plot=plot)))
  }
}

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
  ## separate image arguments from others
  imlist <- argh[isim]
  rest   <- argh[!isim]
  ## determine image names for plotting
  imnames <- argh$labels %orifnull% names(imlist)
  if(length(imnames) != nim || !all(nzchar(imnames))) {
    #' names not given explicitly
    callednames <- paste(cl)[c(FALSE, isim, FALSE)]
    backupnames <- paste0("V", seq_len(nim))
    if(length(callednames) != nim) {
      callednames <- backupnames
    } else if(any(toolong <- (nchar(callednames) > 15))) {
      callednames[toolong] <- backupnames[toolong]
    }
    imnames <- good.names(imnames, good.names(callednames, backupnames))
  }
  ## 
  if(nim == 1) {
    ## one image: plot histogram
    Z <- imlist[[1L]]
    xname <- imnames[1L]
    do.call(hist,
            resolve.defaults(list(x=Z, plot=plot),
                             rest, 
                             list(xlab=xname,
                                  main=paste("Histogram of", xname))))
    ## save pixel values
    pixvals <- list(Z[])
    names(pixvals) <- xname
  } else {
    ## extract pixel rasters and reconcile them
    imwins <- lapply(imlist, as.owin)
    names(imwins) <- NULL
    rasta    <- do.call(intersect.owin, imwins)
    ## extract image pixel values on common raster
    pixvals <- lapply(imlist, "[.im", i=rasta, raster=rasta, drop=TRUE)
  }
  ## combine into data frame
  pixdf <- do.call(data.frame, pixvals)
  ## pairs plot
  if(plot && nim > 1)
    do.call(pairs, resolve.defaults(list(x=pixdf),
                                      rest,
                                      list(labels=imnames, pch=".")))
  labels <- resolve.defaults(rest, list(labels=imnames))$labels
  colnames(pixdf) <- labels
  class(pixdf) <- c("plotpairsim", class(pixdf))
  return(invisible(pixdf))
}

plot.plotpairsim <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  x <- as.data.frame(x)
  if(ncol(x) == 1) {
    do.call(hist.default,
            resolve.defaults(list(x=x[,1]),
                             list(...),
                             list(main=xname, xlab=xname)))
  } else {
    do.call(pairs.default,
            resolve.defaults(list(x=x),
                             list(...),
                             list(pch=".")))
  }
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
  do.call(contour,
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
  do.call(rect,
          resolve.defaults(list(xleft   = breaks[-nB],
                                ybottom = 0,
                                xright  = breaks[-1],
                                ytop    = y),
                           list(...),
                           list(col="grey")))
}

  
  
  
