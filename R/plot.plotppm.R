#
# plot.plotppm.R
#
# engine of plot method for ppm
#
# $Revision: 1.15 $  $Date: 2014/10/14 07:46:06 $
#
#

plot.plotppm <- function(x,data=NULL,trend=TRUE,cif=TRUE,se=TRUE,
                         pause=interactive(),
                         how=c("persp","image","contour"), ...,
                         pppargs=list())
{
  verifyclass(x,"plotppm")
  
  # determine main plotting actions
  superimposed <- !is.null(data)
  if(!missing(trend) && (trend & is.null(x[["trend"]])))
    stop("No trend to plot.\n")
  trend <- trend & !is.null(x[["trend"]])
  if(!missing(cif) && (cif & is.null(x[["cif"]])))
    stop("No cif to plot.\n")
  cif <- cif & !is.null(x[["cif"]])
  if(!missing(se) && (se & is.null(x[["se"]])))
    stop("No SE to plot.\n")
  se <- se & !is.null(x[["se"]])
  surftypes <- c("trend", "cif", "se")[c(trend, cif, se)]

  # marked point process?
  mrkvals <- attr(x,"mrkvals")
  marked <- (length(mrkvals) > 1)
  if(marked)
    data.marks <- marks(data)
  if(marked & superimposed) {
    data.types <- levels(data.marks)
    if(any(sort(data.types) != sort(mrkvals)))
      stop(paste("Data marks are different from mark",
                 "values for argument x.\n"))
  }

  # plotting style
  howmat <- outer(how, c("persp", "image", "contour"), "==")
  howmatch <- apply(howmat, 1, any)
  if (any(!howmatch)) 
    stop(paste("unrecognised option", how[!howmatch]))

  # start plotting
  if(pause)
    oldpar <- par(ask = TRUE)
  on.exit(if(pause) par(oldpar))

  
  for(ttt in surftypes) {
    xs <- x[[ttt]]
    for (i in seq_along(mrkvals)) {
      level <- mrkvals[i]
      main <- paste(if(ttt == "se") "Estimated" else "Fitted",
                    ttt, 
                    if(marked) paste("\n mark =", level) else NULL)
      for (style in how) {
        switch(style,
               persp = {
                 do.call("persp",
                         resolve.defaults(list(xs[[i]]),
                                          list(...), 
                                          spatstat.options("par.persp"),
                                          list(xlab="x", zlab=ttt, main=main)))
               },
               image = {
                 do.call("image",
                         resolve.defaults(list(xs[[i]]),
                                          list(...),
                                          list(main=main)))
                 if(superimposed) {
                   X <- if(marked) data[data.marks == level] else data
                   do.call(plot.ppp, append(list(x=X, add=TRUE), pppargs))
                 }
               },
               contour = {
                 do.call("contour",
                         resolve.defaults(list(xs[[i]]),
                                          list(...),
                                          list(main=main)))
                 if(superimposed) {
                   X <- if(marked) data[data.marks == level] else data
                   do.call(plot.ppp, append(list(x=X, add=TRUE), pppargs))
                 }
               },
               {
                 stop(paste("Unrecognised plot style", style))
               })
      }
    }
  }
  return(invisible())
}

print.plotppm <- function(x, ...) {
  verifyclass(x, "plotppm")
  trend   <- x$trend
  cif     <- x$cif
  mrkvals <- attr(x, "mrkvals")
  ntypes  <- length(mrkvals)
  unmarked <- (ntypes == 1 )
  cat(paste("Object of class", sQuote("plotppm"), "\n"))
  if(unmarked)
    cat("Computed for an unmarked point process\n")
  else {
    cat("Computed for a marked point process, with mark values:\n")
    print(mrkvals)
  }
  cat("Contains the following components:\n")
  if(!is.null(trend)) {
    cat("\n$trend:\tFitted trend.\n")
    if(unmarked) {
      cat("A list containing 1 image\n")
      print(trend[[1]], ...)
    } else {
      cat(paste("A list of", ntypes, "images\n"))
      cat("Typical details:\n")
      print(trend[[1]], ...)
    }
  }
  if(!is.null(cif)) {
    cat("\n$cif:\tFitted conditional intensity.\n")
    if(unmarked) {
      cat("A list containing 1 image\n")
      print(cif[[1]], ...)
    } else {
      cat(paste("A list of", ntypes, "images\n"))
      cat("Typical details:\n")
      print(cif[[1]], ...)
    }
  }
  invisible(NULL)
}

