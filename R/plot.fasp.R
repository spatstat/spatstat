#
#   plot.fasp.R
#
#   $Revision: 1.29 $   $Date: 2016/02/11 10:17:12 $
#
plot.fasp <- function(x, formule=NULL, ..., subset=NULL,
                      title=NULL, banner=TRUE,
                      transpose=FALSE,
                      samex=FALSE, samey=FALSE,
                      mar.panel=NULL,
                      outerlabels=TRUE, cex.outerlabels=1.25,
                      legend=FALSE) {
  # plot dimensions
  which <- x$which
  if(transpose)
    which <- t(which)
  nrows  <- nrow(which)
  ncols  <- ncol(which)
  
# Determine the overall title of the plot
  if(banner) {
    if(!is.null(title)) overall <- title
    else if(!is.null(x$title)) overall <- x$title
    else {
      if(prod(dim(which)) > 1)
        overall <- "Array of diagnostic functions"
      else
        overall <- "Diagnostic function"
      if(is.null(x$dataname)) overall <- paste(overall,".",sep="")
      else overall <- paste(overall," for ",x$dataname,".",sep="")
    }
    if(length(overall) > 1)
      overall <- paste(overall, collapse="\n")
    nlines <-
      if(!is.character(overall)) 1 else length(unlist(strsplit(overall, "\n")))
  } 

# If no formula is given, look for a default formula in x:
  defaultplot <- is.null(formule)
  if(defaultplot && !is.null(x$default.formula))
    formule <- x$default.formula

  if(!is.null(formule)) {
    # ensure formulae are given as character strings.
    formule <- FormatFaspFormulae(formule, "formule")
    # Number of formulae should match number of functions.
    nf <- length(formule)
    nfun <- length(x$fns)
    if(nf == 1 && nfun > 1)
      formule <- rep.int(formule, nfun)
    else if(nf != nfun)
      stop(paste("Wrong number of entries in", sQuote("formule")))
  }
  
# Check on the length of the subset argument.
  ns <- length(subset)
  if(ns > 1) {
    if(ns != length(x$fns))
      stop("Wrong number of entries in subset argument.\n")
    msub <- TRUE
  } else msub <- FALSE

# compute common x, y axis limits for all plots ?
  xlim <- ylim <- NULL
  if(samex || samey) {
    cat("Computing limits\n")
    # call plot.fv to determine plot limits of each panel
    for(i in 1:nrows) {
      for(j in 1:ncols) {
        k <- which[i,j]
        if(!is.na(k)) {
          fun <- as.fv(x$fns[[k]])
          fmla <- if(!defaultplot) formule[k] else NULL
          sub <- if(msub) subset[[k]] else subset
          lims <- plot(fun, fmla, subset=sub, limitsonly=TRUE)
          # update the limits
          if(samex) xlim <- range(xlim, lims$xlim)
          if(samey) ylim <- range(ylim, lims$ylim)
        }
      }
    }
  } 

#############################################################  
# Set up the plot layout
  n <- nrows * ncols
# panels 1..n = plot panels
  codes <- matrix(seq_len(n), byrow=TRUE, ncol=ncols, nrow=nrows)
  heights <- rep.int(1, nrows)
  widths  <- rep.int(1, ncols)
# annotation as chosen
  if(outerlabels) {
    # column headings
    colhead.codes <- max(codes) + (1:ncols)
    colhead.height <- 0.2
    codes <- rbind(colhead.codes, codes)
    heights <- c(colhead.height, heights)
    # row headings
    rowhead.codes <- max(codes) + (1:nrows)
    rowhead.width <- 0.2
    codes <- cbind(c(0,rowhead.codes), codes)
    widths <- c(rowhead.width, widths)
  }
  if(banner) {
    # overall banner
    top.code <- max(codes) + 1
    top.height <- 0.1 * (1+nlines)
    codes <- rbind(top.code, codes)
    heights <- c(top.height, heights)
  }

# declare layout  
  layout(codes, widths=widths, heights=heights)

############################################################  
# Plot the function panels 
#
# determine annotation
  colNames <- colnames(which)
  rowNames <- rownames(which)
  nrc <- max(nrows, ncols)
  ann.def <- par("ann") && (nrc <= 3)
# determine margin around each panel
  if(is.null(mar.panel)) 
    mar.panel <- if(nrc > 3 && outerlabels) rep.int(1/nrc, 4) else par("mar")
  opa <- par(mar=mar.panel, xpd=TRUE)
#
# plot each function  
  for(i in 1:nrows) {
    for(j in 1:ncols) {
      k <- which[i,j]
      if(is.na(k)) plot(0,0,type='n',xlim=c(0,1),
                        ylim=c(0,1),axes=FALSE,xlab='',ylab='', ...)
      else {
        fun <- as.fv(x$fns[[k]])
        fmla <- if(!defaultplot) formule[k] else NULL
        sub <- if(msub) subset[[k]] else subset
        main <- if(outerlabels) "" else
                if(nrows == 1) colNames[j] else
                if(ncols == 1) rowNames[i] else 
                paren(paste(rowNames[i], colNames[j], sep=","))
        do.call(plot,
                resolve.defaults(list(x=fun, fmla=fmla, subset=sub),
                                 list(...),
                                 list(xlim=xlim, ylim=ylim,
                                      main=main, legend=legend),
                                 list(ann=ann.def, axes=ann.def,
                                      frame.plot=TRUE)))
      }
    }
  }
############################################################
#
# Annotation as selected
  if(outerlabels) {
    par(mar=rep.int(0,4), xpd=TRUE)
    # Plot the column headers
    for(j in 1:ncols) {
      plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
           xlim=c(-1,1),ylim=c(-1,1))
      text(0,0,colNames[j], cex=cex.outerlabels)
    }
    # Plot the row labels
    for(i in 1:nrows) {
      plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
           xlim=c(-1,1),ylim=c(-1,1))
      text(0,0,rowNames[i], srt=90, cex=cex.outerlabels)
    }
  }
  if(banner) {
    par(mar=rep.int(0,4), xpd=TRUE)
    # plot the banner
    plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
         xlim=c(-1,1),ylim=c(-1,1))
    cex <- resolve.defaults(list(...), list(cex.title=2))$cex.title
    text(0,0, overall, cex=cex)
  }
  
  # revert
  layout(1)
  par(opa)
  return(invisible(NULL))
}
