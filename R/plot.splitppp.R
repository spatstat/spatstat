#
plot.listof <- plot.splitppp <-
  local({

  # auxiliary functions
  extraplot <- function(nnn, ..., panel.args=NULL, plotcommand="plot") {
    if(is.null(panel.args)) {
      do.call(plotcommand, list(...))
    } else {
      xtra <- if(is.function(panel.args)) panel.args(nnn) else panel.args
      if(!is.list(xtra)) stop("panel.args should be a list")
      do.call(plotcommand, append(list(...), xtra))
    }
  }

  exec.or.plot <- function(cmd, i, xi, ...) {
    if(is.null(cmd)) return(NULL)
    if(is.function(cmd)) {
      do.call(cmd, resolve.defaults(list(i, xi, ...)))
    } else {
      do.call(plot, resolve.defaults(list(cmd, ...)))
    }
  }
  
  plot.splitppp <- function(x, ..., main, arrange=TRUE,
                            nrows=NULL, ncols=NULL,
                            main.panel=NULL,
                            mar.panel=c(2,1,1,2),
                            panel.begin=NULL,
                            panel.end=NULL,
                            panel.args=NULL,
                            plotcommand="plot",
                            adorn.left=NULL,
                            adorn.right=NULL,
                            adorn.top=NULL,
                            adorn.bottom=NULL,
                            adorn.size=0.2) {
    xname <- short.deparse(substitute(x))

    # `boomerang despatch'
    cl <- match.call()
    if(missing(plotcommand) && all(unlist(lapply(x, is.im)))) {
      cl[[1]] <- as.name("image.listof")
      parenv <- sys.parent()
      return(eval(cl, envir=parenv))
    }
            
    n <- length(x)
    names(x) <- good.names(names(x), "Component_", 1:n)
    if(is.null(main.panel))
      main.panel <- names(x)
    else {
      stopifnot(is.character(main.panel) || is.expression(main.panel))
      nmp <- length(main.panel)
      if(nmp == 1)
        main.panel <- rep.int(main.panel, n)
      else if(nmp != n)
        stop("Incorrect length for main.panel")
    }

    if(!arrange) {
      # sequence of plots
      for(i in 1:n) {
        xi <- x[[i]]
        exec.or.plot(panel.begin, i, xi, main=main.panel[i])
        extraplot(i, xi, ...,
                  add=!is.null(panel.begin),
                  main=main.panel[i],
                  panel.args=panel.args, plotcommand=plotcommand)
        exec.or.plot(panel.end, i, xi, add=TRUE)
      }
      if(!is.null(adorn.left))
        warning("adorn.left was ignored because arrange=FALSE")
      if(!is.null(adorn.right))
        warning("adorn.right was ignored because arrange=FALSE")
      if(!is.null(adorn.top))
        warning("adorn.top was ignored because arrange=FALSE")
      if(!is.null(adorn.bottom))
        warning("adorn.bottom was ignored because arrange=FALSE")
      return(invisible(NULL))
    }

    # ARRAY of plots
    # decide whether to plot a main header
    main <- if(!missing(main) && !is.null(main)) main else xname
    if(!is.character(main)) {
      # main title could be an expression
      nlines <- 1
      banner <- TRUE
    } else {
      # main title is character string/vector, possibly ""
      banner <- any(nzchar(main))
      if(length(main) > 1)
        main <- paste(main, collapse="\n")
      nlines <- length(unlist(strsplit(main, "\n")))
    }
    # determine arrangement of plots
    # arrange like mfrow(nrows, ncols) plus a banner at the top
    if(is.null(nrows) && is.null(ncols)) {
      nrows <- as.integer(floor(sqrt(n)))
      ncols <- as.integer(ceiling(n/nrows))
    } else if(!is.null(nrows) && is.null(ncols))
      ncols <- as.integer(ceiling(n/nrows))
    else if(is.null(nrows) && !is.null(ncols))
      nrows <- as.integer(ceiling(n/ncols))
    else stopifnot(nrows * ncols >= length(x))
    nblank <- ncols * nrows - n
    # determine approximate relative dimensions for equal scale plots
    boxes <- try(lapply(x, as.rectangle), silent=TRUE)
    sizes.known <- !inherits(boxes, "try-error")
    # set up layout
    mat <- matrix(c(seq_len(n), integer(nblank)),
                  byrow=TRUE, ncol=ncols, nrow=nrows)
    if(sizes.known) {
      xwidths <- unlist(lapply(boxes, function(z) { diff(z$xrange) }))
      xheights <- unlist(lapply(boxes, function(z) { diff(z$yrange) }))
      heights <- apply(mat, 1, function(j,h) { max(h[j[j>0]]) }, h=xheights)
      widths <- apply(mat, 2, function(i,w) { max(w[i[i>0]]) }, w=xwidths)
    } else {
      heights <- rep.int(1, nrows)
      widths <- rep.int(1, ncols)
    }
    meanheight <- mean(heights)
    meanwidth  <- mean(widths)
    nall <- n
    if(!is.null(adorn.left)) {
      # add margin at left, of width adorn.size * meanwidth
      nall <- i.left <- n+1
      mat <- cbind(i.left, mat)
      widths <- c(adorn.size * meanwidth, widths)
    } 
    if(!is.null(adorn.right)) {
      # add margin at right, of width adorn.size * meanwidth
      nall <- i.right <- nall+1
      mat <- cbind(mat, i.right)
      widths <- c(widths, adorn.size * meanwidth)
    } 
    if(!is.null(adorn.bottom)) {
      # add margin at bottom, of height adorn.size * meanheight
      nall <- i.bottom <- nall+1
      mat <- rbind(mat, i.bottom)
      heights <- c(heights, adorn.size * meanheight)
    } 
    if(!is.null(adorn.top)) {
      # add margin at top, of height adorn.size * meanheight
      nall <- i.top <- nall + 1
      mat <- rbind(i.top, mat)
      heights <- c(adorn.size * meanheight, heights)
    } 
    if(banner) {
      # Increment existing panel numbers
      # New panel 1 is the banner
      panels <- (mat > 0)
      mat[panels] <- mat[panels] + 1
      mat <- rbind(1, mat)
      heights <- c(0.1 * meanheight * (1 + nlines), heights)
    }
    # declare layout
    layout(mat, heights=heights, widths=widths, respect=sizes.known)
    # start output .....
    # .... plot banner
    if(banner) {
      opa <- par(mar=rep.int(0,4), xpd=TRUE)
      plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
           xlim=c(-1,1),ylim=c(-1,1))
      cex <- resolve.defaults(list(...), list(cex.title=2))$cex.title
      text(0,0,main, cex=cex)
    }
    # plot panels
    npa <- par(mar=mar.panel)
    if(!banner) opa <- npa
    for(i in 1:n) {
      xi <- x[[i]]
      exec.or.plot(panel.begin, i, xi, main=main.panel[i])
      extraplot(i, xi, ...,
                add = !is.null(panel.begin), 
                main = main.panel[i],
                panel.args=panel.args, plotcommand=plotcommand)
      exec.or.plot(panel.end, i, xi, add=TRUE)
    }
    # adornments
    if(nall > n) {
      par(mar=rep.int(0,4), xpd=TRUE)
      if(!is.null(adorn.left))
        adorn.left()
      if(!is.null(adorn.right))
        adorn.right()
      if(!is.null(adorn.bottom))
        adorn.bottom()
      if(!is.null(adorn.top))
        adorn.top()
    }
    # revert
    layout(1)
    par(opa)
    return(invisible(NULL))
  }

  plot.splitppp
})

density.splitppp <- function(x, ...) {
  as.listof(lapply(x, density, ...))
}
