#
# interactive plot 
#
#   $Revision: 1.11 $   $Date: 2015/05/07 01:56:47 $
#
#

iplot.default <- function(x, ..., xname) {
 if(missing(xname))
    xname <- short.deparse(substitute(x))
 x <- as.layered(x)
 iplot(x, ..., xname=xname)
}

iplot.layered <- local({

  CommitAndRedraw <- function(panel) {
    ## hack to ensure that panel is immediately updated in rpanel
    kraever("rpanel")
    ## This is really a triple-colon!
    rpanel:::rp.control.put(panel$panelname, panel)
    ## now redraw it
    redraw.iplot.layered(panel)
  }
  
  faster.layers <- function(x, visible) {
    if(any(islinnet <- unlist(lapply(x, inherits, what="linnet")))) {
      # convert linnet layers to psp, for efficiency
      x[islinnet] <- lapply(x[islinnet], as.psp)
    }
    repeat{
      islpp <- unlist(lapply(x, inherits, what="lpp"))
      if(!any(islpp))
        break
      # convert an lpp layer to two layers: psp and ppp, for efficiency
      ii <- min(which(islpp))
      pl <- layerplotargs(x)
      n <- length(x)
      xpre <- if(ii == 1) NULL else x[1:ii]
      xpost <- if(ii == n) NULL else x[(ii+1):n]
      ppre <- if(ii == 1) NULL else pl[1:ii]
      ppost <- if(ii == n) NULL else pl[(ii+1):n]
      a <- as.psp(as.linnet(x[[ii]]))
      b <- as.ppp(x[[ii]])
      x <- layered(LayerList=c(xpre, list(a, b), xpost),
                   plotargs=unname(c(ppre, pl[ii], pl[ii], ppost)))
      visible <- visible[if(ii == 1) c(1, seq_len(n)) else
                         if(ii == n) c(seq_len(n), n) else
                         c(1:(ii-1), ii, ii, (ii+1):n)]
    }
    attr(x, "visible") <- visible
    return(x)
  }
  
iplot.layered <- function(x, ..., xname, visible) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  verifyclass(x, "layered")

  if(missing(visible) || is.null(visible)) {
    visible <- rep(TRUE, length(x))
  } else if(length(visible) == 1) {
    visible <- rep(visible, length(x))
  } else stopifnot(length(visible) == length(x))

  kraever("rpanel")

  x <- faster.layers(x, visible)
  visible <- attr(x, "visible")

  x <- freeze.colourmaps(x)
  bb <- as.rectangle(as.owin(x))
  bbmid <- unlist(centroid.owin(bb))


  lnames <- names(x)
  if(sum(nzchar(lnames)) != length(x))
    lnames <- paste("Layer", seq_len(length(x)))
  ##
  p <- rpanel::rp.control(paste("iplot(", xname, ")", sep=""), 
                          x=x,
                          w=as.owin(x),
                          xname=xname,
                          layernames=lnames,
                          bb=bb,
                          bbmid=bbmid,
                          zoomfactor=1,
                          zoomcentre=bbmid,
                          which = visible,
                          size=c(700, 400))

# Split panel into three
# Left: plot controls
# Middle: data
# Right: navigation/zoom
  rpanel::rp.grid(p, "gcontrols", pos=list(row=0,column=0))
  rpanel::rp.grid(p, "gdisplay",  pos=list(row=0,column=1))
  rpanel::rp.grid(p, "gnavigate", pos=list(row=0,column=2))

#----- Data display ------------

  # This line is to placate the package checker
  mytkr <- NULL

  # Create data display panel 
  rpanel::rp.tkrplot(p, mytkr, plotfun=do.iplot.layered,
                     action=click.iplot.layered,
                     pos=list(row=0,column=0,grid="gdisplay"))

  
#----- Plot controls ------------
  nextrow <- 0
  pozzie <- function(n=nextrow, ...)
    append(list(row=n,column=0,grid="gcontrols"), list(...))
  
# main title
  rpanel::rp.textentry(p, xname, action=redraw.iplot.layered,
                       title="Plot title",
                       pos=pozzie(0))
  nextrow <- 1

# select some layers
  nx <- length(x)
  which <- rep.int(TRUE, nx)
  if(nx > 1) {
    rpanel::rp.checkbox(p, which, labels=lnames,
                        action=redraw.iplot.layered,
                        title="Select layers to plot",
                        pos=pozzie(nextrow), sticky="")
    nextrow <- nextrow + 1
  }
  
# button to print a summary at console
  rpanel::rp.button(p, title="Print summary information",
                    pos=pozzie(nextrow),
                    action=function(panel) {
                        lapply(panel$x, function(z) print(summary(z)))
                        return(panel)
                    })
#  
#----- Navigation controls ------------
  nextrow <- 0
  navpos <- function(n=nextrow,cc=0, ...)
    append(list(row=n,column=cc,grid="gnavigate"), list(...))

  rpanel::rp.button(p, title="Up", pos=navpos(nextrow,1,sticky=""),
                    action=function(panel) {
                        zo <- panel$zoomfactor
                        ce <- panel$zoomcentre
                        bb <- panel$bb
                        height <- sidelengths(bb)[2]
                        stepsize <- (height/4)/zo
                        panel$zoomcentre <- ce + c(0, stepsize)
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1
  rpanel::rp.button(p, title="Left", pos=navpos(nextrow,0,sticky="w"),
                    action=function(panel) {
                        zo <- panel$zoomfactor
                        ce <- panel$zoomcentre
                        bb <- panel$bb
                        width <- sidelengths(bb)[1]
                        stepsize <- (width/4)/zo
                        panel$zoomcentre <- ce - c(stepsize, 0)
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  rpanel::rp.button(p, title="Right", pos=navpos(nextrow,2,sticky="e"),
                    action=function(panel) {
                        zo <- panel$zoomfactor
                        ce <- panel$zoomcentre
                        bb <- panel$bb
                        width <- sidelengths(bb)[1]
                        stepsize <- (width/4)/zo
                        panel$zoomcentre <- ce + c(stepsize, 0)
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1
  rpanel::rp.button(p, title="Down", pos=navpos(nextrow,1,sticky=""),
                    action=function(panel) {
                        zo <- panel$zoomfactor
                        ce <- panel$zoomcentre
                        bb <- panel$bb
                        height <- sidelengths(bb)[2]
                        stepsize <- (height/4)/zo
                        panel$zoomcentre <- ce - c(0, stepsize)
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1

  rpanel::rp.button(p, title="Zoom In", pos=navpos(nextrow,1,sticky=""),
                    action=function(panel) {
                        panel$zoomfactor <- panel$zoomfactor * 2
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1
  rpanel::rp.button(p, title="Zoom Out", pos=navpos(nextrow,1,sticky=""),
                    action=function(panel) {
                        panel$zoomfactor <- panel$zoomfactor / 2
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1
  rpanel::rp.button(p, title="Reset", pos=navpos(nextrow,1,sticky=""),
                    action=function(panel) {
                        panel$zoomfactor <- 1
                        panel$zoomcentre <- panel$bbmid
                        CommitAndRedraw(panel)
                        return(panel)
                    })
  nextrow <- nextrow + 1
  rpanel::rp.button(p, title="Redraw", pos=navpos(nextrow,1,sticky=""),
                    action=redraw.iplot.layered)
  nextrow <- nextrow+1
# quit button 
  rpanel::rp.button(p, title="Quit", quitbutton=TRUE,
            pos=navpos(nextrow, 1, sticky=""),
            action= function(panel) { panel })

  invisible(NULL)
}


  # Function to redraw the whole shebang
  redraw.iplot.layered <- function(panel) {
    rpanel::rp.tkrreplot(panel, mytkr)
    panel
  }


# Function executed when data display is clicked

  click.iplot.layered <- function(panel, x, y) {
    panel$zoomcentre <- panel$zoomcentre +
      (c(x,y) - panel$bbmid)/panel$zoomfactor
    CommitAndRedraw(panel)
    return(panel)
  }

# function that updates the plot when the control panel is operated

do.iplot.layered <- function(panel) { 
  # scale and clip the pattern
  x <- panel$x[panel$which]
  w     <- panel$w
  z     <- panel$zoomfactor
  if(is.null(z)) z <- 1
  ce    <- panel$zoomcentre
  bb    <- panel$bb
  bbmid <- panel$bbmid
  scalex <- shift(scalardilate(shift(x, -ce), z), bbmid)
  scalew <- shift(scalardilate(shift(w, -ce), z), bbmid)
  scalex <- scalex[, bb]
  scalew <- intersect.owin(scalew, bb, fatal=FALSE)
  # determine what is plotted under the clipped pattern
  blankargs <- list(type="n")
  dashargs  <- list(lty=3, border="red")
  panel.begin <- 
    if(is.null(scalew)) {
      # empty intersection; just create the plot space
      layered(bb,          plotargs=list(blankargs))
    } else if(identical(bb, scalew)) {
      if(z == 1) {
        # original state
        # window is rectangular 
        # plot the data window as a solid black rectangle
        layered(bb, scalew,  plotargs=list(blankargs, list(lwd=2)))
      } else {
        # zoom view is entirely inside window
        # plot the clipping region as a red dashed rectangle
        layered(bb, plotargs=list(dashargs))
      }
    } else {
      # field of view is not a subset of window
      # plot the clipping region as a red dashed rectangle
      # Then add the data window
      layered(bb, scalew, plotargs=list(dashargs, list(invert=TRUE)))
    }
  
  # draw it
  opa <- par(ask=FALSE)
  plot(panel.begin, main=panel$xname)
  plot(scalex, add=TRUE)
  par(opa)
  panel
}

freeze.colourmaps <- function(x) {
  # tweak a layered object to ensure that
  # the colours of image layers don't change with zoom/pan
  isim <- unlist(lapply(x, is.im))
  if(any(isim)) {
    # ensure there are plotargs
    pl <- attr(x, "plotargs")
    if(is.null(pl))
      pl <- rep.int(list(list()), length(x))
    # make sure the plotargs include 'zlim'
    for(i in which(isim)) {
      x.i <- x[[i]]
      if(x.i$type %in% c("integer", "real")) 
        pl[[i]] <- resolve.defaults(pl[[i]], list(zlim=range(x.i)))
    }
    # put back
    attr(x, "plotargs") <- pl
  }
  return(x) 
}

iplot.layered
})
