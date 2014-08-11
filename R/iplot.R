#
# interactive plot for ppp objects using rpanel
#
#   $Revision: 1.14 $   $Date: 2014/03/22 07:03:23 $
#
#

# Effect:
# when the user types
#                 iplot(x)
# a pop-up panel displays a standard plot of x and
# buttons allowing control of the plot parameters.

# Coding:
# The panel 'p' contains the following internal variables
#      x          Original point pattern
#      w          Window of point pattern
#      xname      Name of x (for main title)
#      mtype      Type of marks of x
#      bb         frame of x 
#      bbmid      midpoint of frame
# The following variables in 'p' are controlled by panel buttons etc
#      split      Logical: whether to split multitype pattern
#      pointmap   Plot character, or "marks" indicating that marks are used
#      zoomfactor Zoom factor 
#      zoomcentre Centre point for zoom
#      charsize   Character expansion factor cex
#      markscale  Mark scale factor markscale
#      

iplot <- function(x, ...) {
  UseMethod("iplot")
}

iplot.ppp <- local({

iplot.ppp <- function(x, ..., xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  verifyclass(x, "ppp")
  require(rpanel)
  
  if(markformat(x) %in% c("hyperframe", "listof")) 
    marks(x) <- as.data.frame(as.hyperframe(marks(x)))
  if(markformat(x) == "dataframe" && ncol(marks(x)) > 1) {
    warning("Using only the first column of marks")
    marks(x) <- marks(x)[,1]
  }
  mtype <- if(is.multitype(x)) "multitype" else if(is.marked(x)) "marked" else "unmarked"

  bb <- as.rectangle(as.owin(x))
  bbmid <- unlist(centroid.owin(bb))
  ##
  p <- rp.control(paste("iplot(", xname, ")", sep=""), 
                  x=x,
                  w=as.owin(x),
                  xname=xname,
                  mtype=mtype,
                  bb=bb,
                  bbmid=bbmid,
                  split=FALSE,
                  pointmap=if(is.marked(x)) "marks" else "o",
                  zoomfactor=1,
                  zoomcentre=bbmid,
                  size=c(700, 400))

# Split panel into three
# Left: plot controls
# Middle: data
# Right: navigation/zoom
  rp.grid(p, "gcontrols", pos=list(row=0,column=0))
  rp.grid(p, "gdisplay",  pos=list(row=0,column=1))
  rp.grid(p, "gnavigate", pos=list(row=0,column=2))

#----- Data display ------------

  # This line is to placate the package checker
  mytkr <- NULL

  # Create data display panel 
  rp.tkrplot(p, mytkr, plotfun=do.iplot.ppp, action=click.iplot.ppp,
             pos=list(row=0,column=0,grid="gdisplay"))

  
#----- Plot controls ------------
  nextrow <- 0
  pozzie <- function(n=nextrow, ...)
    append(list(row=n,column=0,grid="gcontrols"), list(...))
  
# main title
  rp.textentry(p, xname, action=redraw.iplot.ppp, title="Plot title",
               pos=pozzie(0))
  nextrow <- 1

# split ?
  if(mtype == "multitype") {
    rp.checkbox(p, split, initval=FALSE, 
                title="Split according to marks", action=redraw.iplot.ppp,
                pos=pozzie(1))
    nextrow <- 2
  }

# plot character or mark style
  ptvalues <- c("o", "bullet", "plus")
  ptlabels <- c("open circles", "filled circles", "crosshairs")
  if(is.marked(x)) {
    ptvalues <- c("marks", ptvalues)
    ptlabels <- if(mtype == "multitype")
      c("Symbols depending on mark", ptlabels)
    else c("Circles proportional to mark", ptlabels)
  }
  pointmap <- ptvalues[1]
  rp.radiogroup(p, pointmap, vals=ptvalues, labels=ptlabels,
   			  title="how to plot points", action=redraw.iplot.ppp,
                pos=pozzie(nextrow))
  nextrow <- nextrow+1

# plot character size
  charsize <- 1
  rp.slider(p, charsize, 0, 5, action=redraw.iplot.ppp, 
            title="symbol expansion factor (cex)", initval=1, showvalue=TRUE,
            pos=pozzie(nextrow, sticky=""))
  nextrow <- nextrow+1
  
# mark scale
  if(mtype == "marked") {
    marx <- x$marks
    marx <- marx[is.finite(marx)]
    scal <- mark.scale.default(marx, x$window)
    markscale <- scal
    rp.slider(p, markscale, from=scal/10, to = 10*scal, action=redraw.iplot.ppp,
              initval=scal,
              title="mark scale factor (markscale)", showvalue=TRUE,
              pos=pozzie(nextrow))
    nextrow <- nextrow+1
  }

# button to print a summary at console
  rp.button(p, title="Print summary information",
            pos=pozzie(nextrow),
            action=function(panel) { print(summary(panel$x)); panel} )
#  
#----- Navigation controls ------------
  nextrow <- 0
  navpos <- function(n=nextrow,cc=0, ...)
    append(list(row=n,column=cc,grid="gnavigate"), list(...))

  rp.button(p, title="Up", pos=navpos(nextrow,1,sticky=""),
            action=function(panel) {
              zo <- panel$zoomfactor
              ce <- panel$zoomcentre
              bb <- panel$bb
              height <- sidelengths(bb)[2]
              stepsize <- (height/4)/zo
              panel$zoomcentre <- ce + c(0, stepsize)
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1
  rp.button(p, title="Left", pos=navpos(nextrow,0,sticky="w"),
            action=function(panel) {
              zo <- panel$zoomfactor
              ce <- panel$zoomcentre
              bb <- panel$bb
              width <- sidelengths(bb)[1]
              stepsize <- (width/4)/zo
              panel$zoomcentre <- ce - c(stepsize, 0)
              redraw.iplot.ppp(panel)
              return(panel)
            })
  rp.button(p, title="Right", pos=navpos(nextrow,2,sticky="e"),
            action=function(panel) {
              zo <- panel$zoomfactor
              ce <- panel$zoomcentre
              bb <- panel$bb
              width <- sidelengths(bb)[1]
              stepsize <- (width/4)/zo
              panel$zoomcentre <- ce + c(stepsize, 0)
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1
  rp.button(p, title="Down", pos=navpos(nextrow,1,sticky=""),
            action=function(panel) {
              zo <- panel$zoomfactor
              ce <- panel$zoomcentre
              bb <- panel$bb
              height <- sidelengths(bb)[2]
              stepsize <- (height/4)/zo
              panel$zoomcentre <- ce - c(0, stepsize)
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1

  rp.button(p, title="Zoom In", pos=navpos(nextrow,1,sticky=""),
            action=function(panel) {
              panel$zoomfactor <- panel$zoomfactor * 2
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1
  rp.button(p, title="Zoom Out", pos=navpos(nextrow,1,sticky=""),
            action=function(panel) {
              panel$zoomfactor <- panel$zoomfactor / 2
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1
  rp.button(p, title="Reset", pos=navpos(nextrow,1,sticky=""),
            action=function(panel) {
              panel$zoomfactor <- 1
              panel$zoomcentre <- panel$bbmid
              redraw.iplot.ppp(panel)
              return(panel)
            })
  nextrow <- nextrow + 1
  rp.button(p, title="Redraw", pos=navpos(nextrow,1,sticky=""),
            action=redraw.iplot.ppp)
  nextrow <- nextrow+1
# quit button 
  rp.button(p, title="Quit", quitbutton=TRUE,
            pos=navpos(nextrow, 1, sticky=""),
            action= function(panel) { panel })

  invisible(NULL)
}


  # Function to redraw the whole shebang
  redraw.iplot.ppp <- function(panel) {
    rp.tkrreplot(panel, mytkr)
    panel
  }


# Function executed when data display is clicked

  click.iplot.ppp <- function(panel, x, y) {
    if(panel$split) {
      cat("Mouse interaction is not supported when the point pattern is split\n")
    } else {
      panel$zoomcentre <- panel$zoomcentre +
        (c(x,y) - panel$bbmid)/panel$zoomfactor
      redraw.iplot.ppp(panel)
    }
    return(panel)
  }

# function that updates the plot when the control panel is operated

do.iplot.ppp <- function(panel) { 
  use.marks <- TRUE
  pch <- 16
  switch(panel$pointmap,
         marks={
           use.marks <- TRUE
           pch <- NULL
         }, 
         o = {
           use.marks <- FALSE
           pch <- 1
         }, 
         bullet = {
           use.marks <- FALSE
           pch <- 16
         },
         plus = {
           use.marks <- FALSE
           pch <- 3
         })
  # scale and clip the pattern
  x <- panel$x
  w     <- panel$w
  z     <- panel$zoomfactor
  if(is.null(z)) z <- 1
  ce    <- panel$zoomcentre
  bb    <- panel$bb
  bbmid <- panel$bbmid
  scalex <- shift(affine(shift(x, -ce), diag(c(z,z))), bbmid)
  scalew <- shift(affine(shift(w, -ce), diag(c(z,z))), bbmid)
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
  if(panel$mtype == "multitype" && panel$split) {
    scalex <- split(scalex, un=(panel$pointmap != "marks"))
    plot(scalex, main=panel$xname, 
         use.marks=use.marks, pch=pch, cex=panel$charsize,
         panel.begin=panel.begin)
  } else {
    # draw scaled & clipped window
    plot(panel.begin, main=panel$xname)
    # add points
    if(panel$mtype == "marked" && panel$pointmap == "marks") {
      plot(scalex, add=TRUE, use.marks=use.marks, markscale=panel$markscale)
    } else {
      plot(scalex, add=TRUE, use.marks=use.marks, pch=pch, cex=panel$charsize)
    }
  }
  par(opa)
  panel
}

iplot.ppp
})
