##
##     texture.R
##
##     Texture plots and texture maps
##
##  $Revision: 1.8 $ $Date: 2015/04/18 18:11:52 $

### .................. basic graphics .............................

## put hatching in a window
add.texture <- function(W, texture=4, spacing=NULL, ...) {
  if(is.data.frame(texture)) {
    ## texture = f(x) where f is a texturemap
    out <- do.call("add.texture",
                   resolve.defaults(list(W=W, spacing=spacing),
                                    list(...),
                                    as.list(texture)))
    return(out)
  }
  ## texture should be an integer
  stopifnot(is.owin(W))
  stopifnot(texture %in% 1:8)
  if(is.null(spacing)) {
    spacing <- diameter(as.rectangle(W))/50
  } else {
    check.1.real(spacing)
    stopifnot(spacing > 0)
  }
  P <- L <- NULL
  switch(texture,
         {
           ## texture 1: graveyard
           P <- rsyst(W, dx=3*spacing)
         },
         {
           ## texture 2: vertical lines
           L <- rlinegrid(90, spacing, W)[W]
         },
         {
           ## texture 3: horizontal lines
           L <- rlinegrid(0, spacing, W)[W]
         },
         {
           ## texture 4: forward slashes
           L <- rlinegrid(45, spacing, W)[W]
         },
         {
           ## texture 5: back slashes
           L <- rlinegrid(135, spacing, W)[W]
         },
         {
           ## texture 6: horiz/vert grid
           L0 <- rlinegrid(0, spacing, W)[W]
           L90 <- rlinegrid(90, spacing, W)[W]
           L <- superimpose(L0, L90, W=W, check=FALSE)
         },
         {
           ## texture 7: diagonal grid
           L45 <- rlinegrid(45, spacing, W)[W]
           L135 <- rlinegrid(135, spacing, W)[W]
           L <- superimpose(L45, L135, W=W, check=FALSE)
         },
         {
           ## texture 8: hexagons
           H <- hextess(W, spacing, offset=runifpoint(1, W))
           H <- intersect.tess(H, W)
           do.call.matched("plot.tess",
                           resolve.defaults(list(x=H, add=TRUE),
                                            list(...)))
         })
  if(!is.null(P))
    do.call.matched("plot.ppp",
                    resolve.defaults(list(x=P, add=TRUE),
                                     list(...),
                                     list(chars=3, cex=0.2)),
                    extrargs=c("lwd", "col", "cols", "pch"))
  if(!is.null(L))
    do.call.matched("plot.psp",
                    resolve.defaults(list(x=L, add=TRUE),
                                     list(...)),
                    extrargs=c("lwd","lty","col"))
  return(invisible(NULL))
}

## .................. texture maps ................................

## create a texture map

texturemap <- function(inputs, textures, ...) {
  argh <- list(...)
  isnul <- unlist(lapply(argh, is.null))
  argh <- argh[!isnul]
  df <- do.call("data.frame",
                append(list(input=inputs, texture=textures), argh))
  f <- function(x) {
    df[match(x, df$input), -1, drop=FALSE]
  }
  class(f) <- c("texturemap", class(f))
  attr(f, "df") <- df
  return(f)
}

print.texturemap <- function(x, ...) {
  cat("Texture map\n")
  print(attr(x, "df"))
  return(invisible(NULL))
}

## plot a texture map

plot.texturemap <- local({

  ## recognised additional arguments to and axis()
  axisparams <- c("cex", 
                  "cex.axis", "cex.lab",
                  "col.axis", "col.lab",
                  "font.axis", "font.lab",
                  "las", "mgp", "xaxp", "yaxp",
                  "tck", "tcl", "xpd")

  # rules to determine the map dimensions when one dimension is given
  widthrule <- function(heightrange, separate, n, gap) {
    if(separate) 1 else diff(heightrange)/10
  }
  heightrule <- function(widthrange, separate, n, gap) {
    (if(separate) (n + (n-1)*gap) else 10) * diff(widthrange) 
  }

  plot.texturemap <- function(x, ..., main,
                             xlim=NULL, ylim=NULL, vertical=FALSE, axis=TRUE,
                             labelmap=NULL, gap=0.25,
                             spacing=NULL, add=FALSE) {
    if(missing(main))
      main <- short.deparse(substitute(x))
    df <- attr(x, "df")
#    textures <- df$textures
    n   <- nrow(df)
    check.1.real(gap, "In plot.texturemap")
    explain.ifnot(gap >= 0, "In plot.texturemap")
    separate <- (gap > 0)
    if(is.null(labelmap)) {
      labelmap <- function(x) x
    } else stopifnot(is.function(labelmap))
    ## determine rectangular window for display
    rr <- c(0, n + (n-1)*gap)
    if(is.null(xlim) && is.null(ylim)) {
      u <- widthrule(rr, separate, n, gap)
      if(!vertical) {
        xlim <- rr
        ylim <- c(0,u)
      } else {
        xlim <- c(0,u)
        ylim <- rr
      }
    } else if(is.null(ylim)) {
      if(!vertical) 
        ylim <- c(0, widthrule(xlim, separate, n, gap))
      else 
        ylim <- c(0, heightrule(xlim, separate, n, gap))
    } else if(is.null(xlim)) {
      if(!vertical) 
        xlim <- c(0, heightrule(ylim, separate, n, gap))
      else 
        xlim <- c(0, widthrule(ylim, separate, n, gap))
    } 
    width <- diff(xlim)
    height <- diff(ylim)
    ## determine boxes to be filled with textures,
    if(vertical) {
      boxheight <- min(width, height/(n + (n-1) * gap))
      vgap   <- (height - n * boxheight)/(n-1)
      boxes <- list()
      for(i in 1:n) boxes[[i]] <-
        owin(xlim, ylim[1] + c(i-1, i) * boxheight + (i-1) * vgap)
    } else {
      boxwidth <- min(height, width/(n + (n-1) * gap))
      hgap   <- (width - n * boxwidth)/(n-1)
      boxes <- list()
      for(i in 1:n) boxes[[i]] <-
        owin(xlim[1] + c(i-1, i) * boxwidth + (i-1) * hgap, ylim)
    }
    boxsize <- shortside(boxes[[1]])
    if(is.null(spacing))
      spacing <- 0.1 * boxsize
    
    # .......... initialise plot ...............................
    if(!add)
      do.call.matched("plot.default",
                      resolve.defaults(list(x=xlim, y=ylim,
                                            type="n", main=main,
                                            axes=FALSE, xlab="", ylab="",
                                            asp=1.0),
                                       list(...)))
    
    ## ................ plot texture blocks .................
    for(i in 1:n) {
      dfi <- df[i,,drop=FALSE]
      add.texture(W=boxes[[i]], texture=dfi, ..., spacing=spacing)
      plot(boxes[[i]], add=TRUE)
    }

    if(axis) {
      # ................. draw annotation ..................
      la <- paste(labelmap(df$input))
      if(!vertical) {
        ## add horizontal axis/annotation
        at <- lapply(lapply(boxes, centroid.owin), "getElement", name="x")
        # default axis position is below the ribbon (side=1)
        sidecode <- resolve.1.default("side", list(...), list(side=1))
        if(!(sidecode %in% c(1,3)))
          warning(paste("side =", sidecode,
                        "is not consistent with horizontal orientation"))
        pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
        # don't draw axis lines if plotting separate blocks
        lwd0 <- if(separate) 0 else 1
        # draw axis
        do.call.matched("axis",
                        resolve.defaults(list(...),
                                         list(side = 1, pos = pos, at = at),
                                         list(labels=la, lwd=lwd0)),
                        extrargs=axisparams)
      } else {
        ## add vertical axis
        at <- lapply(lapply(boxes, centroid.owin), "getElement", name="y")
        # default axis position is to the right of ribbon (side=4)
        sidecode <- resolve.1.default("side", list(...), list(side=4))
        if(!(sidecode %in% c(2,4)))
          warning(paste("side =", sidecode,
                        "is not consistent with vertical orientation"))
        pos <- c(ylim[1], xlim[1], ylim[2], xlim[2])[sidecode]
        # don't draw axis lines if plotting separate blocks
        lwd0 <- if(separate) 0 else 1
        # draw labels horizontally if plotting separate blocks
        las0 <- if(separate) 1 else 0
        # draw axis
        do.call.matched("axis",
                        resolve.defaults(list(...),
                                         list(side=4, pos=pos, at=at),
                                         list(labels=la, lwd=lwd0, las=las0)),
                        extrargs=axisparams)
      }
    }
    invisible(NULL)
  }

  plot.texturemap
})

## plot a pixel image using textures

textureplot <- local({

  textureplot <- function(x, ..., main, add=FALSE, clipwin=NULL, do.plot=TRUE,
                          border=NULL, col=NULL, lwd=NULL, lty=NULL,
                          spacing=NULL, textures=1:8,
                          legend=TRUE,
                          leg.side=c("right", "left", "bottom", "top"),
                          legsep=0.1, legwid=0.2) {
    if(missing(main))
      main <- short.deparse(substitute(x))
    if(!(is.im(x) || is.tess(x))) {
      x <- try(as.tess(x), silent=TRUE)
      if(inherits(x, "try-error"))
        stop("x must be a pixel image or a tessellation", call.=FALSE)
    }
    leg.side <- match.arg(leg.side)
    if(!is.null(clipwin))
      x <- x[clipwin, drop=FALSE]
    if(is.im(x)) {
      if(x$type != "factor")
        x <- eval.im(factor(x))
      levX <- levels(x)
    } else {
      tilX <- tiles(x)
      levX <- names(tilX)
    }
    n <- length(levX)
    if(n > 8)
      stop("Too many factor levels or tiles: maximum is 8")
    ## determine texture map
    if(inherits(textures, "texturemap")) {
      tmap <- textures
    } else {
      stopifnot(all(textures %in% 1:8))
      stopifnot(length(textures) >= n)
      mono <- spatstat.options("monochrome")
      col <- enforcelength(col, n, if(mono) 1 else 1:8)
      lwd <- if(is.null(lwd)) NULL else enforcelength(lwd, n, 1)
      lty <- if(is.null(lty)) NULL else enforcelength(lwd, n, 1)
      tmap <- texturemap(inputs=levX, textures=textures[1:n],
                         col=col, lwd=lwd, lty=lty)
    }
    ## determine plot region
    bb <- as.rectangle(x)
    if(!legend) {
      bb.all <- bb
    } else {
      Size <- max(sidelengths(bb))
      bb.leg <-
        switch(leg.side,
               right={
                 ## legend to right of image
                 owin(bb$xrange[2] + c(legsep, legsep+legwid) * Size,
                      bb$yrange)
               },
               left={
                 ## legend to left of image
                 owin(bb$xrange[1] - c(legsep+legwid, legsep) * Size,
                      bb$yrange)
               },
               top={
                 ## legend above image
                 owin(bb$xrange,
                      bb$yrange[2] + c(legsep, legsep+legwid) * Size)
               },
               bottom={
                 ## legend below image
                 owin(bb$xrange,
                      bb$yrange[1] - c(legsep+legwid, legsep) * Size)
           })
      iside <- match(leg.side, c("bottom", "left", "top", "right"))
      bb.all <- boundingbox(bb.leg, bb)
    }
    ## 
    result <- tmap
    attr(result, "bbox") <- bb
    ##
    if(do.plot) {
      ## Plot textures
      if(!add) {
        plot(bb.all, type="n", main="")
        fakemaintitle(bb, main, ...)
      }
      if(is.null(spacing)) spacing <- diameter(as.rectangle(x))/50
      areas <- if(is.im(x)) table(x$v) else tile.areas(x)
      for(i in which(areas > 0)) {
        Zi <- if(is.tess(x)) tilX[[i]] else levelset(x, levX[i], "==")
        Zi <- as.polygonal(Zi)
        if(is.null(border) || !is.na(border))
          plot(Zi, add=TRUE, border=border)
        add.texture(Zi, texture=tmap(levX[i]), spacing=spacing, ...)
      }
      vertical <- leg.side %in% c("left", "right")
      if(legend)
        do.call("plot.texturemap",
                resolve.defaults(list(x=tmap, add=TRUE,
                                      vertical=vertical,
                                      side=iside,
                                      xlim=bb.leg$xrange,
                                      ylim=bb.leg$yrange,
                                      spacing=spacing),
                                 list(...)))
    }
    return(invisible(result))
  }

  enforcelength <- function(x, n, x0) {
    if(is.null(x)) x <- x0
    if(length(x) < n) x <- rep(x, n)
    return(x[1:n])
  }

  textureplot
})



  
