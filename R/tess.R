#
# tess.R
#
# support for tessellations
#
#   $Revision: 1.71 $ $Date: 2015/10/21 09:06:57 $
#
tess <- function(..., xgrid=NULL, ygrid=NULL, tiles=NULL, image=NULL,
                 window=NULL, marks=NULL, keepempty=FALSE,
                 unitname=NULL, check=TRUE) {
  uname <- unitname
  if(!is.null(window)) {
    window <- as.owin(window)
    if(is.null(uname)) uname <- unitname(window) 
  }
  isrect <- !is.null(xgrid) && !is.null(ygrid)
  istiled <- !is.null(tiles)
  isimage <- !is.null(image)
  if(isrect + istiled + isimage != 1)
    stop("Must specify either (xgrid, ygrid) or tiles or img")
  if(isrect) {
    stopifnot(is.numeric(xgrid) && all(diff(xgrid) > 0))
    stopifnot(is.numeric(ygrid) && all(diff(ygrid) > 0))
    if(!is.null(window))
      warning("Argument 'window' ignored, because xgrid, grid are given")
    window <- owin(range(xgrid), range(ygrid), unitname=uname)
    ntiles <- (length(xgrid)-1) * (length(ygrid)-1)
    out <- list(type="rect", window=window, xgrid=xgrid, ygrid=ygrid, n=ntiles)
  } else if(istiled) {
    stopifnot(is.list(tiles))
    if(check) {
      if(!all(sapply(tiles, is.owin)))
        stop("Tiles must be a list of owin objects")
      if(!is.null(uname)) {
        ## attach new unit name to each tile
        tiles <- lapply(tiles, "unitname<-", value=uname)
      } else {
        ## extract unit names from tiles, check agreement, use as unitname
        uu <- unique(lapply(tiles, unitname))
        uu <- uu[!sapply(uu, is.null)]
        nun <- length(uu)
        if(nun > 1)
          stop("Tiles have inconsistent names for the unit of length")
        if(nun == 1) {
          ## use this unit name
          uname <- uu[[1]]
          if(!is.null(window))
            unitname(window) <- uname
        }
      }
    }
    if(!keepempty && check) {
      # remove empty tiles
      isempty <- sapply(tiles, is.empty)
      if(all(isempty))
        stop("All tiles are empty")
      if(any(isempty))
        tiles <- tiles[!isempty]
    }
    ntiles <- length(tiles)
    nam <- names(tiles)
    lev <- if(!is.null(nam) && all(nzchar(nam))) nam else 1:ntiles
    if(is.null(window)) 
      window <- do.call(union.owin, unname(tiles))
    if(is.mask(window) || any(unlist(lapply(tiles, is.mask)))) {
      # convert to pixel image tessellation
      window <- as.mask(window)
      ima <- as.im(window)
      ima$v[] <- NA
      for(i in 1:ntiles)
        ima[tiles[[i]]] <- i
      ima <- ima[window, drop=FALSE]
      ima <- eval.im(factor(ima, levels=1:ntiles))
      levels(ima) <- lev
      out <- list(type="image",
                  window=window, image=ima, n=length(lev))
    } else {
      # tile list
      window <- rescue.rectangle(window)
      out <- list(type="tiled", window=window, tiles=tiles, n=length(tiles))
    }
  } else if(isimage) {
    # convert to factor valued image
    image <- as.im(image)
    if(!is.null(uname)) unitname(image) <- uname
    switch(image$type,
           logical={
             # convert to factor
             if(keepempty) 
               image <- eval.im(factor(image, levels=c(FALSE,TRUE)))
             else
               image <- eval.im(factor(image))
           },
           factor={
             # eradicate unused levels
             if(!keepempty) 
               image <- eval.im(factor(image))
           },
           {
             # convert to factor
             image <- eval.im(factor(image))
           })
               
    if(is.null(window)) window <- as.owin(image)
    out <- list(type="image", window=window, image=image, n=length(levels(image)))
  } else stop("Internal error: unrecognised format")
  ## add marks!
  if(!is.null(marks)) {
    marks <- as.data.frame(marks)
    if(nrow(marks) != out$n)
      stop(paste("wrong number of marks:",
                 nrow(marks), "should be", out$n),
           call.=FALSE)
    out$marks <- marks
  }
  class(out) <- c("tess", class(out))
  return(out)
}

is.tess <- function(x) { inherits(x, "tess") }

print.tess <- function(x, ..., brief=FALSE) {
  full <- !brief
  if(full) cat("Tessellation\n")
  win <- x$window
  switch(x$type,
         rect={
           if(full) {
             unitinfo <- summary(unitname(win))
             if(equispaced(x$xgrid) && equispaced(x$ygrid)) 
               splat("Tiles are equal rectangles, of dimension",
                     signif(mean(diff(x$xgrid)), 5),
                     "x",
                     signif(mean(diff(x$ygrid)), 5),
                     unitinfo$plural, " ", unitinfo$explain)
             else
               splat("Tiles are unequal rectangles")
           }
           splat(length(x$xgrid)-1, "by", length(x$ygrid)-1, "grid of tiles")
         },
         tiled={
           if(full) {
             if(win$type == "polygonal")
               splat("Tiles are irregular polygons")
             else
               splat("Tiles are windows of general type")
           }
           splat(length(x$tiles), "tiles (irregular windows)")
         },
         image={
           nlev <- length(levels(x$image))
           if(full) {
             splat("Tessellation is determined by a factor-valued image with",
                   nlev, "levels")
           } else splat(nlev, "tiles (levels of a pixel image)")
         })
  if(!is.null(marx <- x$marks)) {
    m <- dim(marx)[2] %orifnull% 1
    if(m == 1) splat("Tessellation is marked") else
    splat("Tessellation has", m, "columns of marks:",
          commasep(sQuote(colnames(marx))))
  }
  if(full) print(win)
  invisible(NULL)
}

unitname.tess <- function(x) unitname(x$window)

"unitname<-.tess" <- function(x, value) {
  unitname(x$window) <- value
  switch(x$type,
         rect={},
         tiled={
           x$tiles <- lapply(x$tiles, "unitname<-", value)
         },
         image={
           unitname(x$image) <- value
         })
  return(x)
}

plot.tess <- local({

  plotem <- function(z, ..., col=NULL) {
    if(is.null(col))
      plot(z, ..., add=TRUE)
    else if(z$type != "mask")
      plot(z, ..., border=col, add=TRUE)
    else plot(z, ..., col=col, add=TRUE)
  }

  plotpars <- c("sub", "lty", "lwd",
                "cex.main", "col.main", "font.main",
                "cex.sub", "col.sub", "font.sub", "border")

  plot.tess <- function(x, ..., main, add=FALSE, show.all=!add, col=NULL,
                        do.plot=TRUE,
                        do.labels=FALSE, labels=tilenames(x),
                        labelargs=list()) {
    if(missing(main) || is.null(main))
      main <- short.deparse(substitute(x))
    switch(x$type,
           rect={
             win <- x$window
             result <-
               do.call.matched("plot.owin",
                               resolve.defaults(list(x=win, main=main,
                                                     add=add,
                                                     show.all=show.all,
                                                     do.plot=do.plot),
                                                list(...)),
                               extrargs=plotpars)
             if(do.plot) {
               xg <- x$xgrid
               yg <- x$ygrid
               do.call.matched("segments",
                               resolve.defaults(list(x0=xg, y0=win$yrange[1],
                                                     x1=xg, y1=win$yrange[2]),
                                                list(col=col),
                                                list(...),
                                                .StripNull=TRUE))
               do.call.matched("segments",
                               resolve.defaults(list(x0=win$xrange[1], y0=yg,
                                                     x1=win$xrange[2], y1=yg),
                                                list(col=col),
                                                list(...),
                                                .StripNull=TRUE))
             }
           },
           tiled={
             result <-
               do.call.matched("plot.owin",
                               resolve.defaults(list(x=x$window, main=main,
                                                     add=add,
                                                     show.all=show.all,
                                                     do.plot=do.plot),
                                                list(...)),
                               extrargs=plotpars)
             if(do.plot) {
               til <- tiles(x)
               lapply(til, plotem, ..., col=col)
             }
           },
           image={
             result <-
               do.call("plot",
                       resolve.defaults(list(x$image, add=add, main=main,
                                             show.all=show.all,
                                             do.plot=do.plot),
                                        list(...),
                                        list(valuesAreColours=FALSE)))
           })
    if(do.plot && do.labels) {
      labels <- paste(as.vector(labels))
      til <- tiles(x)
      incircles <- lapply(til, incircle)
      x0 <- sapply(incircles, getElement, name="x")
      y0 <- sapply(incircles, getElement, name="y")
      do.call.matched("text.default",
                      resolve.defaults(list(x=x0, y = y0),
                                       list(labels=labels),
                                       labelargs))
    }
    return(invisible(result))
  }

  plot.tess
})


"[<-.tess" <- function(x, i, ..., value) {
  switch(x$type,
         rect=,
         tiled={
           til <- tiles(x)
           til[i] <- value
           ok <- !unlist(lapply(til, is.null))
           x <- tess(tiles=til[ok])
         },
         image={
           stop("Cannot assign new values to subsets of a pixel image")
         })
  return(x)
}
  
"[.tess" <- function(x, i, ...) {
  trap.extra.arguments(..., .Context="in [.tess")
  if(missing(i)) return(x)
  if(is.owin(i))
    return(intersect.tess(x, i))
  switch(x$type,
         rect=,
         tiled={
           til <- tiles(x)[i]
           return(tess(tiles=til))
         },
         image={
           img <- x$image
           oldlev <- levels(img)
           newlev <- unique(oldlev[i])
           img <- eval.im(factor(img, levels=newlev))
           return(tess(image=img))
         })
}

tiles <- function(x) {
  switch(x$type,
         rect={
           out <- list()
           xg <- x$xgrid
           yg <- x$ygrid
           nx <- length(xg) - 1
           ny <- length(yg) - 1
           for(j in rev(seq_len(ny)))
             for(i in seq_len(nx)) {
               winij <- owin(xg[c(i,i+1)], yg[c(j,j+1)])
               dout <- list(winij)
               names(dout) <- paste("Tile row ", ny-j+1, ", col ", i,
                                    sep="")
               out <- append(out, dout)
             }
         },
         tiled={
           out <- x$tiles
           if(is.null(names(out)))
             names(out) <- paste("Tile", seq_along(out))
         },
         image={
           out <- list()
           ima <- x$image
           lev <- levels(ima)
           for(i in seq_along(lev))
             out[[i]] <- solutionset(ima == lev[i])
           names(out) <- paste(lev)
         })
  out <- as.solist(out)
  return(out)
}

tiles.empty <- function(x) {
  stopifnot(is.tess(x))
  switch(x$type,
         rect = {
           nx <- length(x$xgrid) - 1
           ny <- length(x$ygrid) - 1
           ans <- rep(FALSE, nx * ny)
         },
         tiled = {
           ans <- sapply(x$tiles, is.empty)
         },
         image = {
           ans <- (table(x$image[]) == 0)
         })
  return(ans)
}
           
tilenames <- function(x) {
  stopifnot(is.tess(x))
  switch(x$type,
         rect={
           if(!is.null(x$tilenames)) {
             out <- x$tilenames
           } else {
             nx <- length(x$xgrid) - 1
             ny <- length(x$ygrid) - 1
             ij <- expand.grid(1:nx, 1:ny)
             out <- paste0("Tile row ", ij[,2], ", col ", ij[,1])
           }
         },
         tiled={
           out <- names(x$tiles)
           if(sum(nzchar(out)) != x$n)
             out <- paste("Tile", seq_len(x$n))
         },
         image={
           out <- levels(x$image)
         }
         )
  return(as.character(out))
}

"tilenames<-" <- function(x, value) {
  stopifnot(is.tess(x))
  if(!is.null(value)) {
    ## validate length
    value <- as.character(value)
    nv <- length(value)
    switch(x$type,
           rect = {
             nx <- length(x$xgrid) - 1
             ny <- length(x$ygrid) - 1
             n <- nx * ny
           },
           tiled = { n <- length(x$tiles) },
           image = { n <- length(levels(x$image)) })
    if(nv != n)
      stop("Replacement value has wrong length",
           paren(paste(nv, "instead of", n)))
  }
  switch(x$type,
         rect={
           x$tilenames <- value
         },
         tiled={
           names(x$tiles) <- value
         },
         image={
           levels(x$image) <- value %orifnull% (1:n)
         }
         )
  return(x)
}

marks.tess <- function(x, ...) {
  stopifnot(is.tess(x))
  return(x$marks)
}

"marks<-.tess" <- function(x, ..., value) {
  stopifnot(is.tess(x))
  if(!is.null(value)) {
    value <- as.data.frame(value)
    if(nrow(value) != x$n)
      stop(paste("replacement value for marks has wrong length:",
                 nrow(value), "should be", x$n),
           call.=FALSE)
    rownames(value) <- NULL
    if(ncol(value) == 1) colnames(value) <- "marks"
  }
  x$marks <- value
  return(x)
}

unmark.tess <- function(X) { marks(X) <- NULL; return(X) }

tile.areas <- function(x) {
  stopifnot(is.tess(x))
  switch(x$type,
         rect={
           xg <- x$xgrid
           yg <- x$ygrid
#           nx <- length(xg) - 1 
#           ny <- length(yg) - 1
           a <- outer(rev(diff(yg)), diff(xg), "*")
           a <- as.vector(t(a))
           names(a) <- as.vector(t(tilenames(x)))
         },
         tiled={
           a <- unlist(lapply(x$tiles, area))
         },
         image={
           z <- x$image
           a <- table(z$v) * z$xstep * z$ystep
         })
  return(a)
}

         
as.im.tess <- function(X, W=NULL, ...,
                       eps=NULL, dimyx=NULL, xy=NULL,
                       na.replace=NULL) {
  # if W is present, it may have to be converted
  if(!is.null(W)) {
    stopifnot(is.owin(W))
    if(W$type != "mask")
      W <- as.mask(W, eps=eps, dimyx=dimyx, xy=xy)
  } 
  switch(X$type,
         image={
           out <- as.im(X$image, W=W, eps=eps, dimyx=dimyx, xy=xy,
                        na.replace=na.replace)
         },
         tiled={
           if(is.null(W))
             W <- as.mask(as.owin(X), eps=eps, dimyx=dimyx, xy=xy)
           til <- X$tiles
           ntil <- length(til)
           nama <- names(til)
           if(is.null(nama) || !all(nzchar(nama)))
             nama <- paste(seq_len(ntil))
           xy <- list(x=W$xcol, y=W$yrow)
           for(i in seq_len(ntil)) {
             indic <- as.mask(til[[i]], xy=xy)
             tag <- as.im(indic, value=i)
             if(i == 1) {
               out <- tag
               outv <- out$v
             } else {
               outv <- pmin.int(outv, tag$v, na.rm=TRUE)
             }
           }
           out <- im(factor(outv, levels=seq_len(ntil), labels=nama),
                     out$xcol, out$yrow)
           unitname(out) <- unitname(W)
         },
         rect={
           if(is.null(W))
             out <- as.im(as.rectangle(X), eps=eps, dimyx=dimyx, xy=xy)
           else
             out <- as.im(W)
           xg <- X$xgrid
           yg <- X$ygrid
           nrows <- length(yg) - 1
           ncols <- length(xg) - 1
           jx <- findInterval(out$xcol, xg, rightmost.closed=TRUE)
           iy <- findInterval(out$yrow, yg, rightmost.closed=TRUE)
           M <- as.matrix(out)
           Jcol <- jx[col(M)]
           Irow <- nrows - iy[row(M)] + 1
           Ktile <- Jcol + ncols * (Irow - 1)
           Ktile <- factor(Ktile, levels=seq_len(nrows * ncols))
           out <- im(Ktile, xcol=out$xcol, yrow=out$yrow,
                     unitname=unitname(W))
         }
         )
  return(out)
}

as.function.tess <- function(x, ...) {
  V <- x
  f <- function(x,y) { tileindex(x,y,V) }
  g <- funxy(f, Window(V))
  return(g)
}

tileindex <- function(x, y, Z) {
  stopifnot(is.tess(Z))
  stopifnot(length(x) == length(y))
  switch(Z$type,
         rect={
           jx <- findInterval(x, Z$xgrid, rightmost.closed=TRUE)
           iy <- findInterval(y, Z$ygrid, rightmost.closed=TRUE)
           nrows <- length(Z$ygrid) - 1
           ncols <- length(Z$xgrid) - 1
           iy[iy < 1 | iy > nrows] <- NA
           jx[jx < 1 | jx > ncols] <- NA
           jcol <- jx
           irow <- nrows - iy + 1
           ktile <- jcol + ncols * (irow - 1)
           m <- factor(ktile, levels=seq_len(nrows*ncols))
           ij <- expand.grid(j=seq_len(ncols),i=seq_len(nrows))
           levels(m) <- paste("Tile row ", ij$i, ", col ", ij$j, sep="")
         },
         tiled={
           n <- length(x)
           todo <- seq_len(n)
           nt <- length(Z$tiles)
           m <- integer(n)
           for(i in 1:nt) {
             ti <- Z$tiles[[i]]
             hit <- inside.owin(x[todo], y[todo], ti)
             if(any(hit)) {
               m[todo[hit]] <- i
               todo <- todo[!hit]
             }
             if(length(todo) == 0)
               break
           }
           m[m == 0] <- NA
           nama <- names(Z$tiles)
           lev <- seq_len(nt)
           lab <- if(!is.null(nama) && all(nzchar(nama))) nama else paste("Tile", lev)
           m <- factor(m, levels=lev, labels=lab)
         },
         image={
           Zim <- Z$image
           m <- factor(Zim[list(x=x, y=y), drop=FALSE], levels=levels(Zim))
         }
         )
  return(m)
}
  
as.tess <- function(X) {
  UseMethod("as.tess")
}

as.tess.tess <- function(X) {
  fields <- 
    switch(X$type,
           rect={ c("xgrid", "ygrid") },
           tiled={ "tiles" },
           image={ "image" },
           stop(paste("Unrecognised tessellation type", sQuote(X$type))))
  fields <- c(c("type", "window", "n", "marks"), fields)
  X <- unclass(X)[fields]
  class(X) <- c("tess", class(X))
  return(X)
}

as.tess.im <- function(X) {
  return(tess(image = X))
}

as.tess.list <- function(X) {
  W <- lapply(X, as.owin)
  return(tess(tiles=W))
}

as.tess.owin <- function(X) {
  return(tess(tiles=list(X)))
}

domain.tess <- Window.tess <- function(X, ...) { as.owin(X) } 

intersect.tess <- function(X, Y, ..., keepmarks=FALSE) {
  X <- as.tess(X)
  if(is.owin(Y) && Y$type == "mask") {
    # special case
    # convert to pixel image 
    result <- as.im(Y)
    Xtiles <- tiles(X)
    for(i in seq_along(Xtiles)) {
      tilei <- Xtiles[[i]]
      result[tilei] <- i
    }
    result <- result[Y, drop=FALSE]
    out <- tess(image=result, window=Y)
    if(keepmarks) marks(out) <- marks(X)
    return(out)
  }
  if(is.owin(Y)) {
    # efficient code when Y is a window, retaining names of tiles of X
    Ztiles <- lapply(tiles(X), intersect.owin, B=Y, ..., fatal=FALSE)
    isempty <- sapply(Ztiles, is.empty)
    Ztiles <- Ztiles[!isempty]
    Xwin <- as.owin(X)
    Ywin <- Y
    if(keepmarks) {
      marksX <- marks(X)
      if(!is.null(marksX))
        marx <- as.data.frame(marksX)[!isempty, ]
    }
  } else {
    # general case
    Y <- as.tess(Y)
    Xtiles <- tiles(X)
    Ytiles <- tiles(Y)
    Ztiles <- list()
    namesX <- tilenames(X)
    namesY <- tilenames(Y)
    if(keepmarks) {
      Xmarks <- as.data.frame(marks(X))
      Ymarks <- as.data.frame(marks(Y))
      gotXmarks <- (ncol(Xmarks) > 0)
      gotYmarks <- (ncol(Ymarks) > 0)
      if(gotXmarks && gotYmarks) {
        colnames(Xmarks) <- paste0("X", colnames(Xmarks))
        colnames(Ymarks) <- paste0("Y", colnames(Ymarks))
      }
      if(gotXmarks || gotYmarks) {
        marx <- if(gotXmarks && gotYmarks) {
          cbind(Xmarks[integer(0), , drop=FALSE],
                Ymarks[integer(0), , drop=FALSE])
        } else if(gotXmarks) {
          Xmarks[integer(0), , drop=FALSE]
        } else {
          Ymarks[integer(0), , drop=FALSE]
        }
      } else keepmarks <- FALSE
    }
    for(i in seq_along(Xtiles)) {
      Xi <- Xtiles[[i]]
      Ti <- lapply(Ytiles, intersect.owin, B=Xi, ..., fatal=FALSE)
      isempty <- sapply(Ti, is.empty)
      nonempty <- !isempty
      if(any(nonempty)) {
        Ti <- Ti[nonempty]
        names(Ti) <- paste(namesX[i], namesY[nonempty], sep="x")
        Ztiles <- append(Ztiles, Ti)
        if(keepmarks) {
          extra <- if(gotXmarks && gotYmarks) {
            data.frame(X=Xmarks[i, ,drop=FALSE],
                       Y=Ymarks[nonempty, ,drop=FALSE],
                       row.names=NULL)
          } else if(gotYmarks) {
            Ymarks[nonempty, ,drop=FALSE]
          } else {
            Xmarks[rep(i, sum(nonempty)), ,drop=FALSE]
          }
          marx <- rbind(marx, extra)
        }
      }
    }
    Xwin <- as.owin(X)
    Ywin <- as.owin(Y)
  }
  Zwin <- intersect.owin(Xwin, Ywin)
  out <- tess(tiles=Ztiles, window=Zwin)
  if(keepmarks) 
    marks(out) <- marx
  return(out)
}


bdist.tiles <- local({

  vdist <- function(x,w) {
    z <- as.ppp(vertices(x), W=w, check=FALSE)
    min(bdist.points(z))
  }
  edist <- function(x,b) {
    xd <- crossdist(edges(x, check=FALSE), b, type="separation")
    min(xd)
  }

  bdist.tiles <-  function(X) {
    if(!is.tess(X))
      stop("X must be a tessellation")
    W <- as.owin(X)
    switch(X$type,
           rect=,
           tiled={
             tt <- tiles(X)
             if(is.convex(W)) {
               # distance is minimised at a tile vertex
               d <- sapply(tt, vdist, w=W)
             } else {
               # coerce everything to polygons
               W  <- as.polygonal(W)
               tt <- lapply(tt, as.polygonal)
               # compute min dist from tile edges to window edges
               d <- sapply(tt, edist, b=edges(W))
             }
           },
           image={
             Xim <- X$image
             # compute boundary distance for each pixel
             bd <- bdist.pixels(as.owin(Xim), style="image")
             bd <- bd[W, drop=FALSE]
             # split over tiles
             bX <- split(bd, X)
             # compute minimum distance over each level of factor
             d <- sapply(bX, function(z) { summary(z)$min })
           }
           )
    return(d)
  }
  bdist.tiles
})


## ......... geometrical transformations ..................

shift.tess <- function(X, ...) {
  Y <- X
  Y$window <- wY <- shift(X$window, ...)
  vec <- getlastshift(wY)
  switch(X$type,
         rect={
           Y$xgrid <- Y$xgrid + vec[1]
           Y$ygrid <- Y$ygrid + vec[2]
         },
         tiled={
           Y$tiles <- lapply(Y$tiles, shift, vec=vec)
         },
         image = {
           Y$image <- shift(Y$image, vec)
         })
  attr(Y, "lastshift") <- vec
  return(Y)
}

affine.tess <- function(X, mat=diag(c(1,1)), vec=c(0,0), ...) {
  Y <- X
  Y$window <- affine(X$window, mat=mat, vec=vec, ...)
  switch(Y$type,
         rect = {
           if(all(mat == diag(diag(mat)))) {
             ## result is rectangular
             Y$xgrid <- sort(mat[1,1] * X$xgrid + vec[1])
             Y$ygrid <- sort(mat[2,2] * X$ygrid + vec[2])
           } else {
             ## shear transformation; treat rectangles as general tiles
             Y <- tess(tiles=tiles(X), window=Y$window)
             Y$tiles <- lapply(Y$tiles, affine, mat=mat, vec=vec, ...)
           }
         },
         tiled={
           Y$tiles <- lapply(Y$tiles, affine, mat=mat, vec=vec, ...)
         },
         image = {
           Y$image <- affine(Y$image, mat=mat, vec=vec, ...)
         })
  return(Y)
}

reflect.tess <- function(X) {
  Y <- X
  Y$window <- reflect(Y$window)
  switch(X$type,
         rect = {
           Y$xgrid <- rev(- Y$xgrid)
           Y$ygrid <- rev(- Y$ygrid)
         },
         tiled = {
           Y$tiles <- lapply(Y$tiles, reflect)
         },
         image = {
           Y$image <- reflect(Y$image)
         })
  return(Y)
}

scalardilate.tess <- function(X, f, ...) {
  Y <- X
  Y$window <- scalardilate(X$window, f, ...)
  switch(X$type,
         rect = {
           Y$xgrid <- f * Y$xgrid
           Y$ygrid <- f * Y$ygrid
         },
         tiled = {
           Y$tiles <- lapply(Y$tiles, scalardilate, f=f, ...)
         },
         image = {
           Y$image <- scalardilate(Y$image, f=f, ...)
         })
  return(Y)
}

rotate.tess <- function(X, angle=pi/2, ..., centre=NULL) {
  if(angle %% (2 * pi) == 0) return(X)
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  Y <- X
  Y$window <- rotate(X$window, angle=angle, ...)
  switch(X$type,
         rect = {
           if(angle %% (pi/2) == 0) {
             ## result is rectangular
             co <- round(cos(angle))
             si <- round(sin(angle))
             Y$xgrid <- sort((if(co == 0) 0 else (co * X$xgrid)) -
                             (if(si == 0) 0 else (si * X$ygrid)))
             Y$ygrid <- sort((if(si == 0) 0 else (si * X$xgrid)) +
                             (if(co == 0) 0 else (co * X$ygrid)))
           } else {
             ## general tessellation
             Y <- tess(tiles=lapply(tiles(X), rotate, angle=angle, ...),
                       window=Y$window)
           }
         },
         tiled = {
           Y$tiles <- lapply(X$tiles, rotate, angle=angle, ...)
         },
         image = {
           Y$image <- rotate(X$image, angle=angle, ...)
         })
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}
  
