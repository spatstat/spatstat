#' Calculating Laslett's transform
#' Original by Kassel Hingee
#' Adapted by Adrian Baddeley
#' Copyright (C) 2016 Kassel Hingee and Adrian Baddeley

# $Revision: 1.9 $  $Date: 2020/03/16 10:28:51 $

laslett <- function(X, ...,
                    verbose=FALSE, plotit=TRUE,
                    discretise=FALSE,
                    type = c("lower", "upper", "left", "right")){
  #' validate X and convert to a logical matrix
  type <- match.arg(type)
  oldX <- X

  if(is.im(X)) {
    X <- solutionset(X != 0)
  } else if(!is.owin(X)) 
    stop("X should be an image or a window", call.=FALSE)

  if(type != "lower") {
    nrot <- match(type, c("right", "upper", "left"))
    theta <- nrot * pi/2
    X <- rotate(X, angle=-theta)
  }
  
  if(!discretise && (is.polygonal(X) || is.rectangle(X))) {
    result <- polyLaslett(X, ..., oldX=oldX, verbose=verbose, plotit=FALSE)
  } else {
    result <- maskLaslett(X, ..., oldX=oldX, verbose=verbose, plotit=FALSE)
  }

  if(type != "lower") {
    #' rotate back
    prods <- c("TanOld", "TanNew", "Rect")
    result[prods] <- lapply(result[prods], rotate, angle=theta)
  }

  if(plotit)
    plot(result, ...)

  result$type <- type
  
  return(result)
}

maskLaslett <- local({
  
  sumtoright <- function(x) { rev(cumsum(rev(x))) - x }

  maskLaslett <- function(X, ...,
                          eps=NULL, dimyx=NULL, xy=NULL, 
                          oldX=X, verbose=FALSE, plotit=TRUE) {
    if(is.null(oldX)) oldX <- X
    X <- as.mask(X, eps=eps, dimyx=dimyx, xy=xy)
    unitX <- unitname(X)
    if(is.empty(X))
      stop("Empty window!")
    M <- as.matrix(X)
    #' ....... Compute transformed set ...................
    #' Total width of transformed set on each row
    TotFalse <- rowSums(!M)
    ## compute transformed set
    Laz <- (col(M) <= TotFalse[row(M)])
    Laz <- owin(mask=Laz, xrange=X$xrange, yrange=X$yrange, unitname=unitX)
    #' Largest sub-rectangle of transformed set
    width <- min(TotFalse) * X$xstep
    Rect <- owin(X$xrange[1L] + c(0, width), X$yrange, unitname=unitX)
    #' Along each horizontal line (row),
    #' compute a running count of FALSE pixels.
    #' This is the mapping for the set transform
    #' (the value at any pixel gives the new column number
    #' for the transformed pixel)
    CumulFalse <- t(apply(!M, 1L, cumsum))
    #' discard one column for consistency with other matrices below
    CumulFalse <- CumulFalse[,-1L,drop=FALSE]

    #' ....... Find lower tangent points .................

    #' compute discrete gradient in x direction
    G <- t(apply(M, 1, diff))
    #' detect entries, exits, changes
    Exit <- (G == -1)
    Enter <- (G == 1)
    Change <- Exit | Enter
    #' form a running total of the number of pixels inside X
    #' to the **right** of the current pixel
    FutureInside <- t(apply(M, 1, sumtoright))[,-1L,drop=FALSE]
    #' find locations of changes 
    loc <- which(Change, arr.ind=TRUE)
    #' don't consider entries/exits in the bottom row
    ok <- (loc[,"row"] > 1) 
    loc <- loc[ok, , drop=FALSE]
    #' corresponding locations on horizontal line below current line
    below <- cbind(loc[,"row"]-1L, loc[,"col"])
    #' look up data at these locations
    df <- data.frame(row=loc[,"row"],
                     col=loc[,"col"],
                     newcol=CumulFalse[loc],
                     Exit=Exit[loc],
                     Enter=Enter[loc],
                     InsideBelow=M[below],
                     FutureInsideBelow=FutureInside[below])
    #' identify candidates for tangents
    df$IsCandidate <-
      with(df, Enter & !InsideBelow & (newcol < TotFalse[row]))
    #' collect data for each horizontal line (row)
    #' then sort by increasing x (column) within each line.
    oo <- with(df, order(row, col))
    df <- df[oo, , drop=FALSE]
    #' divide data into one piece for each hztal line
    g <- split(df, df$row)
    #' Initialise empty list of tangent points
    tangents <- data.frame(row=integer(0), col=integer(0), newcol=integer(0))
    #' process each hztal line
    for(p in g) {
      tangents <-
        with(p, {
          candidates <- which(IsCandidate) # indices are row numbers in 'p'
          if(verbose) print(p)
          exits <- which(Exit)
          for(i in candidates) {
            if(verbose) cat(paste("candidate", i, "\n"))
            if(any(found <- (exits > i))) {
              j <- exits[min(which(found))]
              if(verbose) cat(paste("next exit:", j, "\n"))
              #' check no pixels inside X in row below between i and j
              if(FutureInsideBelow[i] == FutureInsideBelow[j]) {
                if(verbose)
                  cat(paste("Tangent (1) at row=", row[i],
                            "col=", col[i], "\n"))
                tangents <- rbind(tangents,
                                  data.frame(row=row[i],
                                             col=col[i],
                                             newcol=newcol[i]))
              }
            } else {
              #' no exits on this row
              if(verbose)
                cat("no subsequent exit\n")
              if(FutureInsideBelow[i] == 0) {
                if(verbose)
                  cat(paste("Tangent (2) at row=", row[i],
                            "col=", col[i], "\n"))
                tangents <- rbind(tangents,
                                  data.frame(row=row[i],
                                             col=col[i],
                                             newcol=newcol[i]))
              }
            }
          }
          if(verbose) cat("====\n")
          tangents
        })
    }
    tangents$oldx <- X$xcol[tangents$col]
    tangents$newx <- X$xcol[tangents$newcol]
    tangents$y    <- X$yrow[tangents$row]
    TanOld <- with(tangents, ppp(oldx, y, window=Frame(X), unitname=unitX))
    TanNew <- with(tangents, ppp(newx, y, window=Laz), unitname=unitX)
    result <- list(oldX=oldX,
                   TanOld=TanOld, TanNew=TanNew, Rect=Rect,
                   df=tangents)
    class(result) <- c("laslett", class(result))
    if(plotit)
      plot(result, ...)
    return(result)
  }

  maskLaslett
})

print.laslett <- function(x, ...) {
  cat("Laslett Transform\n")
  cat("\nOriginal object:\n")
  print(x$oldX)
  cat("\nTransformed set:\n")
  W <- Window(x$TanNew)
  print(W)
  unitinfo <- summary(unitname(W))
  cat("\nTransformed area:", area.owin(W),
      "square", unitinfo$plural, unitinfo$explain,
      fill=TRUE)
  cat("\n")
  type <- x$type %orifnull% "lower"
  cat(npoints(x$TanNew), type, "tangent points found.", fill=TRUE)
  return(invisible(NULL))
}
  
plot.laslett <- function(x, ...,
                         Xpars=list(box=TRUE, col="grey"),
                         pointpars=list(pch=3, cols="blue"),
                         rectpars=list(lty=3, border="green")) {
  Display <- with(x,
                  solist(Original=
                         layered(oldX,
                                 TanOld,
                                 plotargs=list(Xpars, pointpars)),
                         Transformed=
                         layered(TanNew,
                                 Rect, 
                                 plotargs=list(pointpars, rectpars))))

  #' ignore arguments intended for as.mask
  argh <- list(...)
  if(any(bad <- names(argh) %in% c("eps", "dimyx", "xy")))
    argh <- argh[!bad]

  do.call(plot,
          resolve.defaults(list(x=Display),
                           argh,
                           list(main="", mar.panel=0, hsep=1,
                                equal.scales=TRUE)))
  return(invisible(NULL))
}

polyLaslett <- function(X, ..., oldX=X, verbose=FALSE, plotit=TRUE) {
  X <- as.polygonal(X)
  if(is.empty(X))
    stop("Empty window!")
  unitX <- unitname(X)
  # expand frame slightly
  B <- Frame(X)
  B <- grow.rectangle(B, max(sidelengths(B))/8)
  x0 <- B$xrange[1L]
  x1 <- B$xrange[2L]
  # extract vertices
  v <- vertices(X)
  nv <- length(v$x)
  # ..........  compute transformed set .....................
  # make horizontal segments from each vertex to sides of box
  left <- with(v, psp(rep(x0,nv), y, x, y, window=B, marks=1:nv, check=FALSE))
  right <- with(v, psp(x, y, rep(x1,nv), y, window=B, marks=1:nv, check=FALSE))
  # intersect each horizontal segment with the window
  if(verbose) cat("Processing", nv, "polygon vertices... ")
  clipleft <- clip.psp(left, X)
  clipright <- clip.psp(right, X)
  if(verbose) cat("Done.\n")
  # calculate lengths of clipped segments, and group by vertex.
  # marks indicate which hztal segment was the parent of each piece.
  lenleft <- tapply(lengths_psp(clipleft),
                    factor(marks(clipleft), levels=1:nv),
                    sum)
  lenright <- tapply(lengths_psp(clipright),
                     factor(marks(clipright), levels=1:nv),
                     sum)
  lenleft[is.na(lenleft)] <- 0
  lenright[is.na(lenright)] <- 0
  emptylenleft <- lengths_psp(left) - lenleft
  emptylenright <- lengths_psp(right) - lenright
  # The transformed polygon 
  isrightmost <- (lenright == 0)
  yright <- v$y[isrightmost]
  xright <- x0 + (emptylenleft+emptylenright)[isrightmost]
  minxright <- min(xright) # right margin of largest rectangle
  ord <- order(yright)
  Ty <- yright[ord]
  Tx <- xright[ord]
  nT <- length(Ty)
  if(Tx[nT] > x0) {
    Ty <- c(Ty, Ty[nT])
    Tx <- c(Tx, x0)
  }
  if(Tx[1L] > x0) {
    Ty <- c(Ty[1L], Ty)
    Tx <- c(x0,    Tx)
  }
  TX <- owin(B$xrange, B$yrange, poly=list(x=Tx, y=Ty), check=FALSE)
  TX <- TX[Frame(X)]
  # ..........  identify lower tangents .....................
  V <- as.ppp(v, W=Frame(X), unitname=unitX)
  is.candidate <- is.tangent <- logical(nv)
  # apply simple criteria for ruling in or out
  Plist <- X$bdry
  cumnv <- 0
  for(i in seq_along(Plist)) {
    P <- Plist[[i]]
    xx <- P$x
    yy <- P$y
    nn <- length(xx)
#    xnext <- c(xx[-1L], xx[1L])
    ynext <- c(yy[-1L], yy[1L])
#    xprev <- c(xx[nn], xx[-nn])
    yprev <- c(yy[nn], yy[-nn])
    is.candidate[cumnv + seq_len(nn)] <- 
      if(!is.hole.xypolygon(P)) {
        (yprev > yy & ynext >= yy)
      } else {
        (yprev >= yy & ynext > yy)
      }
    cumnv <- cumnv + nn
  }
  ##  was.candidate <- is.candidate
  
  #' reject candidates lying too close to boundary
  tooclose <- (bdist.points(V[is.candidate]) < diameter(Frame(V))/1000)
  is.candidate[is.candidate][tooclose] <- FALSE
  #' evaluate candidate points
  #' make tiny boxes around vertex
  candidates <- which(is.candidate)
  nc <- length(candidates)
  nnd <- nndist(V)
  if(verbose) {
    cat(paste("Processing", nc, "tangent candidates ... "))
    pstate <- list()
  }
  tiny <- .Machine$double.eps
  for(j in 1:nc) {
    i <- candidates[j]
    eps <- nnd[i]/16
    xi <- v$x[i]
    yi <- v$y[i]
    Below <- owin(xi + c(-eps,eps), yi + c(-eps, 0))
#    Above <- owin(xi + c(-eps, eps), yi + c(0, eps))
    UpLeft <- owin(xi + c(-eps, 0), yi + c(0, eps))
    is.tangent[i] <- (overlap.owin(X, Below) <= tiny) &&
                     (overlap.owin(X, UpLeft) < eps^2)
    if(verbose)
      pstate <- progressreport(j, nc, state=pstate)
  }
  if(verbose) cat(paste("Found", sum(is.tangent), "tangents\n"))
  TanOld <- V[is.tangent]
  ynew <- TanOld$y
  xnew <- x0 + emptylenleft[is.tangent]
  TanNew <- ppp(xnew, ynew, window=TX, check=FALSE, unitname=unitX)
  #  maximal rectangle
  Rect <- owin(c(X$xrange[1L], minxright), X$yrange, unitname=unitX)
  #
  df <- data.frame(xold=TanOld$x, xnew=TanNew$x, y=TanNew$y)
  #
  result <- list(oldX=oldX,
                 TanOld=TanOld, TanNew=TanNew, Rect=Rect,
                 df=df)
  class(result) <- c("laslett", class(result))
  if(plotit)
    plot(result, ...)
  return(result)
}

