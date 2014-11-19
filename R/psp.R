#
#  psp.R
#
#  $Revision: 1.78 $ $Date: 2014/11/19 07:19:42 $
#
# Class "psp" of planar line segment patterns
#
#
#################################################
# creator
#################################################
psp <- function(x0, y0, x1, y1, window, marks=NULL,
                check=spatstat.options("checksegments")) {
  stopifnot(is.numeric(x0))
  stopifnot(is.numeric(y0))
  stopifnot(is.numeric(x1))
  stopifnot(is.numeric(y1))
  stopifnot(is.vector(x0))
  stopifnot(is.vector(y0))
  stopifnot(is.vector(x1))
  stopifnot(is.vector(y1))
  stopifnot(length(x0) == length(y0))
  stopifnot(length(x1) == length(y1))
  stopifnot(length(x0) == length(x1))
  ends <- data.frame(x0=x0,y0=y0,x1=x1,y1=y1)
  if(!missing(window))
    verifyclass(window,"owin")
  if(check) {
    ok <- inside.owin(x0,y0, window) & inside.owin(x1,y1,window)
    if((nerr <- sum(!ok)) > 0)
      stop(paste(nerr, ngettext(nerr, "segment does not", "segments do not"),
                 "lie entirely inside the window.\n"), call.=FALSE)
  }
  out <- list(ends=ends,
              window=window,
              n = nrow(ends))

# add marks if any
  if(!is.null(marks)) {
    if(is.matrix(marks))
      marks <- as.data.frame(marks)
    if(is.data.frame(marks)) {
      omf <- "dataframe"
      nmarks <- nrow(marks)
      rownames(marks) <- seq_len(nmarks)
      whinge <- "The number of rows of marks"
    } else {
      omf <- "vector"
      names(marks) <- NULL
      nmarks <- length(marks)
      whinge <- "The length of the marks vector"
    }
    if(nmarks != out$n) stop(paste(whinge, "!= length of x and y.\n"))
    out$marks <- marks
    out$markformat <- omf
  } else {
    out$markformat <- "none"
  }

  class(out) <- c("psp", class(out))
  return(out)
}

######################################################
#  conversion
######################################################

is.psp <- function(x) { inherits(x, "psp") }

as.psp <- function(x, ..., from=NULL, to=NULL) {
  # special case: two point patterns
  if(is.null(from) != is.null(to))
    stop(paste("If one of", sQuote("from"), "and", sQuote("to"),
               "is specified, then both must be specified.\n"))
  if(!is.null(from) && !is.null(to)) {
    verifyclass(from, "ppp")
    verifyclass(to, "ppp")
    if(from$n != to$n)
      stop(paste("The point patterns", sQuote("from"), "and", sQuote("to"),
                 "have different numbers of points.\n"))
    uni <- union.owin(from$window, to$window)
    Y <- do.call("psp",
                 resolve.defaults(list(from$x, from$y, to$x, to$y),
                                  list(...),
                                  list(window=uni)))
    return(Y)
  }
  UseMethod("as.psp")
}

as.psp.psp <- function(x, ..., check=FALSE, fatal=TRUE) {
  if(!verifyclass(x, "psp", fatal=fatal))
    return(NULL)
  ends <- x$ends
  psp(ends$x0, ends$y0, ends$x1, ends$y1, window=x$window,
      marks=x$marks, check=check)
}

as.psp.data.frame <- function(x, ..., window=NULL, marks=NULL,
                              check=spatstat.options("checksegments"), fatal=TRUE) {
  window <- suppressWarnings(as.owin(window,fatal=FALSE))
  if(!is.owin(window)) {
    if(fatal) stop("Cannot interpret \"window\" as an object of class owin.\n")
    return(NULL)
  }

  if(checkfields(x,"marks")) {
    if(is.null(marks)) marks <- x$marks
    else warning(paste("Column named \"marks\" ignored;\n",
                       "argument named \"marks\" has precedence.\n",sep=""))
    x$marks <- NULL
  }

  if(checkfields(x, c("x0", "y0", "x1", "y1"))) {
    out <- psp(x$x0, x$y0, x$x1, x$y1, window=window,
               check=check)
    x <- x[-match(c("x0","y0","x1","y1"),names(x))]
  }
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    rr <- x$length/2
    dx <- cos(x$angle) * rr
    dy <- sin(x$angle) * rr
    bb <- boundingbox(window)
    rmax <- max(rr)
    bigbox <- owin(bb$xrange + c(-1,1) * rmax, bb$yrange + c(-1,1) * rmax)
    pattern <- psp(x$xmid - dx, x$ymid - dy, x$xmid + dx, x$ymid + dy,
                   window=bigbox,check=FALSE)
    out <- pattern[window]
    x <- x[-match(c("xmid","ymid","length","angle"),names(x))]
  }
  else if(ncol(x) >= 4) {
    out <- psp(x[,1], x[,2], x[,3], x[,4], window=window,
               check=check)
    x <- x[-(1:4)]
  }
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern.", call.=FALSE)
  else out <- NULL

  if(!is.null(out)) {
    if(is.null(marks) & ncol(x) > 0) marks <- x
    if(is.null(marks)) {
       out$markformat <- "none"
    } else {
       out$marks <- marks
       out$markformat <- if(is.data.frame(marks)) "dataframe" else "vector"
       out <- as.psp(out,check=FALSE)
    }
  }
  return(out)
}

as.psp.matrix <- function(x, ..., window=NULL, marks=NULL,
                          check=spatstat.options("checksegments"), fatal=TRUE) {
   x <- as.data.frame(x)
   as.psp(x,...,window=window,marks=marks,check=check,fatal=fatal)
}

as.psp.default <- function(x, ..., window=NULL, marks=NULL,
                           check=spatstat.options("checksegments"), fatal=TRUE) {
  if(checkfields(x,"marks")) {
	if(is.null(marks)) marks <- x$marks
	else warning(paste("Component of \"x\" named \"marks\" ignored;\n",
                             "argument named \"marks\" has precedence.\n",sep=""))
  }
  if(checkfields(x, c("x0", "y0", "x1", "y1")))
    return(psp(x$x0, x$y0, x$x1, x$y1, window=window, marks=marks,
               check=check))
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    rr <- x$length/2
    dx <- cos(x$angle) * rr
    dy <- sin(x$angle) * rr
    window <- as.owin(window)
    bb <- boundingbox(window)
    rmax <- max(rr)
    bigbox <- owin(bb$xrange + c(-1,1) * rmax, bb$yrange + c(-1,1) * rmax)
    pattern <- psp(x$x - dx, x$y - dy, x$x + dx, x$y + dy,
                   window=bigbox, marks=marks, check=FALSE)
    clipped <- pattern[window]
    return(clipped)
  }
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern")
  return(NULL)
}

as.psp.owin <- function(x, ..., window=NULL,
                        check=spatstat.options("checksegments"), fatal=TRUE) {
  .Deprecated("edges", package="spatstat")
  edges(x, ..., window=window, check=check)
}

edges <- function(x, ...,
                  window=NULL, check=FALSE) {
  x <- as.owin(x)
  if(is.null(window)) window <- as.rectangle(x)
  x <- as.polygonal(x)
  x0 <- y0 <- x1 <- y1 <- numeric(0)
  bdry <- x$bdry
  for(i in seq_along(bdry)) {
    po <- bdry[[i]]
    ni <- length(po$x)
    nxt <- c(2:ni, 1)
    x0 <- c(x0, po$x)
    y0 <- c(y0, po$y)
    x1 <- c(x1, po$x[nxt])
    y1 <- c(y1, po$y[nxt])
  }
  out <- psp(x0, y0, x1, y1,  window=window, check=check)
  return(out)
}


#################

as.data.frame.psp <- function(x, row.names=NULL, ...) {
  df <- as.data.frame(x$ends, row.names=row.names)
  if(is.marked(x))
    df <- cbind(df, if(x$markformat=="dataframe") marks(x)
                    else data.frame(marks=marks(x)))
  return(df)
}

#######  manipulation ##########################

append.psp <- function(A,B) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  stopifnot(identical(A$window, B$window))
  marks <- marks(A) %mapp% marks(B)
  ends <- rbind(A$ends, B$ends)
  out  <- as.psp(ends,window=A$window,marks=marks,check=FALSE)
  return(out)
}

rebound.psp <- function(x, rect) {
  verifyclass(x, "psp")
  x$window <- rebound.owin(x$window, rect)
  return(x)
}


#################################################
#  marks
#################################################

is.marked.psp <- function(X, ...) {
  marx <- marks(X, ...)
  return(!is.null(marx))
}

marks.psp <- function(x, ..., dfok = TRUE) {
  # data frames of marks are as of 19/March 2011 implemented for psp
    ma <- x$marks
    if ((is.data.frame(ma) || is.matrix(ma)) && !dfok) 
        stop("Sorry, not implemented when the marks are a data frame.\n")
    return(ma)
}

"marks<-.psp" <- function(x, ..., value) {
  stopifnot(is.psp(x))
  if(is.null(value)) {
    return(unmark(x))
  }
  m <- value
  if(!(is.vector(m) || is.factor(m) || is.data.frame(m) || is.matrix(m)))
    stop("Incorrect format for marks")

    if (is.hyperframe(m)) 
        stop("Hyperframes of marks are not supported in psp objects.\n")
    nseg <- nsegments(x)
    if (!is.data.frame(m) && !is.matrix(m)) {
        if (length(m) == 1) 
            m <- rep.int(m, nseg)
        else if (nseg == 0) 
            m <- rep.int(m, 0)
        else if (length(m) != nseg) 
            stop("Number of marks != number of line segments.\n")
        marx <- m
    }
    else {
        m <- as.data.frame(m)
        if (ncol(m) == 0) {
            marx <- NULL
        }
        else {
            if (nrow(m) == nseg) {
                marx <- m
            }
            else {
                if (nrow(m) == 1 || nseg == 0) {
                  marx <- as.data.frame(lapply(as.list(m),function(x,k) {
                    rep.int(x, k)}, k = nseg))
                }
                else stop("Number of rows of data frame != number of points.\n")
            }
        }
    }
    Y <- as.psp(x$ends, window = x$window, marks = marx, check = FALSE)
    return(Y)
}

markformat.psp <- function(x) {
    mf <- x$markformat
    if(is.null(mf)) 
      mf <- markformat(marks(x))
    return(mf)
}

unmark.psp <- function(X) {
  X$marks <- NULL
  X$markformat <- "none"
  return(X)
}

#################################################
#  plot and print methods
#################################################

plot.psp <- function(x, ..., main, add=FALSE, show.all=!add, which.marks=1,
                     ribbon=show.all, ribsep=0.15, ribwid=0.05, ribn=1024,
                     do.plot=TRUE) {
  if(missing(main) || is.null(main))
    main <- short.deparse(substitute(x))
  verifyclass(x, "psp")
  #
  n <- nsegments(x)
  marx <- marks(x)
  #
  use.colour <- !is.null(marx) && (n != 0)
  do.ribbon <- identical(ribbon, TRUE) && use.colour 
  ##
  ## ....   initialise plot; draw observation window  ......
  owinpars <- setdiff(graphicsPars("owin"), "col")
  if(!do.ribbon) {
    ## window of x only
    bb.all <- as.rectangle(as.owin(x))
    if(do.plot && show.all)
      do.call.plotfun("plot.owin", 
                      resolve.defaults(list(x=x$window, main=main,
                                            add=add, show.all=show.all),
                                       list(...)),
                      extrargs=owinpars)
  } else {
    ## enlarged window with room for colour ribbon
    ## x at left, ribbon at right
    bb <- as.rectangle(as.owin(x))
    xwidth <- diff(bb$xrange)
    xheight <- diff(bb$yrange)
    xsize <- max(xwidth, xheight)
    bb.rib <- owin(bb$xrange[2] + c(ribsep, ribsep+ribwid) * xsize,
                   bb$yrange)
    bb.all <- boundingbox(bb.rib, bb)
    if(do.plot) {
      pt <- prepareTitle(main)
      ## establish coordinate system
      if(!add)
      do.call.plotfun("plot.owin",
                      resolve.defaults(list(x=bb.all,
                                            type="n",
                                            main=pt$blank),
                                       list(...)),
                      extrargs=owinpars)
      ## now plot window of x
      ## with title centred on this window
      if(show.all) {
        do.call.plotfun("plot.owin", 
                        resolve.defaults(list(x=x$window,
                                              add=TRUE,
                                              main=main,
                                              show.all=TRUE),
                                         list(...)),
                        extrargs=owinpars)
        ## title done. 
        main <- ""
      }
    }
  }

  # plot segments
  if(n == 0) {
    result <- symbolmap()
    attr(result, "bbox") <- bb.all
    return(invisible(result))
  }
  
  # determine colours if any
  if(!use.colour) {
    # black
    col <- colmap <- NULL
  } else {
    # multicoloured 
    marx <- as.data.frame(marx)[, which.marks]
    if(is.character(marx) || length(unique(marx)) == 1)
      marx <- factor(marx)
    if(is.factor(marx)) {
      lev <- levels(marx)
      colmap <- colourmap(col=rainbow(length(lev)), inputs=factor(lev))
    } else {
      if(!all(is.finite(marx)))
        warning("Some mark values are infinite or NaN or NA")
      colmap <- colourmap(col=rainbow(ribn), range=range(marx, finite=TRUE))
    }
    col <- colmap(marx)
  }

  ## convert to greyscale?
  if(spatstat.options("monochrome")) {
    col <- to.grey(col)
    colmap <- to.grey(colmap)
  }

  if(do.plot) {
    ## plot segments
    do.call.plotfun("segments",
                    resolve.defaults(as.list(x$ends),
                                     list(...),
                                     list(col=col),
                                     .StripNull=TRUE),
                    extrargs=names(par()))
    ## plot ribbon
    if(do.ribbon) 
      plot(colmap, vertical=TRUE, add=TRUE,
           xlim=bb.rib$xrange, ylim=bb.rib$yrange)
  }
  
  # return colour map
  result <- colmap %orifnull% colourmap()
  attr(result, "bbox") <- bb.all
  return(invisible(result))
}

print.psp <- function(x, ...) {
  verifyclass(x, "psp")
  n <- x$n
  ism <- is.marked(x, dfok = TRUE)
  splat(if(ism) "marked" else NULL,
        "planar line segment pattern:",
        n, ngettext(n, "line segment", "line segments"))
  if(ism) {
    mks <- marks(x, dfok = TRUE)
    if(is.data.frame(mks)) {
      splat("Mark variables: ",
            paste(names(mks), collapse = ", "))
    } else {
      if(is.factor(mks)) {
        splat("multitype, with levels =",
              paste(levels(mks), collapse = "\t"))
      } else {
        splat("marks are",
              if(is.numeric(mks)) "numeric," else NULL,
              "of type", sQuote(typeof(mks)))
      }
    }
  }
  print(x$window)
  return(invisible(NULL))
}

unitname.psp <- function(x) {
  return(unitname(x$window))
}

"unitname<-.psp" <- function(x, value) {
  w <- x$window
  unitname(w) <- value
  x$window <- w
  return(x)
}

####################################################
#    summary information
####################################################

endpoints.psp <- function(x, which="both") {
  verifyclass(x, "psp")
  ends <- x$ends
  n <- x$n
  switch(which,
         both={
           first <- second <- rep.int(TRUE, n)
         },
         first={
           first <- rep.int(TRUE, n)
           second <- rep.int(FALSE, n)
         },
         second={
           first <- rep.int(FALSE, n)
           second <- rep.int(TRUE, n)
         },
         left={
           first <- (ends$x0 < ends$x1)
           second <- !first
         },
         right={
           first <- (ends$x0 > ends$x1)
           second <- !first
         },
         lower={
           first <- (ends$y0 < ends$y1)
           second <- !first
         },
         upper={
           first <- (ends$y0 > ends$y1)
           second <- !first
         },
         stop(paste("Unrecognised option: which=", sQuote(which)))
         )
  ok <- rbind(first, second)
  xmat <- rbind(ends$x0, ends$x1)
  ymat <- rbind(ends$y0, ends$y1)
  idmat <- col(ok)
  xx <- as.vector(xmat[ok])
  yy <- as.vector(ymat[ok])
  id <- as.vector(idmat[ok])
  result <- ppp(xx, yy, window=x$window, check=FALSE)
  attr(result, "id") <- id
  return(result)
}

midpoints.psp <- function(x) {
  verifyclass(x, "psp")
  xm <- eval(expression((x0+x1)/2), envir=x$ends)
  ym <- eval(expression((y0+y1)/2), envir=x$ends)
  win <- x$window
  ok <- inside.owin(xm, ym, win)
  if(any(!ok)) {
    warning(paste("Some segment midpoints lie outside the original window;",
                  "window replaced by bounding box"))
    win <- boundingbox(win)
  }
  ppp(x=xm, y=ym, window=win, check=FALSE)
}

lengths.psp <- function(x) {
  verifyclass(x, "psp")
  eval(expression(sqrt((x1-x0)^2 + (y1-y0)^2)), envir=x$ends)
}

angles.psp <- function(x, directed=FALSE) {
  verifyclass(x, "psp")
  a <- eval(expression(atan2(y1-y0, x1-x0)), envir=x$ends)
  if(!directed) 
    a <- a %% pi
  return(a)
}

summary.psp <- function(object, ...) {
  verifyclass(object, "psp")
  len <- lengths.psp(object)
  out <- list(n = object$n,
              len = summary(len),
              totlen = sum(len),
              ang= summary(angles.psp(object)),
              w = summary.owin(object$window),
              marks=if(is.null(object$marks)) NULL else summary(object$marks),
              unitinfo=summary(unitname(object)))
  class(out) <- c("summary.psp", class(out))
  return(out)
}

print.summary.psp <- function(x, ...) {
  cat(paste(x$n, "line segments\n"))
  cat("Lengths:\n")
  print(x$len)
  unitblurb <- paste(x$unitinfo$plural, x$unitinfo$explain)
  cat(paste("Total length:", x$totlen, unitblurb, "\n"))
  cat(paste("Length per unit area:", x$totlen/x$w$area, "\n"))
  cat("Angles (radians):\n")
  print(x$ang)
  print(x$w)
  if(!is.null(x$marks)) {
    cat("Marks:\n")
    print(x$marks)
  }
  return(invisible(NULL))
}

  
########################################################
#  subsets
########################################################

"[.psp" <-
  function(x, i, j, drop, ...) {

    verifyclass(x, "psp")
    
    if(missing(i) && missing(j))
      return(x)
        
    if(!missing(i)) {
      style <- if(inherits(i, "owin")) "window" else "index"
      switch(style,
             window={
               x <- clip.psp(x, window=i, check=FALSE)
             },
             index={
               enz <- x$ends[i, ]
               win <- x$window
               marx <- marksubset(x$marks, i, markformat(x))
               x <- with(enz, psp(x0, y0, x1, y1, window=win, marks=marx,
                                  check=FALSE))
             })
    }

    if(!missing(j))
      x <- x[j] # invokes code above
    
    return(x)
 }
  


####################################################
# affine transformations
####################################################

affine.psp <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "psp")
  W <- affine.owin(X$window, mat=mat, vec=vec, ...)
  E <- X$ends
  ends0 <- affinexy(list(x=E$x0,y=E$y0), mat=mat, vec=vec)
  ends1 <- affinexy(list(x=E$x1,y=E$y1), mat=mat, vec=vec)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE),
      check=FALSE)
}

shift.psp <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "psp")
  if(!is.null(origin)) {
    stopifnot(is.character(origin))
    if(!missing(vec))
      warning("Argument vec ignored; argument origin has precedence.\n")
    origin <- pickoption("origin", origin, c(centroid="centroid",
                                             midpoint="midpoint",
                                             bottomleft="bottomleft"))
    W <- as.owin(X)
    locn <- switch(origin,
                   centroid={ unlist(centroid.owin(W)) },
                   midpoint={ c(mean(W$xrange), mean(W$yrange)) },
                   bottomleft={ c(W$xrange[1], W$yrange[1]) })
    return(shift(X, -locn))
  }
  # perform shift
  W <- shift.owin(X$window, vec=vec, ...)
  E <- X$ends
  ends0 <- shiftxy(list(x=E$x0,y=E$y0), vec=vec, ...)
  ends1 <- shiftxy(list(x=E$x1,y=E$y1), vec=vec, ...)
  Y <- psp(ends0$x, ends0$y, ends1$x, ends1$y,
           window=W, marks=marks(X, dfok=TRUE),
           check=FALSE)
  # tack on shift vector
  attr(Y, "lastshift") <- vec
  return(Y)
}

rotate.psp <- function(X, angle=pi/2, ..., centre=NULL) {
  verifyclass(X, "psp")
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  W <- rotate.owin(X$window, angle=angle, ...)
  E <- X$ends
  ends0 <- rotxy(list(x=E$x0,y=E$y0), angle=angle)
  ends1 <- rotxy(list(x=E$x1,y=E$y1), angle=angle)
  Y <- psp(ends0$x, ends0$y, ends1$x, ends1$y,
           window=W, marks=marks(X, dfok=TRUE),
           check=FALSE)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}

is.empty.psp <- function(x) { return(x$n == 0) } 

identify.psp <- function(x, ..., labels=seq_len(nsegments(x)), n=nsegments(x), plot=TRUE) {
  Y <- x
  W <- as.owin(Y)
  mids <- midpoints.psp(Y)
  if(!(is.numeric(n) && (length(n) == 1) && (n %% 1 == 0) && (n >= 0)))
    stop("n should be a single integer")
  out <- integer(0)
  while(length(out) < n) {
    xy <- locator(1)
    # check for interrupt exit
    if(length(xy$x) == 0)
      return(out)
    # find nearest segment
    X <- ppp(xy$x, xy$y, window=W)
    ident <- project2segment(X, Y)$mapXY
    # add to list
    if(ident %in% out) {
      cat(paste("Segment", ident, "already selected\n"))
    } else {
      if(plot) {
        # Display
        mi <- mids[ident]
        li <- labels[ident]
        text(mi$x, mi$y, labels=li)
      }
      out <- c(out, ident)
    }
  }
  # exit if max n reached
  return(out)
}

nsegments <- function(x) {
	UseMethod("nsegments")
}

nobjects.psp <- nsegments.psp <- function(x) {
   x$n
}

as.ppp.psp <- function (X, ..., fatal=TRUE) 
{
  Y <- endpoints.psp(X, which="both")
  m  <- marks(X)
  marks(Y) <- markappend(m, m)
  return(Y)
}

domain.psp <- Window.psp <- function(X, ...) { as.owin(X) }

"Window<-.psp" <- function(X, ..., value) {
  verifyclass(value, "owin")
  X[value]
}
