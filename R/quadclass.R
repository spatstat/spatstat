#
#	quadclass.S
#
#	Class 'quad' to define quadrature schemes
#	in (rectangular) windows in two dimensions.
#
#	$Revision: 4.24 $	$Date: 2014/10/24 00:22:30 $
#
# An object of class 'quad' contains the following entries:
#
#	$data:	an object of class 'ppp'
#		defining the OBSERVATION window, 
#		giving the locations (& marks) of the data points.
#
#	$dummy:	object of class 'ppp'
#		defining the QUADRATURE window, 
#		giving the locations (& marks) of the dummy points.
#	
#	$w: 	vector giving the nonnegative weights for the
#		data and dummy points (data first, followed by dummy)
#
#		w may also have an attribute attr(w, "zeroes")
#               equivalent to (w == 0). If this is absent
#               then all points are known to have positive weights.
#
#       $param:
#               parameters that were used to compute the weights
#               and possibly to create the dummy points (see below).
#              
#       The combined (data+dummy) vectors of x, y coordinates of the points, 
#       and their weights, are extracted using standard functions 
#       x.quad(), y.quad(), w.quad() etc.
#
# ----------------------------------------------------------------------
#  Note about parameters:
#
#       If the quadrature scheme was created by quadscheme(),
#       then $param contains
#
#           $param$weight
#                list containing the values of all parameters
#                actually used to compute the weights.
#
#           $param$dummy
#                list containing the values of all parameters
#                actually used to construct the dummy pattern
#                via default.dummy();
#                or NULL if the dummy pattern was provided externally
#
#           $param$sourceid
#                vector mapping the quadrature points to the
#                original data and dummy points.
#
#   If you constructed the quadrature scheme manually, this
#   structure may not be present.
#
#-------------------------------------------------------------

quad <- function(data, dummy, w, param=NULL) {
  
  data <- as.ppp(data)
  dummy <- as.ppp(dummy)

  n <- data$n + dummy$n
	
  if(missing(w))
    w <- rep.int(1, n)
  else {
    w <- as.vector(w)
    if(length(w) != n)
      stop("length of weights vector w is not equal to total number of points")
  }

  if(is.null(attr(w, "zeroes")) && any( w == 0))
	attr(w, "zeroes") <- (w == 0)

  Q <- list(data=data, dummy=dummy, w=w, param=param)
  class(Q) <- "quad"

  invisible(Q)
}

# ------------------ extractor functions ----------------------

x.quad <- function(Q) {
  verifyclass(Q, "quad")
  c(Q$data$x, Q$dummy$x)
}

y.quad <- function(Q) {
  verifyclass(Q, "quad")
  c(Q$data$y, Q$dummy$y)
}

w.quad <- function(Q) {
  verifyclass(Q, "quad")
  Q$w
}

param.quad <- function(Q) {
  verifyclass(Q, "quad")
  Q$param
}
 
n.quad <- function(Q) {
  verifyclass(Q, "quad")
  Q$data$n + Q$dummy$n
}

marks.quad <- function(x, dfok=FALSE, ...) {
  verifyclass(x, "quad")
  dat <- x$data
  dum <- x$dummy
  if(dfok) warning("ignored dfok = TRUE; not implemented")
  mdat <- marks(dat, dfok=FALSE, ...)
  mdum <- marks(dum, dfok=FALSE, ...)
  if(is.null(mdat) && is.null(mdum))
    return(NULL)
  if(is.null(mdat))
    mdat <- rep.int(NA_integer_, dat$n)
  if(is.null(mdum))
    mdum <- rep.int(NA_integer_, dum$n)
  if(is.factor(mdat) && is.factor(mdum)) {
    mall <- cat.factor(mdat, mdum)
  } else mall <- c(mdat, mdum)
  return(mall)
}

is.marked.quad <- function(X, na.action="warn", ...) {
  marx <- marks(X, ...)
  if(is.null(marx))
    return(FALSE)
  if(any(is.na(marx)))
    switch(na.action,
           warn = {
             warning(paste("some mark values are NA in the point pattern",
                           short.deparse(substitute(X))))
           },
           fatal = {
             return(FALSE)
           },
           ignore = {}
           )
  return(TRUE)
}

is.multitype.quad <- function(X, na.action="warn", ...) {
  marx <- marks(X, ...)
  if(is.null(marx))
    return(FALSE)
  if(any(is.na(marx)))
    switch(na.action,
           warn = {
             warning(paste("some mark values are NA in the point pattern",
                           short.deparse(substitute(X))))
           },
           fatal = {
             return(FALSE)
           },
           ignore = {}
           )
  return(!is.data.frame(marx) && is.factor(marx))
}

is.data <- function(Q) {
  verifyclass(Q, "quad")
  return(c(rep.int(TRUE, Q$data$n),
	   rep.int(FALSE, Q$dummy$n)))
}

equals.quad <- function(Q) {
    # return matrix E such that E[i,j] = (X[i] == U[j])
    # where X = Q$data and U = union.quad(Q)
    n <- Q$data$n
    m <- Q$dummy$n
    E <- matrix(FALSE, nrow=n, ncol=n+m)
    diag(E) <- TRUE
    E
}

equalsfun.quad <- function(Q) {
  stopifnot(inherits(Q, "quad"))
  return(function(i,j) { i == j })
}

equalpairs.quad <- function(Q) {
  # return two-column matrix E such that
  #     X[E[i,1]] == U[E[i,2]] for all i
  # where X = Q$data and U = union.quad(Q)
  n <- Q$data$n
  return(matrix(rep.int(seq_len(n),2), ncol=2))
}
      
union.quad <- function(Q) {
  verifyclass(Q, "quad")
  ppp(x= c(Q$data$x, Q$dummy$x),
      y= c(Q$data$y, Q$dummy$y),
      window=Q$dummy$window,
      marks=marks.quad(Q),
      check=FALSE)
}
	
#
#   Plot a quadrature scheme
#
#
plot.quad <- function(x, ..., main, add=FALSE, dum=list(), tiles=FALSE) {
  if(missing(main) || is.null(main)) 
    main <- short.deparse(substitute(x))
  verifyclass(x, "quad")
  data <- x$data
  dummy <- x$dummy
  # determine plot parameters for dummy points
  dum <- resolve.defaults(dum, list(pch=".", add=TRUE))
  tt <- NULL
  if(tiles) {
    # show tiles that determined the weights
    wp <- x$param$weight
    tt <- NULL
    if(is.null(wp) || is.null(wp$method)) {
      warning("Tile information is not available")
    } else {
      switch(wp$method,
             grid = {
               ntile <- wp$ntile
               tt <- quadrats(as.owin(x), ntile[1], ntile[2])
             },
             dirichlet = {
               U <- union.quad(x)
               if(wp$exact) {
                 tt <- dirichlet(U)
               } else {
                 win <- as.mask(as.owin(U))
                 tileid <- image(exactdt(U)$i,
                                 win$xcol, win$yrow, win$xrange, win$yrange)
                 tt <- tess(image=tileid[win, drop=FALSE])
               }
             },
             warning("Unrecognised 'method' for tile weights")
             )
    }
  }
  pixeltiles <- !is.null(tt) && tt$type == "image"
  tileargs <- resolve.defaults(list(x=tt, main=main, add=add),
                               list(...),
                               if(!pixeltiles) list(col="grey") else NULL)
  if(!is.marked(data)) {
    if(!is.null(tt)) {
      do.call("plot", tileargs)
      add <- TRUE
    }
    plot(data, main=main, add=add, ...)
    do.call("plot", append(list(x=dummy), dum))
  } else if(is.multitype(data) && !add) {
    oldpar <- par(ask = interactive() &&
                  (.Device %in% c("X11", "GTK", "windows", "Macintosh")))
    on.exit(par(oldpar))
    data.marks <- marks(data)
    dummy.marks <- marks(dummy)
    types <- levels(data.marks)
    for(k in types) {
      add <- FALSE
      if(!is.null(tt)) {
        do.call("plot", tileargs)
        add <- TRUE
      }
      maink <- paste(main, "\n mark = ", k, sep="")
      plot(unmark(data[data.marks == k]), main=maink, add=add, ...)
      do.call("plot", append(list(x=unmark(dummy[dummy.marks == k])),
                             dum))
    }
  } else {
    if(!is.null(tt)) {
      do.call("plot", tileargs)
      add <- TRUE
    }
    plot(data, ..., main=main, add=add)
    do.call("plot", append(list(x=dummy), dum))
  }
  invisible(NULL)
}

# subset operator

"[.quad" <- function(x, ...) {
  U <- union.quad(x)
  Z <- is.data(x)
  w <- w.quad(x)
  # determine serial numbers of points to be included
  V <- U %mark% seq_len(U$n)
  i <- marks(V[...])
  # extract corresponding subsets of vectors
  Z <- Z[i]
  w <- w[i]
  # take subset of points, using any type of subset index
  U <- U[...]
  # stick together
  quad(U[Z], U[!Z], w)
}

domain.quad <- Window.quad <- function(X, ...) { as.owin(X) }

"Window<-.quad" <- function(X, ..., value) {
  verifyclass(value, "owin")
  return(X[value])
}

unitname.quad <- function(x) {
  return(unitname(x$data))
}

"unitname<-.quad" <- function(x, value) {
  unitname(x$data) <- value
  unitname(x$dummy) <- value
  return(x)
}


