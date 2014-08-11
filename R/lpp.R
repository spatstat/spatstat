#
# lpp.R
#
#  $Revision: 1.22 $   $Date: 2013/10/20 00:51:27 $
#
# Class "lpp" of point patterns on linear networks

lpp <- function(X, L) {
  stopifnot(inherits(L, "linnet"))
  localnames <- c("seg", "tp")
  if(checkfields(X, c("x", "y", localnames))) {
    # includes spatial and local coordinates
    X <- as.data.frame(X)
    # local coords
    lo <- X[ , localnames, drop=FALSE]
    # spatial coords and marks
    df <- X[, !(names(X) %in% localnames), drop=FALSE]
    # validate local coordinates
    if(nrow(X) > 0) {
      nedge <- nobjects(as.psp(L))
      if(with(lo, any(seg < 1 || seg > nedge)))
        stop("Segment index coordinate 'seg' exceeds bounds")
      if(with(lo, any(tp < 0 || tp > 1)))
        stop("Local coordinate 'tp' outside [0,1]")
    }
  } else {
    # local coordinates must be computed
    if(!is.ppp(X))
      X <- as.ppp(X, W=L$window)
    # project to segment
    pro <- project2segment(X, as.psp(L))
    # projected points (spatial coordinates and marks)
    df  <- as.data.frame(pro$Xproj)
    # local coordinates
    lo  <- data.frame(seg=pro$mapXY, tp=pro$tp)
  }
  # combine spatial, local, marks
  nmark <- ncol(df) - 2
  if(nmark == 0) {
    df <- cbind(df, lo)
    ctype <- c(rep("s", 2), rep("l", 2))
  } else {
    df <- cbind(df[,1:2], lo, df[, -(1:2), drop=FALSE])
    ctype <- c(rep("s", 2), rep("l", 2), rep("m", nmark))
  }
  out <- ppx(data=df, domain=L, coord.type=ctype)
  class(out) <- c("lpp", class(out))
  return(out)
}

print.lpp <- function(x, ...) {
  stopifnot(inherits(x, "lpp"))
  cat("Point pattern on linear network\n")
  sd <- summary(x$data)
  np <- sd$ncases
  nama <- sd$col.names
  cat(paste(np, ngettext(np, "point", "points"), "\n"))
  if(any(iscoord <- (x$ctype == "spatial")))
    cat(paste(sum(iscoord), "-dimensional space coordinates ",
              paren(paste(nama[iscoord], collapse=",")), "\n", sep=""))
  if(any(istime <- (x$ctype == "temporal")))
    cat(paste(sum(istime), "-dimensional time coordinates ",
              paren(paste(nama[istime], collapse=",")), "\n", sep=""))
  if(any(islocal <- (x$ctype == "local"))) 
    cat(paste(sum(islocal), ngettext(sum(islocal), "column", "columns"),
              "of local coordinates:",
              commasep(sQuote(nama[islocal])), "\n"))
  if(any(ismark <- (x$ctype == "mark"))) 
    cat(paste(sum(ismark), ngettext(sum(ismark), "column", "columns"),
              "of marks:",
              commasep(sQuote(nama[ismark])), "\n"))
  print(x$domain, ...)
  invisible(NULL)
}

# plot.lpp removed: plot.ppx sufficient

summary.lpp <- function(object, ...) {
  stopifnot(inherits(object, "lpp"))
  L <- object$domain
  npoints <- nrow(object$data)
  totlen <-  sum(lengths.psp(L$lines))
  marx <- marks(object)
  summarx <- if(is.null(marx)) NULL else summary(marx)
  out <- list(npoints=npoints,
              totlength=totlen,
              intensity=npoints/totlen,
              nvert=L$vertices$n,
              nedge=L$lines$n,
              unitinfo=summary(unitname(L)),
              marks=summarx)
  class(out) <- "summary.lpp"
  return(out)
}

print.summary.lpp <- function(x, ...) {
  cat("Point pattern on linear network\n")
  cat(paste(x$npoints, "points\n"))
  cat(paste("Linear network with",
            x$nvert, "vertices and",
            x$nedge, "edges\n"))
  u <- x$unitinfo
  cat(paste("Total edge length", x$totlength, u$plural, u$explain, "\n"))
  cat(paste("Average intensity", x$intensity,
            "points per", if(u$vanilla) "unit length" else u$singular, "\n"))
  if(!is.null(x$marks)) {
    cat("Marks:\n")
    print(x$marks)
  }
  invisible(NULL)
}

intensity.lpp <- function(X, ...) {
  len <- sum(lengths.psp(as.psp(as.linnet(X))))
  if(is.multitype(X)) table(marks(X))/len else npoints(X)/len
}

is.lpp <- function(x) {
  inherits(x, "lpp")
}

as.lpp <- function(x, y=NULL, seg=NULL, tp=NULL, ...,
                   marks=NULL, L=NULL, check=FALSE) {
  nomore <- is.null(y) && is.null(seg) && is.null(tp)
  if(inherits(x, "lpp") && nomore) {
    X <- x
  } else {
    if(!inherits(L, "linnet"))
      stop("L should be a linear network")
    if(is.ppp(x) && nomore) {
      X <- lpp(x, L)
    } else {
      xy <- xy.coords(x,y)[c("x", "y")]
      if(!is.null(seg) && !is.null(tp)) {
        # add segment map information
        xy <- append(xy, list(seg=seg, tp=tp))
      } else {
        # convert to ppp, typically suppressing check mechanism
        xy <- as.ppp(xy, W=as.owin(L), check=check)
      }
      X <- lpp(xy, L)
    }
  }
  if(!is.null(marks))
    marks(X) <- marks
  return(X)
}

as.ppp.lpp <- function(X, ..., fatal=TRUE) {
  verifyclass(X, "lpp", fatal=fatal)
  L <- X$domain
  Y <- as.ppp(coords(X, temporal=FALSE, local=FALSE),
              W=L$window, check=FALSE)
  marks(Y) <- marks(X)
  return(Y)
}

as.owin.lpp <- function(W,  ..., fatal=TRUE) {
  as.owin(as.ppp(W, ..., fatal=fatal))
}

as.linnet.lpp <- function(X, ..., fatal=TRUE) {
  verifyclass(X, "lpp", fatal=fatal)
  X$domain
}
  
"[.lpp" <- function (x, i, ...) {
  # invoke [.ppx
  y <- NextMethod("[")
  class(y) <- c("lpp", class(y))
  return(y)
}

unitname.lpp <- function(x) {
  u <- unitname(x$domain)
  return(u)
}

"unitname<-.lpp" <- function(x, value) {
  w <- x$domain
  unitname(w) <- value
  x$domain <- w
  return(x)
}

"marks<-.lpp" <- function(x, ..., value) {
  Y <- NextMethod("marks<-")
  class(Y) <- c("lpp", class(Y))
  Y
}
  
unmark.lpp <- function(X) {
  Y <- NextMethod("unmark")
  class(Y) <- c("lpp", class(Y))
  Y
}

as.psp.lpp <- function(x, ..., fatal=TRUE){
  verifyclass(x, "lpp", fatal=fatal)
  return(x$domain$lines)
}

local2lpp <- function(L, seg, tp, X=NULL) {
  stopifnot(inherits(L, "linnet"))
  if(is.null(X)) {
    # map to (x,y)
    Ldf <- as.data.frame(L$lines)
    dx <- with(Ldf, x1-x0)
    dy <- with(Ldf, y1-y0)
    x <- with(Ldf, x0[seg] + tp * dx[seg])
    y <- with(Ldf, y0[seg] + tp * dy[seg])
  } else {
    x <- X$x
    y <- X$y
  }
  # compile into data frame
  data <- data.frame(x=x, y=y, seg=seg, tp=tp)
  ctype <- c("s", "s", "l", "l")
  out <- ppx(data=data, domain=L, coord.type=ctype)
  class(out) <- c("lpp", class(out))
  return(out)
}

####################################################
# subset extractor
####################################################

"[.lpp" <- function (x, i, j, ...) {
  if(!missing(i) && !is.null(i)) {
    if(is.owin(i)) {
      # spatial domain: call code for 'j'
      xi <- x[,i]
    } else {
      # usual row-type index
      da <- x$data
      daij <- da[i, , drop=FALSE]
      xi <- ppx(data=daij, domain=x$domain, coord.type=as.character(x$ctype))
      class(xi) <- c("lpp", class(xi))
    }
    x <- xi
  } 
  if(missing(j) || is.null(j))
    return(x)
  stopifnot(is.owin(j))
  W <- j
  L <- x$domain
  da <- x$data
  # Find vertices that lie inside 'j'
  okvert <- inside.owin(L$vertices, w=W)
  # find segments whose endpoints both lie in 'upper'
  okedge <- okvert[L$from] & okvert[L$to]
  # assign new serial numbers to vertices, and recode 
  newserial <- cumsum(okvert)
  newfrom <- newserial[L$from[okedge]]
  newto   <- newserial[L$to[okedge]]
  # make new linear network
  Lnew <- linnet(L$vertices[W], edges=cbind(newfrom, newto))
  # find data points that lie on accepted segments
  coo <- coords(x)
  okxy <- okedge[coo$seg]
  cook <- coo[okxy,]
  # make new lpp object
  dfnew <- data.frame(x=cook$x,
                      y=cook$y,
                      seg=cook$seg,
                      tp=cook$tp)
  ctype <- c(rep("spatial", 2), rep("local", 2))
  xj <- ppx(data=dfnew, domain=Lnew, coord.type=ctype)
  class(xj) <- c("lpp", class(xj))
  marks(xj) <- marks(x[okxy])
  return(xj)
}

####################################################
# affine transformations
####################################################

scalardilate.lpp <- function(X, f, ...) {
  trap.extra.arguments(..., .Context="In scalardilate(X,f)")
  check.1.real(f, "In scalardilate(X,f)")
  stopifnot(is.finite(f) && f > 0)
  Y <- X
  Y$data$x <- f * X$data$x
  Y$data$y <- f * X$data$y
  Y$domain <- scalardilate(X$domain, f)
  return(Y)
}

affine.lpp <- function(X,  mat=diag(c(1,1)), vec=c(0,0), ...) {
  verifyclass(X, "lpp")
  Y <- X
  Y$data[, c("x","y")] <- affinexy(X$data[, c("x","y")], mat=mat, vec=vec)
  Y$domain <- affine(X$domain, mat=mat, vec=vec, ...)
  return(Y)
}

shift.lpp <- function(X, ...) {
  verifyclass(X, "lpp")
  Y <- X
  Y$domain <- shift(X$domain, ...)
  vec <- attr(Y$domain, "lastshift")
  Y$data[, c("x","y")] <- shiftxy(X$data[, c("x","y")], vec=vec)
  # tack on shift vector
  attr(Y, "lastshift") <- vec
  return(Y)
}

rotate.lpp <- function(X, angle=pi/2, ...) {
  verifyclass(X, "lpp")
  Y <- X
  Y <- X
  Y$data[, c("x","y")] <- rotxy(X$data[, c("x","y")], angle=angle)
  Y$domain <- rotate(X$domain, angle=angle, ...)
  return(Y)
}

rescale.lpp <- function(X, s) {
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- scalardilate(X, f=1/s)
  unitname(Y) <- rescale(unitname(X), s)
  return(Y)
}
