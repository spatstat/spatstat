#
# lpp.R
#
#  $Revision: 1.19 $   $Date: 2013/01/24 04:02:18 $
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
  
