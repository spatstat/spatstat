#
# lpp.R
#
#  $Revision: 1.39 $   $Date: 2015/11/25 02:56:11 $
#
# Class "lpp" of point patterns on linear networks

lpp <- function(X, L, ...) {
  stopifnot(inherits(L, "linnet"))
  localnames <- c("seg", "tp")
  if(is.matrix(X)) X <- as.data.frame(X)
  if(checkfields(X, c("x", "y", localnames))) {
    # includes spatial and local coordinates
    X <- as.data.frame(X)
    # local coords
    lo <- X[ , localnames, drop=FALSE]
    # spatial coords and marks
    df <- X[, !(names(X) %in% localnames), drop=FALSE]
    # validate local coordinates
    if(nrow(X) > 0) {
      nedge <- nsegments(L)
      if(with(lo, any(seg < 1 || seg > nedge)))
        stop("Segment index coordinate 'seg' exceeds bounds")
      if(with(lo, any(tp < 0 || tp > 1)))
        stop("Local coordinate 'tp' outside [0,1]")
    }
  } else {
    # local coordinates must be computed
    if(!is.ppp(X))
      X <- as.ppp(X, W=L$window, ...)
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
  splat("Point pattern on linear network")
  sd <- summary(x$data)
  np <- sd$ncases
  nama <- sd$col.names
  splat(np, ngettext(np, "point", "points"))
  ## check for unusual coordinates
  ctype <- x$ctype
  nam.m <- nama[ctype == "mark"]
  nam.t <- nama[ctype == "temporal"]
  nam.c <- setdiff(nama[ctype == "spatial"], c("x","y"))
  nam.l <- setdiff(nama[ctype == "local"], c("seg", "tp"))
  if(length(nam.c) > 0)
    splat("Additional spatial coordinates", commasep(sQuote(nam.c)))
  if(length(nam.l) > 0)
    splat("Additional local coordinates", commasep(sQuote(nam.l)))
  if(length(nam.t) > 0)
    splat("Additional temporal coordinates", commasep(sQuote(nam.t)))
  if((nmarks <- length(nam.m)) > 0) {
    if(nmarks > 1) {
      splat(nmarks, "columns of marks:", commasep(sQuote(nam.m)))
    } else {
      marx <- marks(x)
      if(is.factor(marx)) {
        exhibitStringList("Multitype, with possible types:", levels(marx))
      } else splat("Marks of type", sQuote(typeof(marx)))
    }
  }
  print(x$domain, ...)
  return(invisible(NULL))
}

plot.lpp <- function(x, ..., main, add=FALSE,
                     use.marks=TRUE, which.marks=NULL,
                     show.all=!add, show.window=FALSE,
                     do.plot=TRUE, multiplot=TRUE) {
  if(missing(main))
    main <- short.deparse(substitute(x))
  ## Handle multiple columns of marks as separate plots
  ##  (unless add=TRUE or which.marks selects a single column
  ##   or multipage = FALSE)
  if(use.marks && is.data.frame(mx <- marks(x))) {
    implied.all <- is.null(which.marks)
    want.several <- implied.all || is.data.frame(mx <- mx[,which.marks])
    do.several <- want.several && !add && multiplot
    if(do.several) {
      ## generate one plot for each column of marks
      y <- solapply(mx, setmarks, x=x)
      out <- do.call("plot",
                     c(list(x=y, main=main, do.plot=do.plot,
                            show.window=show.window),
                       list(...)))
      return(invisible(out))
    } 
    if(is.null(which.marks)) {
      which.marks <- 1
      if(do.plot) message("Plotting the first column of marks")
    }
  }
  ## determine space required, including legend
  P <- as.ppp(x)
  a <- plot(P, ..., do.plot=FALSE)
  if(!do.plot) return(a)
  ## initialise graphics space
  if(!add) {
    if(show.window) {
      plot(Window(P), main=main, invert=TRUE, ...)
    } else {
      b <- attr(a, "bbox")
      plot(b, type="n", main=main, ..., show.all=FALSE)
    }
  }
  ## plot linear network
  L <- as.linnet(x)
  do.call.matched("plot.linnet",
                  resolve.defaults(list(x=L, add=TRUE),
                                   list(...)),
                  extrargs=c("lty", "lwd", "col"))
  ## plot points, legend, title
  ans <- do.call.matched("plot.ppp",
                         c(list(x=P, add=TRUE, main=main,
                                show.all=show.all, show.window=FALSE),
                           list(...)),
                         extrargs=c("shape", "size", "pch", "cex",
                           "fg", "bg", "cols", "lty", "lwd", "etch",
                           "cex.main", "col.main", "line", "outer", "sub"))
  return(invisible(ans))
}


summary.lpp <- function(object, ...) {
  stopifnot(inherits(object, "lpp"))
  L <- object$domain
  result <- summary(L)
  np <- npoints(object)
  result$npoints <- np <- npoints(object)
  result$intensity <- np/result$totlength
  result$is.marked <- is.marked(object)
  result$is.multitype <- is.marked(object)
  if(result$is.marked) {
    mks <- marks(object)
    if(result$multiple.marks <- is.data.frame(mks)) {
      result$marknames <- names(mks)
      result$is.numeric <- FALSE
      result$marktype <- "dataframe"
      result$is.multitype <- FALSE
    } else {
      result$is.numeric <- is.numeric(mks)
      result$marknames <- "marks"
      result$marktype <- typeof(mks)
      result$is.multitype <- is.multitype(object)
    }
    if(result$is.multitype) {
      tm <- as.vector(table(mks))
      tfp <- data.frame(frequency=tm,
                        proportion=tm/sum(tm),
                        intensity=tm/result$totlength,
                        row.names=levels(mks))
      result$marks <- tfp
    } else 
      result$marks <- summary(mks)
  }
  class(result) <- "summary.lpp"
  return(result)
}

print.summary.lpp <- function(x, ...) {
  splat("Point pattern on linear network")
  splat(x$npoints, "points")
  splat("Linear network with",
        x$nvert, "vertices and",
        x$nline, "lines")
  u <- x$unitinfo
  dig <- getOption('digits')
  splat("Total length", signif(x$totlength, dig), u$plural, u$explain)
  splat("Average intensity", signif(x$intensity, dig),
        "points per", if(u$vanilla) "unit length" else u$singular)
  if(x$is.marked) {
    if(x$multiple.marks) {
      splat("Mark variables:", commasep(x$marknames, ", "))
      cat("Summary:\n")
      print(x$marks)
    } else if(x$is.multitype) {
      cat("Multitype:\n")
      print(signif(x$marks,dig))
    } else {
      splat("marks are ",
            if(x$is.numeric) "numeric, ",
            "of type ", sQuote(x$marktype),
            sep="")
      cat("Summary:\n")
      print(x$marks)
    }
  }
  print(x$win, prefix="Enclosing window: ")
  invisible(NULL)
}

intensity.lpp <- function(X, ...) {
  len <- sum(lengths.psp(as.psp(as.linnet(X))))
  if(is.multitype(X)) table(marks(X))/len else npoints(X)/len
}

is.lpp <- function(x) {
  inherits(x, "lpp")
}

is.multitype.lpp <- function(X, na.action="warn", ...) {
  marx <- marks(X)
  if(is.null(marx))
    return(FALSE)
  if((is.data.frame(marx) || is.hyperframe(marx)) && ncol(marx) > 1)
    return(FALSE)
  if(!is.factor(marx))
    return(FALSE)
  if((length(marx) > 0) && any(is.na(marx)))
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

as.lpp <- function(x, y=NULL, seg=NULL, tp=NULL, ...,
                   marks=NULL, L=NULL, check=FALSE, sparse) {
  nomore <- is.null(y) && is.null(seg) && is.null(tp) 
  if(inherits(x, "lpp") && nomore) {
    X <- x
    if(!missing(sparse) && !is.null(sparse))
      X$domain <- as.linnet(domain(X), sparse=sparse)
  } else {
    if(!inherits(L, "linnet"))
      stop("L should be a linear network")
    if(!missing(sparse) && !is.null(sparse))
      L <- as.linnet(L, sparse=sparse)
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

Window.lpp <- function(X, ...) { as.owin(X) }

as.owin.lpp <- function(W,  ..., fatal=TRUE) {
  as.owin(as.ppp(W, ..., fatal=fatal))
}

domain.lpp <- function(X, ...) { as.linnet(X) }

as.linnet.lpp <- function(X, ..., fatal=TRUE, sparse) {
  verifyclass(X, "lpp", fatal=fatal)
  L <- X$domain
  if(!missing(sparse))
    L <- as.linnet(L, sparse=sparse)
  return(L)
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

nsegments.lpp <- function(x) {
  return(x$domain$lines$n)
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

"[.lpp" <- function (x, i, j, drop=FALSE, ...) {
  if(!missing(i) && !is.null(i)) {
    if(is.owin(i)) {
      # spatial domain: call code for 'j'
      xi <- x[,i]
    } else {
      # usual row-type index
      da <- x$data
      daij <- da[i, , drop=FALSE]
      xi <- ppx(data=daij, domain=x$domain, coord.type=as.character(x$ctype))
      if(drop)
        xi <- xi[drop=TRUE] # call [.ppx to remove unused factor levels
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

shift.lpp <- function(X, vec=c(0,0), ..., origin=NULL) {
  verifyclass(X, "lpp")
  Y <- X
  Y$domain <- shift(X$domain, vec=vec, ..., origin=origin)
  vec <- getlastshift(Y$domain)
  Y$data[, c("x","y")] <- shiftxy(X$data[, c("x","y")], vec=vec)
  # tack on shift vector
  attr(Y, "lastshift") <- vec
  return(Y)
}

rotate.lpp <- function(X, angle=pi/2, ..., centre=NULL) {
  verifyclass(X, "lpp")
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  Y <- X
  Y$data[, c("x","y")] <- rotxy(X$data[, c("x","y")], angle=angle)
  Y$domain <- rotate(X$domain, angle=angle, ...)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}

rescale.lpp <- function(X, s, unitname) {
  if(missing(unitname)) unitname <- NULL
  if(missing(s)) s <- 1/unitname(X)$multiplier
  Y <- scalardilate(X, f=1/s)
  unitname(Y) <- rescale(unitname(X), s, unitname)
  return(Y)
}

superimpose.lpp <- function(..., L=NULL) {
  objects <- list(...)
  if(!is.null(L) && !inherits(L, "linnet"))
    stop("L should be a linear network")
  if(length(objects) == 0) {
    if(is.null(L)) return(NULL)
    emptyX <- lpp(list(x=numeric(0), y=numeric(0)), L)
    return(emptyX)
  }
  islpp <- unlist(lapply(objects, is.lpp))
  if(is.null(L) && !any(islpp))
    stop("Cannot determine linear network: no lpp objects given")
  nets <- unique(lapply(objects[islpp], as.linnet))
  if(length(nets) > 1)
    stop("Point patterns are defined on different linear networks")
  if(!is.null(L)) {
    nets <- unique(append(nets, list(L)))
    if(length(nets) > 1)
      stop("Argument L is a different linear network")
  }
  L <- nets[[1]]
  ## convert list(x,y) to linear network, etc
  if(any(!islpp))
    objects[!islpp] <- lapply(objects[!islpp], lpp, L=L)
  ## concatenate coordinates 
  locns <- do.call("rbind", lapply(objects, coords))
  ## concatenate marks (or use names of arguments)
  marx <- superimposeMarks(objects, sapply(objects, npoints))
  ## make combined pattern
  Y <- lpp(locns, L)
  marks(Y) <- marx
  return(Y)
}

#
# interactive plot for lpp objects
#

iplot.lpp <- function(x, ..., xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  stopifnot(is.lpp(x))
  ## predigest
  L <- domain(x)
  v <- vertices(L)
  deg <- vertexdegree(L)
  dv <- textstring(v, txt=paste(deg))
  y <- layered(lines=as.psp(L),
               vertices=v,
               degree=dv,
               points=as.ppp(x))
  iplot(y, ..., xname=xname, visible=c(TRUE, FALSE, FALSE, TRUE))
}

identify.lpp <- function(x, ...) {
  verifyclass(x, "lpp")
  P <- as.ppp(x)
  id <- identify(P$x, P$y, ...)
  if(!is.marked(x)) return(id)
  marks <- as.data.frame(P)[id, -(1:2)]
  out <- cbind(data.frame(id=id), marks)
  row.names(out) <- NULL
  return(out)
}
