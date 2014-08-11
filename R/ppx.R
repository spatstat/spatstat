#
#   ppx.R
#
#  class of general point patterns in any dimension
#
#  $Revision: 1.39 $  $Date: 2013/01/15 09:22:59 $
#

ppx <- local({
  
  ctype.table <- c("spatial", "temporal", "local", "mark")
  ctype.real  <- c(TRUE,      TRUE,       FALSE,   FALSE)

  ppx <- function(data, domain=NULL, coord.type=NULL) {
    data <- as.hyperframe(data)
    # columns suitable for spatial coordinates
    suitable <- with(unclass(data),
                     vtype == "dfcolumn" &
                     (vclass == "numeric" | vclass == "integer"))
    if(is.null(coord.type)) {
      # assume all suitable columns of data are spatial coordinates
      # and all other columns are marks.
      ctype <- ifelse(suitable, "spatial", "mark")
    } else {
      stopifnot(is.character(coord.type))
      stopifnot(length(coord.type) == ncol(data))
      ctypeid <- pmatch(coord.type, ctype.table, duplicates.ok=TRUE)
      # validate
      if(any(uhoh <- is.na(ctypeid)))
        stop(paste("Unrecognised coordinate",
                   ngettext(sum(uhoh), "type", "types"),
                   commasep(sQuote(coord.type[uhoh]))))
      if(any(uhoh <- (!suitable & ctype.real[ctypeid]))) {
        nuh <- sum(uhoh)
        stop(paste(ngettext(nuh, "Coordinate", "Coordinates"),
                   commasep(sQuote(names(data)[uhoh])),
                   ngettext(nuh, "does not", "do not"),
                   "contain real numbers"))
      }
      ctype <- ctype.table[ctypeid]
    }
    ctype <- factor(ctype, levels=ctype.table)
    out <- list(data=data, ctype=ctype, domain=domain)
    class(out) <- "ppx"
    return(out)
  }

  ppx
})


is.ppx <- function(x) { inherits(x, "ppx") }

nobjects.ppx <- npoints.ppx <- function(x) { nrow(x$data) }

print.ppx <- function(x, ...) {
  cat("Multidimensional point pattern\n")
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
    cat(paste(sum(ismark), ngettext(sum(ismark), "column", "columns"),
              "of local coordinates:",
              commasep(sQuote(nama[ismark])), "\n"))
  if(any(ismark <- (x$ctype == "mark"))) 
    cat(paste(sum(ismark), ngettext(sum(ismark), "column", "columns"),
              "of marks:",
              commasep(sQuote(nama[ismark])), "\n"))
  if(!is.null(x$domain)) {
    cat("Domain:\n\t")
    print(x$domain)
  }
  invisible(NULL)
}

summary.ppx <- function(object, ...) { print(object, ...) }

plot.ppx <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  coo <- coords(x, local=FALSE)
  dom <- x$domain
  m <- ncol(coo)
  if(m == 1) {
    coo <- coo[,1]
    ran <- diff(range(coo))
    ylim <- c(-1,1) * ran/20
    do.call("plot.default",
            resolve.defaults(list(coo, rep(0, length(coo))),
                             list(...),
                             list(asp=1, ylim=ylim,
                                  axes=FALSE, xlab="", ylab="")))
    axis(1, pos=ylim[1])
  } else if(m == 2) {
    if(is.null(dom)) {
      # plot x, y coordinates only
      nama <- names(coo)
      do.call.matched("plot.default",
                      resolve.defaults(list(x=coo[,1], y=coo[,2], asp=1),
                                       list(...),
                                       list(main=xname),
                                       list(xlab=nama[1], ylab=nama[2])))
    } else {
      add <- resolve.defaults(list(...), list(add=FALSE))$add
      if(!add) {
        # plot domain, whatever it is
        do.call("plot", resolve.defaults(list(dom),
                                       list(...),
                                       list(main=xname)))
      }
      # convert to ppp
      x2 <- ppp(coo[,1], coo[,2], window=as.owin(dom),
                marks=as.data.frame(marks(x)), check=FALSE)
      # invoke plot.ppp
      return(do.call("plot", resolve.defaults(list(x2),
                                              list(add=TRUE),
                                              list(...))))
    }
  } else if(m == 3) {
    # convert to pp3
    if(is.null(dom))
      dom <- box3(range(coo[,1]), range(coo[,2]), range(coo[,3]))
    x3 <- pp3(coo[,1], coo[,2], coo[,3], dom)
    # invoke plot.pp3
    nama <- names(coo)
    do.call("plot",
            resolve.defaults(list(x3),
                             list(...),
                             list(main=xname),
                             list(xlab=nama[1], ylab=nama[2], zlab=nama[3])))
  } else stop(paste("Don't know how to plot a general point pattern in",
               ncol(coo), "dimensions"))
  return(invisible(NULL))
}

"[.ppx" <- function (x, i, ...) {
  da <- x$data
  daij <- da[i, , drop=FALSE]
  out <- list(data=daij, ctype=x$ctype, domain=x$domain)
  class(out) <- "ppx"
  return(out)
}

coords <- function(x, ...) {
  UseMethod("coords")
}

coords.ppx <- function(x, ..., spatial=TRUE, temporal=TRUE, local=TRUE) {
  ctype <- x$ctype
  chosen <- (ctype == "spatial" & spatial) |
            (ctype == "temporal" & temporal) | 
            (ctype == "local" & local) 
  as.data.frame(x$data[, chosen])
}

coords.ppp <- function(x, ...) { data.frame(x=x$x,y=x$y) }

"coords<-" <- function(x, ..., value) {
  UseMethod("coords<-")
}

"coords<-.ppp" <- function(x, ..., value) {
  win <- x$window
  if(is.null(value)) {
    # empty pattern
    return(ppp(window=win))
  }
  value <- as.data.frame(value)
  if(ncol(value) != 2)
    stop("Expecting a 2-column matrix or data frame, or two vectors")
  result <- as.ppp(value, win)
  marks(result) <- marks(x)
  return(result)
}

"coords<-.ppx" <- function(x, ..., spatial=TRUE, temporal=TRUE, local=TRUE, value) {
  ctype <- x$ctype
  chosen <- (ctype == "spatial" & spatial) |
            (ctype == "temporal" & temporal) | 
            (ctype == "local" & local) 
  x$data[, chosen] <- value
  return(x)
}

as.hyperframe.ppx <- function(x, ...) { x$data }

as.data.frame.ppx <- function(x, ...) { as.data.frame(x$data, ...) } 

as.matrix.ppx <- function(x, ...) { as.matrix(as.data.frame(x, ...)) }

marks.ppx <- function(x, ..., drop=TRUE) {
  ctype <- x$ctype
  chosen <- (ctype == "mark")
  if(!any(chosen)) return(NULL)
  x$data[, chosen, drop=drop]
}

"marks<-.ppx" <- function(x, ..., value) {
  ctype <- x$ctype
  retain <- (ctype != "mark")
  coorddata <- x$data[, retain, drop=TRUE]
  if(is.null(value)) {
    newdata <- coorddata
    newctype <- ctype[retain]
  } else {
    if(is.matrix(value) && nrow(value) == nrow(x$data)) {
      # assume matrix is to be treated as data frame
      value <- as.data.frame(value)
    }
    if(!is.data.frame(value) && !is.hyperframe(value)) 
      value <- hyperframe(marks=value)
    if(is.hyperframe(value) || is.hyperframe(coorddata)) {
      value <- as.hyperframe(value)
      coorddata <- as.hyperframe(coorddata)
    }
    if(ncol(value) == 0) {
      newdata <- coorddata
      newctype <- ctype[retain]
    } else {
      newdata <- cbind(coorddata, value)
      newctype <- factor(c(as.character(ctype[retain]),
                           rep("mark", ncol(value))),
                         levels=levels(ctype))
    }
  }
  out <- list(data=newdata, ctype=newctype, domain=x$domain)
  class(out) <- "ppx"
  return(out)
}

unmark.ppx <- function(X) {
  marks(X) <- NULL
  return(X)
}

markformat.ppx <- function(x) {
  mf <- x$markformat
  if(is.null(mf)) 
    mf <- markformat(marks(x))
  return(mf)
}

boxx <- function(..., unitname=NULL) {
  if(length(list(...)) == 0)
    stop("No data")
  ranges <- data.frame(...)
  nama <- names(list(...))
  if(is.null(nama) || !all(nzchar(nama)))
    names(ranges) <- paste("x", 1:ncol(ranges),sep="")
  if(nrow(ranges) != 2)
    stop("Data should be vectors of length 2")
  if(any(unlist(lapply(ranges, diff)) <= 0))
    stop("Illegal range: Second element <= first element")
  out <- list(ranges=ranges, units=as.units(unitname))
  class(out) <- "boxx"
  return(out)
}

print.boxx <- function(x, ...) {
  m <- ncol(x$ranges)
  cat(paste(m, "-dimensional box:\n", sep=""))
  bracket <- function(z) paste("[",
                               paste(signif(z, 5), collapse=", "),
                               "]", sep="")
  v <- paste(unlist(lapply(x$ranges, bracket)), collapse=" x ")
  s <- summary(unitname(x))
  cat(paste(v, s$plural, s$explain, "\n"))
  invisible(NULL)
}

unitname.boxx <- function(x) { x$units }

"unitname<-.boxx" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

unitname.ppx <- function(x) { unitname(x$domain) }

"unitname<-.ppx" <- function(x, value) {
  d <- x$domain
  unitname(d) <- value
  x$domain <- d
  return(x)
}

sidelengths.boxx <- function(x) {
  stopifnot(inherits(x, "boxx"))
  y <- unlist(lapply(x$ranges, diff))
  return(y)
}
  
volume.boxx <- function(x) {
  prod(sidelengths(x))
}

diameter.boxx <- function(x) {
  d <- sqrt(sum(sidelengths(x)^2))
  return(d)
}

shortside.boxx <- function(x) {
  return(min(sidelengths(x)))
}

eroded.volumes.boxx <- function(x, r) {
  len <- sidelengths(x)
  ero <- sapply(as.list(len), function(z, r) { pmax(0, z - 2 * r)}, r=r)
  apply(ero, 1, prod)
}

runifpointx <- function(n, domain) {
  stopifnot(inherits(domain, "boxx"))
  coo <- lapply(domain$ranges,
                function(ra, n) { runif(n, min=ra[1], max=ra[2]) },
                n=n)
  df <- do.call("data.frame", coo)
  ppx(df, domain)
}

rpoisppx <- function(lambda, domain) {
  stopifnot(inherits(domain, "boxx"))
  vol <- volume.boxx(domain)
  stopifnot(is.numeric(lambda) && length(lambda) == 1 && lambda >= 0)
  n <- rpois(1, lambda * vol)
  runifpointx(n, domain)
}

