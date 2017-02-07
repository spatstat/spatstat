#'
#'    polygood.R
#'
#'   Check validity of polygon data
#'
#'  $Revision: 1.1 $  $Date: 2017/01/08 00:37:54 $
#'

#' check validity of a polygonal owin

owinpolycheck <- function(W, verbose=TRUE) {
  verifyclass(W, "owin")
  stopifnot(W$type == "polygonal")

  # extract stuff
  B <- W$bdry
  npoly <- length(B)
  outerframe <- owin(W$xrange, W$yrange)
  # can't use as.rectangle here; we're still checking validity
  boxarea.mineps <- area.owin(outerframe) * (1 - 0.00001)

  # detect very large datasets
  BS <- object.size(B)
  blowbyblow <- verbose & (BS > 1e4 || npoly > 20)
  #
  
  answer <- TRUE
  notes <- character(0)
  err <- character(0)
  
  # check for duplicated points, self-intersection, outer frame
  if(blowbyblow) {
    cat(paste("Checking", npoly, ngettext(npoly, "polygon...", "polygons...")))
    pstate <- list()
  }

  dup <- self <- is.box <- logical(npoly)
  
  for(i in 1:npoly) {
    if(blowbyblow && npoly > 1L)
      pstate <- progressreport(i, npoly, state=pstate)
    Bi <- B[[i]]
    # check for duplicated vertices
    dup[i] <- as.logical(anyDuplicated(ppp(Bi$x, Bi$y,
                                           window=outerframe, check=FALSE)))
    if(dup[i] && blowbyblow)
      message(paste("Polygon", i, "contains duplicated vertices"))
    # check for self-intersection
    self[i] <- xypolyselfint(B[[i]], proper=TRUE, yesorno=TRUE)
    if(self[i] && blowbyblow)
      message(paste("Polygon", i, "is self-intersecting"))
    # check whether one of the current boundary polygons
    # is the bounding box itself (with + sign)
    is.box[i] <- (length(Bi$x) == 4) && (Area.xypolygon(Bi) >= boxarea.mineps)
  }
  if(blowbyblow)
    cat("done.\n")
  
  if((ndup <- sum(dup)) > 0) {
    whinge <- paste(ngettext(ndup, "Polygon", "Polygons"),
                    if(npoly == 1L) NULL else
                    commasep(which(dup)), 
                    ngettext(ndup, "contains", "contain"),
                    "duplicated vertices")
    notes <- c(notes, whinge)
    err <- c(err, "duplicated vertices")
    if(verbose) 
      message(whinge)
    answer <- FALSE
  }
  
  if((nself <- sum(self)) > 0) {
    whinge <-  paste(ngettext(nself, "Polygon", "Polygons"),
                     if(npoly == 1L) NULL else
                     commasep(which(self)),
                     ngettext(nself, "is", "are"),
                     "self-intersecting")
    notes <- c(notes, whinge)
    if(verbose) 
      message(whinge)
    err <- c(err, "self-intersection")
    answer <- FALSE
  }
  
  if(sum(is.box) > 1L) {
    answer <- FALSE
    whinge <- paste("Polygons",
                    commasep(which(is.box)),
                    "coincide with the outer frame")
    notes <- c(notes, whinge)
    err <- c(err, "polygons duplicating the outer frame")
  }
  
  # check for crossings between different polygons
  cross <- matrix(FALSE, npoly, npoly)
  if(npoly > 1L) {
    if(blowbyblow) {
      cat(paste("Checking for cross-intersection between",
                npoly, "polygons..."))
      pstate <- list()
    }
    P <- lapply(B, xypolygon2psp, w=outerframe, check=FALSE)
    for(i in seq_len(npoly-1L)) {
      if(blowbyblow)
        pstate <- progressreport(i, npoly-1L, state=pstate)
      Pi <- P[[i]]
      for(j in (i+1L):npoly) {
        crosses <- if(is.box[i] || is.box[j]) FALSE else {
          anycrossing.psp(Pi, P[[j]])
        }
        cross[i,j] <- cross[j,i] <- crosses
        if(crosses) {
          answer <- FALSE
          whinge <- paste("Polygons", i, "and", j, "cross over")
          notes <- c(notes, whinge)
          if(verbose) 
            message(whinge)
          err <- c(err, "overlaps between polygons")
        }
      }
    }
    if(blowbyblow)
      cat("done.\n")
  }

  err <- unique(err)
  attr(answer, "notes") <- notes
  attr(answer, "err") <-  err
  return(answer)
}

#' check for self-intersections in an xypolygon

xypolyselfint <- function(p, eps=.Machine$double.eps,
                          proper=FALSE, yesorno=FALSE, checkinternal=FALSE) {
  verify.xypolygon(p)
  n <- length(p$x)
  verbose <- (n > 1000)
  if(verbose)
    cat(paste("[Checking polygon with", n, "edges..."))
  x0 <- p$x
  y0 <- p$y
  dx <- diff(x0[c(1:n,1L)])
  dy <- diff(y0[c(1:n,1L)])
  if(yesorno) {
    # get a yes-or-no answer
    answer <- .C("xypsi",
                 n=as.integer(n),
                 x0=as.double(x0),
                 y0=as.double(y0),
                 dx=as.double(dx),
                 dy=as.double(dy),
                 xsep=as.double(2 * max(abs(dx))),
                 ysep=as.double(2 * max(abs(dy))),
                 eps=as.double(eps),
                 proper=as.integer(proper),
                 answer=as.integer(integer(1L)))$answer
    if(verbose)
      cat("]\n")
    return(answer != 0)
  }
  out <- .C("Cxypolyselfint",
            n=as.integer(n),
            x0=as.double(x0),
            y0=as.double(y0),
            dx=as.double(dx),
            dy=as.double(dy), 
            eps=as.double(eps),
            xx=as.double(numeric(n^2)),
            yy=as.double(numeric(n^2)),
            ti=as.double(numeric(n^2)),
            tj=as.double(numeric(n^2)),
            ok=as.integer(integer(n^2)))

  uhoh <- (matrix(out$ok, n, n) != 0)
  if(proper) {
    # ignore cases where two vertices coincide 
    ti <- matrix(out$ti, n, n)[uhoh]
    tj <- matrix(out$tj, n, n)[uhoh]
    i.is.vertex <- (abs(ti) < eps) | (abs(ti - 1) < eps)
    j.is.vertex <- (abs(tj) < eps) | (abs(tj - 1) < eps)
    dup <- i.is.vertex & j.is.vertex
    uhoh[uhoh] <- !dup
  }
  if(checkinternal && any(uhoh != t(uhoh)))
    warning("Internal error: incidence matrix is not symmetric")
  xx <- matrix(out$xx, n, n)
  yy <- matrix(out$yy, n, n)
  uptri <- (row(uhoh) < col(uhoh))
  xx <- as.vector(xx[uhoh & uptri])
  yy <- as.vector(yy[uhoh & uptri])
  result <- list(x=xx, y=yy)
  if(verbose)
    cat("]\n")
  return(result)
}
  
