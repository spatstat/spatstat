#
#    xypolygon.S
#
#    $Revision: 1.58 $    $Date: 2013/05/01 08:08:39 $
#
#    low-level functions defined for polygons in list(x,y) format
#
verify.xypolygon <- function(p, fatal=TRUE) {
  whinge <- NULL
  if(!is.list(p) || !all(c("x","y") %in% names(p)))
    whinge <- "polygon must be a list with components x and y"
  else if(is.null(p$x) || is.null(p$y) || !is.numeric(p$x) || !is.numeric(p$y))
    whinge <- "components x and y must be numeric vectors"
  else if(length(p$x) != length(p$y))
    whinge <- "lengths of x and y vectors unequal"
  ok <- is.null(whinge)
  if(!ok && fatal)
    stop(whinge)
  return(ok)
}

inside.xypolygon <- function(pts, polly, test01=TRUE, method="C") {
  # pts:  list(x,y) points to be tested
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  # test01: logical - if TRUE, test whether all values in output are 0 or 1
  pts <- xy.coords(pts, NULL)
  verify.xypolygon(polly)
  
  x <- pts$x
  y <- pts$y
  xp <- polly$x
  yp <- polly$y

  full.npts <- npts <- length(x)
  nedges <- length(xp)   # sic

  # Check for points (x,y) that coincide with vertices (xp, yp)
  # Handle them separately
  DUP <- spatstat.options("dupC")
  z <- .C("Cmatchxy",
          na=as.integer(npts),
          xa=as.double(x),
          ya=as.double(y),
          nb=as.integer(nedges),
          xb=as.double(xp),
          yb=as.double(yp),
          match=as.integer(integer(npts)),
          DUP=DUP)
#          PACKAGE="spatstat")
  is.vertex <- (z$match != 0)
  retain <- !is.vertex
  # Remove vertices from subsequent consideration; replace them later
  if(vertices.present <- !all(retain)) {
    x <- x[retain]
    y <- y[retain]
    npts <- sum(retain)
  }

  #-------------   MAIN ALGORITHM -------------------------------
  score <- numeric(npts)
  on.boundary <- rep.int(FALSE, npts)

  if(anyretain<- any(retain)) {
    switch(method,
           C={
             #------------------ call C routine ------------------
             temp <- .C("inxyp",
                        x=as.double(x),
                        y=as.double(y),
                        xp=as.double(xp),
                        yp=as.double(yp),
                        npts=as.integer(npts),
                        nedges=as.integer(nedges),
                        score=as.integer(score),
                        onbndry=as.integer(on.boundary))
#                        PACKAGE="spatstat")
             score <- temp$score/2
             on.boundary <- as.logical(temp$onbndry)
           },
           Fortran={
             #------------------ call Fortran routine ------------------
             temp <- .Fortran("inxyp",
                              x=as.double(x),
                              y=as.double(y),
                              xp=as.double(xp),
                              yp=as.double(yp),
                              npts=as.integer(npts),
                              nedges=as.integer(nedges),
                              score=as.double(score),
                              onbndry=as.logical(on.boundary))
#                              PACKAGE="spatstat")
             score <- temp$score
             on.boundary <- temp$onbndry
           },
           interpreted={
             #----------------- original interpreted code --------------
             for(i in 1:nedges) {
               x0 <- xp[i]
               y0 <- yp[i]
               x1 <- if(i == nedges) xp[1] else xp[i+1]
               y1 <- if(i == nedges) yp[1] else yp[i+1]
               dx <- x1 - x0
               dy <- y1 - y0
               if(dx < 0) {
                 # upper edge
                 xcriterion <- (x - x0) * (x - x1)
                 consider <- (xcriterion <= 0)
                 if(any(consider)) {
                   ycriterion <-
                     y[consider] * dx - x[consider] * dy +  x0 * dy - y0 * dx
                   # closed inequality
                   contrib <- (ycriterion >= 0) *
                     ifelseAB(xcriterion[consider] == 0, 1/2, 1)
                   # positive edge sign
                   score[consider] <- score[consider] + contrib
                   # detect whether any point lies on this segment
                   on.boundary[consider] <-
                     on.boundary[consider] | (ycriterion == 0)
                 }
               } else if(dx > 0) {
                 # lower edge
                 xcriterion <- (x - x0) * (x - x1)
                 consider <- (xcriterion <= 0)
                 if(any(consider)) {
                   ycriterion <-
                     y[consider] * dx - x[consider] * dy + x0 * dy - y0 * dx
                   # open inequality
                   contrib <- (ycriterion < 0) *
                     ifelseAB(xcriterion[consider] == 0, 1/2, 1)
                   # negative edge sign
                   score[consider] <- score[consider] - contrib
                   # detect whether any point lies on this segment
                   on.boundary[consider] <-
                     on.boundary[consider] | (ycriterion == 0)
                 }
               } else {
                 # vertical edge
                 consider <- (x == x0)
                 if(any(consider)) {
                   # zero score
                   # detect whether any point lies on this segment
                   yconsider <- y[consider]
                   ycriterion <- (yconsider - y0) * (yconsider - y1)
                   on.boundary[consider] <-
                     on.boundary[consider] | (ycriterion <= 0)
                 }
               }
             }
           },
           stop(paste("Unrecognised choice for", sQuote("method")))
           )
  }
  
  #------------------- END SWITCH ------------------------------

  # replace any polygon vertices that were temporarily removed
  if(vertices.present) {
    full.score <- numeric(full.npts)
    full.on.boundary <- rep.int(FALSE, full.npts)
    if(anyretain) {
      full.score[retain] <- score
      full.on.boundary[retain] <- on.boundary
    }
    full.score[is.vertex] <- 1
    full.on.boundary[is.vertex] <- TRUE
    score       <- full.score
    on.boundary <- full.on.boundary
    npts        <- full.npts
    
  }
    
  #-------------------------------------------------
  
  # any point recognised as lying on the boundary gets score 1.
  score[on.boundary] <- 1

  if(test01) {
    # check sanity
    if(!all((score == 0) | (score == 1)))
      warning("internal error: some scores are not equal to 0 or 1")
  }

  attr(score, "on.boundary") <- on.boundary
  
  return(score)
}

owinpoly2mask <- function(w, rasta, check=TRUE) {
  if(check) {
    verifyclass(w, "owin")
    stopifnot(w$type == "polygonal")
    verifyclass(rasta, "owin")
    stopifnot(rasta$type == "mask")
  }
  
  bdry <- w$bdry

  x0    <- rasta$xcol[1]
  y0    <- rasta$yrow[1]
  xstep <- rasta$xstep
  ystep <- rasta$ystep
  dimyx <- rasta$dim
  nx    <- dimyx[2]
  ny    <- dimyx[1]

  epsilon <- with(.Machine, double.base^floor(double.ulp.digits/2))

  score <- numeric(nx*ny)
  
  for(i in seq_along(bdry)) {
    p <- bdry[[i]]
    xp <- p$x
    yp <- p$y
    np <- length(p$x)
    # repeat last vertex
    xp <- c(xp, xp[1])
    yp <- c(yp, yp[1])
    np <- np + 1
    # rescale coordinates so that pixels are at integer locations
    xp <- (xp - x0)/xstep
    yp <- (yp - y0)/ystep
    # avoid exact integer locations for vertices
    whole <- (ceiling(xp) == floor(xp))
    xp[whole] <-  xp[whole] + epsilon
    whole <- (ceiling(yp) == floor(yp))
    yp[whole] <-  yp[whole] + epsilon
    # call C
    z <- .C("poly2imI",
            xp=as.double(xp),
            yp=as.double(yp),
            np=as.integer(np),
            nx=as.integer(nx),
            ny=as.integer(ny),
            out=as.integer(integer(nx * ny)))
#            PACKAGE="spatstat")
    if(i == 1)
      score <- z$out
    else 
      score <- score + z$out
  }
  status <- (score != 0)
  out <- owin(rasta$xrange, rasta$yrange, mask=matrix(status, ny, nx))
  return(out)
}

is.hole.xypolygon <- function(polly) {
  h <- polly$hole
  if(is.null(h))
    h <- (area.xypolygon(polly) < 0)
  return(h)
}
  
area.xypolygon <- function(polly) {
  #
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  #

  # area could be pre-calculated
  if(!is.null(pa <- polly$area) && is.numeric(pa) && length(pa)==1)
    return(pa)

  # else calculate
  verify.xypolygon(polly)
  xp <- polly$x
  yp <- polly$y
  
  nedges <- length(xp)   # sic
  
  # place x axis below polygon
  yp <- yp - min(yp) 

  # join vertex n to vertex 1
  nxt <- c(2:nedges, 1)

  # x step, WITH sign
  dx <- xp[nxt] - xp

  # average height 
  ym <- (yp + yp[nxt])/2
  
  -sum(dx * ym)
}

bdrylength.xypolygon <- function(polly) {
  verify.xypolygon(polly)
  xp <- polly$x
  yp <- polly$y
  nedges <- length(xp)
  nxt <- c(2:nedges, 1)
  dx <- xp[nxt] - xp
  dy <- yp[nxt] - yp
  sum(sqrt(dx^2 + dy^2))
}

reverse.xypolygon <- function(p, adjust=FALSE) {
  # reverse the order of vertices
  # (=> change sign of polygon)
  verify.xypolygon(p)
  p$x <- rev(p$x)
  p$y <- rev(p$y)
  if(adjust) {
    if(!is.null(p$hole)) p$hole <- !p$hole
    if(!is.null(p$area)) p$area <- -p$area
  }
  return(p)
}

overlap.xypolygon <- function(P, Q) {
  # compute area of overlap of two simple closed polygons 
  verify.xypolygon(P)
  verify.xypolygon(Q)
  
  xp <- P$x
  yp <- P$y
  np <- length(xp)
  nextp <- c(2:np, 1)

  xq <- Q$x
  yq <- Q$y
  nq <- length(xq)
  nextq <- c(2:nq, 1)

  # adjust y coordinates so all are nonnegative
  ylow <- min(c(yp,yq))
  yp <- yp - ylow
  yq <- yq - ylow

  area <- 0
  for(i in 1:np) {
    ii <- c(i, nextp[i])
    xpii <- xp[ii]
    ypii <- yp[ii]
    for(j in 1:nq) {
      jj <- c(j, nextq[j])
      area <- area +
        overlap.trapezium(xpii, ypii, xq[jj], yq[jj])
    }
  }
  return(area)
}

overlap.trapezium <- function(xa, ya, xb, yb, verb=FALSE) {
  # compute area of overlap of two trapezia
  # which have same baseline y = 0
  #
  # first trapezium has vertices
  # (xa[1], 0), (xa[1], ya[1]), (xa[2], ya[2]), (xa[2], 0).
  # Similarly for second trapezium
  
  # Test for vertical edges
  dxa <- diff(xa)
  dxb <- diff(xb)
  if(dxa == 0 || dxb == 0)
    return(0)

  # Order x coordinates, x0 < x1
  if(dxa > 0) {
    signa <- 1
    lefta <- 1
    righta <- 2
    if(verb) cat("A is positive\n")
  } else {
    signa <- -1
    lefta <- 2
    righta <- 1
    if(verb) cat("A is negative\n")
  }
  if(dxb > 0) {
    signb <- 1
    leftb <- 1
    rightb <- 2
    if(verb) cat("B is positive\n")
  } else {
    signb <- -1
    leftb <- 2
    rightb <- 1
    if(verb) cat("B is negative\n")
  }
  signfactor <- signa * signb # actually (-signa) * (-signb)
  if(verb) cat(paste("sign factor =", signfactor, "\n"))

  # Intersect x ranges
  x0 <- max(xa[lefta], xb[leftb])
  x1 <- min(xa[righta], xb[rightb])
  if(x0 >= x1)
    return(0)
  if(verb) {
    cat(paste("Intersection of x ranges: [", x0, ",", x1, "]\n"))
    abline(v=x0, lty=3)
    abline(v=x1, lty=3)
  }

  # Compute associated y coordinates
  slopea <- diff(ya)/diff(xa)
  y0a <- ya[lefta] + slopea * (x0-xa[lefta])
  y1a <- ya[lefta] + slopea * (x1-xa[lefta])
  slopeb <- diff(yb)/diff(xb)
  y0b <- yb[leftb] + slopeb * (x0-xb[leftb])
  y1b <- yb[leftb] + slopeb * (x1-xb[leftb])
  
  # Determine whether upper edges intersect
  # if not, intersection is a single trapezium
  # if so, intersection is a union of two trapezia

  yd0 <- y0b - y0a
  yd1 <- y1b - y1a
  if(yd0 * yd1 >= 0) {
    # edges do not intersect
    areaT <- (x1 - x0) * (min(y1a,y1b) + min(y0a,y0b))/2
    if(verb) cat(paste("Edges do not intersect\n"))
  } else {
    # edges do intersect
    # find intersection
    xint <- x0 + (x1-x0) * abs(yd0/(yd1 - yd0))
    yint <- y0a + slopea * (xint - x0)
    if(verb) {
      cat(paste("Edges intersect at (", xint, ",", yint, ")\n"))
      points(xint, yint, cex=2, pch="O")
    }
    # evaluate left trapezium
    left <- (xint - x0) * (min(y0a, y0b) + yint)/2
    # evaluate right trapezium
    right <- (x1 - xint) * (min(y1a, y1b) + yint)/2
    areaT <- left + right
    if(verb)
      cat(paste("Left area = ", left, ", right=", right, "\n"))    
  }

  # return area of intersection multiplied by signs 
  return(signfactor * areaT)
}


xypolygon2psp <- function(p, w, check=spatstat.options("checksegments")) {
  verify.xypolygon(p)
  n <- length(p$x)
  nxt <- c(2:n, 1)
  return(psp(p$x, p$y, p$x[nxt], p$y[nxt], window=w, check=check))
}
         
    
xypolyselfint <- function(p, eps=.Machine$double.eps,
                          proper=FALSE, yesorno=FALSE, checkinternal=FALSE) {
  verify.xypolygon(p)
  n <- length(p$x)
  verbose <- (n > 1000)
  if(verbose)
    cat(paste("[Checking polygon with", n, "edges..."))
  x0 <- p$x
  y0 <- p$y
  dx <- diff(x0[c(1:n,1)])
  dy <- diff(y0[c(1:n,1)])
  DUP <- spatstat.options("dupC")
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
                 answer=as.integer(integer(1)),
                 DUP=DUP)$answer
#                 PACKAGE="spatstat")
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
            ok=as.integer(integer(n^2)),
     DUP=DUP)
#     PACKAGE="spatstat")

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
  if(blowbyblow)
    cat(paste("Checking", npoly, ngettext(npoly, "polygon...", "polygons...")))

  dup <- self <- is.box <- logical(npoly)

  for(i in 1:npoly) {
    if(blowbyblow && npoly > 1)
      progressreport(i, npoly)
    Bi <- B[[i]]
    # check for duplicated vertices
    dup[i] <- any(duplicated(ppp(Bi$x, Bi$y, window=outerframe, check=FALSE)))
    if(dup[i] && blowbyblow)
      message(paste("Polygon", i, "contains duplicated vertices"))
    # check for self-intersection
    self[i] <- xypolyselfint(B[[i]], proper=TRUE, yesorno=TRUE)
    if(self[i] && blowbyblow)
      message(paste("Polygon", i, "is self-intersecting"))
    # check whether one of the current boundary polygons
    # is the bounding box itself (with + sign)
    is.box[i] <- (length(Bi$x) == 4) && (area.xypolygon(Bi) >= boxarea.mineps)
  }
  if(blowbyblow)
    cat("done.\n")
  
  if((ndup <- sum(dup)) > 0) {
    whinge <- paste(ngettext(ndup, "Polygon", "Polygons"),
                    if(npoly == 1) NULL else
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
                     if(npoly == 1) NULL else
                     commasep(which(self)),
                     ngettext(nself, "is", "are"),
                     "self-intersecting")
    notes <- c(notes, whinge)
    if(verbose) 
      message(whinge)
    err <- c(err, "self-intersection")
    answer <- FALSE
  }
  
  if((nbox <- sum(is.box)) > 1) {
    answer <- FALSE
    whinge <- paste("Polygons",
                    commasep(which(is.box)),
                    "coincide with the outer frame")
    notes <- c(notes, whinge)
    err <- c(err, "polygons duplicating the outer frame")
  }
  
  # check for crossings between different polygons
  cross <- matrix(FALSE, npoly, npoly)
  if(npoly > 1) {
    if(blowbyblow)
      cat(paste("Checking for cross-intersection between",
                npoly, "polygons..."))
    P <- lapply(B, xypolygon2psp, w=outerframe, check=FALSE)
    for(i in seq_len(npoly-1)) {
      if(blowbyblow)
        progressreport(i, npoly-1)
      Pi <- P[[i]]
      for(j in (i+1):npoly) {
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

simplify.xypolygon <- function(p, dmin) {
  verify.xypolygon(p)
  x <- p$x
  y <- p$y
  n <- length(x)
  if(n <= 3) return(p)
  dmin2 <- dmin^2
  # edge lengths: len2[i] is distance from i to i+1 
  len2 <- (x - c(x[-1], x[1]))^2 + (y - c(y[-1],y[1]))^2
  #
  while(n > 3 && any(len2 < dmin2)) {
    # delete the shortest edge
    kill <- which.min(len2)
    # edge from 'kill' to 'kill+1' will be removed 
    # Replacement vertex is midpoint of segment
    left <- if(kill == 1) n else (kill - 1)
    killplus1 <- if(kill == n) 1 else (kill + 1)
    right <- if(killplus1 == n) 1 else (killplus1 + 1)
    xmid <- (x[kill]+x[killplus1])/2
    ymid <- (y[kill]+y[killplus1])/2
    d2leftmid <- (xmid-x[left])^2+(ymid-y[left])^2
    d2midright <- (xmid-x[right])^2+(ymid-y[right])^2
    # adjust vectors: first replace segment endpoints without deleting
    x[kill] <- xmid
    y[kill] <- ymid
    x[killplus1] <- xmid
    y[killplus1] <- ymid
    len2[left] <- d2leftmid
    len2[kill] <- 0
    len2[killplus1] <- d2midright
    # now delete 
    x <- x[-kill]
    y <- y[-kill]
    len2 <- len2[-kill]
    n <- n-1
  }
  #
  p$x <- x
  p$y <- y
  p$area <- area.xypolygon(p[c("x","y")])
  return(p)
}

inside.triangle <- function(x, y, xx, yy) {
  # test whether points x[i], y[i] lie in triangle xx[1:3], yy[1:3]
  # using barycentric coordinates
  # vector 0 is edge from A to C
  v0x <- xx[3] - xx[1]
  v0y <- yy[3] - yy[1]
  # vector 1 is edge from A to B
  v1x <- xx[2] - xx[1]
  v1y <- yy[2] - yy[1]
  # vector 2 is from vertex A to point P
  v2x <- x - xx[1]
  v2y <- y - yy[1]
  # inner products
  dot00 <- v0x^2 + v0y^2
  dot01 <- v0x * v1x + v0y * v1y
  dot02 <- v0x * v2x + v0y * v2y
  dot11 <- v1x^2 + v1y^2
  dot12 <- v1x * v2x + v1y * v2y
  # unnormalised barycentric coordinates
  Denom <- dot00 * dot11 - dot01 * dot01
  u <- dot11 * dot02 - dot01 * dot12
  v <- dot00 * dot12 - dot01 * dot02
  # test
  return((u > 0) & (v > 0) & (u + v < Denom))
  
  
}
