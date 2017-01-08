#
#    xypolygon.R
#
#    $Revision: 1.67 $    $Date: 2017/01/08 00:38:10 $
#
#    low-level functions defined for polygons in list(x,y) format
#
verify.xypolygon <- function(p, fatal=TRUE) {
  whinge <- NULL
  if(!is.list(p) || !all(c("x","y") %in% names(p)))
    whinge <- "polygon must be a list with components x and y"
  else if(is.null(p$x) || is.null(p$y) || !is.numeric(p$x) || !is.numeric(p$y))
    whinge <- "components x and y must be numeric vectors"
  else if(anyNA(p$x) || anyNA(p$y))
    whinge <- "x and y coordinates must not contain NA values"
  else if(length(p$x) != length(p$y))
    whinge <- "lengths of x and y vectors unequal"
  else if(length(p$x) < 3)
    whinge <- "need at least 3 vertices for a polygon"
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
  z <- .C("Cmatchxy",
          na=as.integer(npts),
          xa=as.double(x),
          ya=as.double(y),
          nb=as.integer(nedges),
          xb=as.double(xp),
          yb=as.double(yp),
          match=as.integer(integer(npts)))
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
             score <- temp$score/2
             on.boundary <- as.logical(temp$onbndry)
           },
# Fortran code removed!           
#           Fortran={
#             #------------------ call Fortran routine ------------------
#             temp <- DOTFortran("inxyp",
#                              x=as.double(x),
#                              y=as.double(y),
#                              xp=as.double(xp),
#                              yp=as.double(yp),
#                              npts=as.integer(npts),
#                              nedges=as.integer(nedges),
#                              score=as.double(score),
#                              onbndry=as.logical(on.boundary))
#             score <- temp$score
#             on.boundary <- temp$onbndry
#           },
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


is.hole.xypolygon <- function(polly) {
  h <- polly$hole
  if(is.null(h))
    h <- (Area.xypolygon(polly) < 0)
  return(h)
}
  
Area.xypolygon <- function(polly) {
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
  p$area <- Area.xypolygon(p[c("x","y")])
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
