#
# infline.R
#
# Infinite lines
#
# $Revision: 1.28 $ $Date: 2017/02/07 07:47:20 $
#

infline <- function(a=NULL, b=NULL, h=NULL, v=NULL, p=NULL, theta=NULL) {
  if(is.null(a) != is.null(b))
    stop("invalid specification of a,b")
  if(is.null(p) != is.null(theta))
    stop("invalid specification of p,theta")
  if(!is.null(h)) 
    out <- data.frame(a=h, b=0, h=h, v=NA, p=h, theta=pi/2)
  else if(!is.null(v)) 
    out <- data.frame(a=NA,b=NA,h=NA,v=v,p=v,theta=ifelseAB(v < 0, pi, 0))
  else if(!is.null(a)) {
    # a, b specified
    z <- data.frame(a=a,b=b)
    a <- z$a
    b <- z$b
    theta <- ifelseAX(b == 0, pi/2, atan(-1/b))
    theta <- theta %% pi
    p <- a * sin(theta)
    out <- data.frame(a=a, b=b,
                      h=ifelseXB(b==0, a, NA),
                      v=NA, p=p, theta=theta)
  } else if(!is.null(p)) {
    # p, theta specified
    z <- data.frame(p=p,theta=theta)
    p <- z$p
    theta <- z$theta
    theta <- theta %% (2*pi)
    if(any(reverse <- (theta >= pi))) {
      theta[reverse] <- theta[reverse] - pi
      p[reverse]     <- -p[reverse]
    }
    vert <- (theta == 0)
    horz <- (cos(theta) == 0)
    gene <- !(vert | horz)
    v <- ifelseXB(vert, p, NA)
    h <- ifelseXB(horz, p, NA)
    a <- ifelseXB(gene, p/sin(theta), NA)
    b <- ifelseXB(gene, -cos(theta)/sin(theta), NA)
    out <- data.frame(a=a,b=b,h=h,v=v,p=p,theta=theta)
  } else stop("No data given!")
  class(out) <- c("infline", class(out))
  return(out)
}

is.infline <- function(x) { inherits(x, "infline") }

plot.infline <- function(x, ...) {
  for(i in seq_len(nrow(x))) {
    xi <- as.list(x[i, 1:4])
    xi[sapply(xi, is.na)] <- NULL
    do.call(abline, append(xi, list(...)))
  }
  return(invisible(NULL))
}

print.infline <- function(x, ...) {
  n <- nrow(x)
  splat(n, "infinite", ngettext(n, "line", "lines"))
  print(as.data.frame(x), ...)
  return(invisible(NULL))
}

clip.infline <- function(L, win) {
  # clip a set of infinite straight lines to a window
  win <- as.owin(win)
  stopifnot(inherits(L, "infline"))
  nL <- nrow(L)
  if(nL == 0)
    return(psp(numeric(0),numeric(0),numeric(0),numeric(0), window=win))
  seqL <- seq_len(nL)
  # determine circumcircle of win
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owin(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
  # convert line coordinates to origin (xmid, ymid)
  p <- L$p
  theta <- L$theta
  co <- cos(theta)
  si <- sin(theta)
  p <- p - xmid * co - ymid * si
  # compute intersection points with circumcircle 
  hit <- (abs(p) < rmax)
  if(!any(hit)) 
    return(psp(numeric(0),numeric(0),numeric(0),numeric(0), window=win))
  p <- p[hit]
  theta <- theta[hit]
  q <- sqrt(rmax^2 - p^2)
  co <- co[hit]
  si <- si[hit]
  id <- seqL[hit]
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
           marks = factor(id, levels=seqL),
           window=boundbox, check=FALSE)
  # clip to window
  X <- X[win]
  return(X)
}
  
chop.tess <- function(X, L) {
  stopifnot(is.infline(L))
  stopifnot(is.tess(X)||is.owin(X))
  X <- as.tess(X)

  if(X$type == "image") {
    Xim <- X$image
    xr <- Xim$xrange
    yr <- Xim$yrange
    # extract matrices of pixel values and x, y coordinates
    Zmat <- as.integer(as.matrix(Xim))
    xmat <- rasterx.im(Xim)
    ymat <- rastery.im(Xim)
    # process lines
    for(i in seq_len(nrow(L))) {
      # line i chops window into two pieces
      if(!is.na(h <- L[i, "h"])) {
        # horizontal line
        if(h > yr[1L] && h < yr[2L]) 
          Zmat <- 2 * Zmat + (ymat > h)
      } else if(!is.na(v <- L[i, "v"])) {
        # vertical line
        if(v > xr[1L] && v < xr[2L])
          Zmat <- 2 * Zmat + (xmat < v)
      } else {
        # generic line y = a + bx
        a <- L[i, "a"]
        b <- L[i, "b"]
        Zmat <- 2 * Zmat + (ymat > a + b * xmat)
      }
    }
    # Now just put back as factor image
    Zim <- im(Zmat, xcol=Xim$xcol, yrow=Xim$yrow, unitname=unitname(Xim))
    Z <- tess(image=Zim)
    return(Z)
  }

  #---- polygonal computation --------
  # get bounding box
  B <- as.rectangle(as.owin(X))
  xr <- B$xrange
  yr <- B$yrange

  # get coordinates
  for(i in seq_len(nrow(L))) {
    # line i chops box B into two pieces
    if(!is.na(h <- L[i, "h"])) {
      # horizontal line
      if(h < yr[1L] || h > yr[2L])
        Z <- NULL
      else {
        lower <- owin(xr, c(yr[1L], h))
        upper <- owin(xr, c(h, yr[2L]))
        Z <- tess(tiles=list(lower,upper), window=B)
      }
    } else if(!is.na(v <- L[i, "v"])) {
      # vertical line
      if(v < xr[1L] || v > xr[2L])
        Z <- NULL
      else {
        left <- owin(c(xr[1L], v), yr)
        right <- owin(c(v, xr[2L]), yr)
        Z <- tess(tiles=list(left,right), window=B)
      }
    } else {
      # generic line
      a <- L[i, "a"]
      b <- L[i, "b"]
      # Intersect with extended left and right sides of B
      yleft <- a + b * xr[1L]
      yright <- a + b * xr[2L]
      ylo <- min(yleft, yright, yr[1L]) - 1
      yhi <- max(yleft, yright, yr[2L]) + 1
      lower <- owin(poly=list(x=xr[c(1L,1L,2L,2L)],
                              y=c(yleft,ylo,ylo,yright)))
      upper <- owin(poly=list(x=xr[c(1L,2L,2L,1L)],
                              y=c(yleft,yright,yhi,yhi)))
      Bplus <- owin(xr, c(ylo, yhi), unitname=unitname(B))
      Z <- tess(tiles=list(lower,upper), window=Bplus)
    }
    # intersect this simple tessellation with X
    if(!is.null(Z)) {
      X <- intersect.tess(X, Z)
      tilenames(X) <- paste("Tile", seq_len(length(tiles(X))))
    }
  }
  return(X)
}

whichhalfplane <- function(L, x, y=NULL) {
  verifyclass(L, "infline")
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  m <- length(x)
  n <- nrow(L)
  Z <- matrix(as.logical(NA_integer_), n, m)
  for(i in seq_len(n)) {
    if(!is.na(h <- L[i, "h"])) {
      #' horizontal line
      Z[i,] <- (y < h)
    } else if(!is.na(v <- L[i, "v"])) {
      #' vertical line
      Z[i,] <- (x < v)
    } else {
      #' generic line y = a + bx
      a <- L[i, "a"]
      b <- L[i, "b"]
      Z[i,] <- (y < a + b * x)
    }
  }
  return(Z)
}

rotate.infline <- function(X, angle=pi/2, ...) {
  if(nrow(X) == 0) return(X)
  Y <- with(X, infline(p = p, theta=theta + angle))
  return(Y)
}

shift.infline <- function(X, vec=c(0,0), ...) {
  if(nrow(X) == 0) return(X)
  vec <- as2vector(vec)
  Y <- with(X, infline(p = p + vec[1L] * cos(theta) + vec[2L] * sin(theta),
                       theta=theta))
  return(Y)
}

reflect.infline <- function(X) {
  if(nrow(X) == 0) return(X)
  Y <- with(X, infline(p = p,
                       theta=(theta + pi) %% (2 * pi)))
  return(Y)
}

flipxy.infline <- function(X) {
  if(nrow(X) == 0) return(X)
  Y <- with(X, infline(p = p,
                       theta=(pi/2 - theta) %% (2 * pi)))
  return(Y)
}

