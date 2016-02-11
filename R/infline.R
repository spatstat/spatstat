#
# infline.R
#
# Infinite lines
#
# $Revision: 1.21 $ $Date: 2016/02/11 10:17:12 $
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
  cat(paste(if(n > 1) n else NULL, "infinite ",
            ngettext(n, "line", "lines"), "\n"))
  print(as.data.frame(x), ...)
  return(invisible(NULL))
}

clip.infline <- function(L, win) {
  # clip a set of infinite straight lines to a window
  win <- as.owin(win)
  stopifnot(inherits(L, "infline"))
  # determine circumcircle of win
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owin(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
  # compute intersection points with circumcircle 
  p <- L$p
  theta <- L$theta
  hit <- (abs(p) < rmax)
  if(!any(hit)) 
    return(psp(numeric(0),numeric(0),numeric(0),numeric(0), window=win))
  p <- p[hit]
  theta <- theta[hit]
  q <- sqrt(rmax^2 - p^2)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
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
        if(h > yr[1] && h < yr[2]) 
          Zmat <- 2 * Zmat + (ymat > h)
      } else if(!is.na(v <- L[i, "v"])) {
        # vertical line
        if(v > xr[1] && v < xr[2])
          Zmat <- 2 * Zmat + (xmat < h)
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
      if(h < yr[1] || h > yr[2])
        Z <- NULL
      else {
        lower <- owin(xr, c(yr[1], h))
        upper <- owin(xr, c(h, yr[2]))
        Z <- tess(tiles=list(lower,upper), window=B)
      }
    } else if(!is.na(v <- L[i, "v"])) {
      # vertical line
      if(v < xr[1] || v > xr[2])
        Z <- NULL
      else {
        left <- owin(c(xr[1], v), yr)
        right <- owin(c(v, xr[2]), yr)
        Z <- tess(tiles=list(left,right), window=B)
      }
    } else {
      # generic line
      a <- L[i, "a"]
      b <- L[i, "b"]
      # Intersect with extended left and right sides of B
      yleft <- a + b * xr[1]
      yright <- a + b * xr[2]
      ylo <- min(yleft, yright, yr[1]) - 1
      yhi <- max(yleft, yright, yr[2]) + 1
      lower <- owin(poly=list(x=xr[c(1,1,2,2)],
                              y=c(yleft,ylo,ylo,yright)))
      upper <- owin(poly=list(x=xr[c(1,2,2,1)],
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



