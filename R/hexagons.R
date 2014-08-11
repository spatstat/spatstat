## hexagons.R
## $Revision: 1.3 $ $Date: 2014/01/09 06:53:01 $

hexgrid <- function(W, s, offset=c(0,0), trim=TRUE) {
  W <- as.owin(W)
  check.1.real(s)
  stopifnot(s > 0)
  R <- grow.rectangle(as.rectangle(W), s)
  p0 <- centroid.owin(R)
  p0 <- shiftxy(p0, offset)
  ## 'even' points
  hstep <- 3 * s
  vstep <- sqrt(3) * s
  xeven <- prolongseq(p0$x, R$xrange, step=hstep)
  yeven <- prolongseq(p0$y, R$yrange, step=vstep)
  xyeven <- expand.grid(x=xeven, y=yeven)
  ## 'odd' points
  xodd <- prolongseq(p0$x + hstep/2, R$xrange, step=hstep)
  yodd <- prolongseq(p0$y + vstep/2, R$yrange, step=vstep)
  xyodd <- expand.grid(x=xodd, y=yodd)
  ##
  xy <- concatxy(xyeven, xyodd)
  XY <- as.ppp(xy, W=R)
  ##
  if(trim) return(XY[W])
  ok <- inside.owin(XY, w=dilation.owin(W, s))
  return(XY[ok])
}

hextess <- function(W, s, offset=c(0,0), trim=TRUE) {
  W <- as.owin(W)
  G <- hexgrid(W, s, offset, trim=FALSE)
  if(trim && is.mask(W)) {
    ## Result is a pixel image tessellation
    ## Determine pixel resolution by extending 'W' to larger domain of 'G'
    rasta <- harmonise.im(as.im(1, W), as.owin(G))[[1]]
    rasta <- as.mask(rasta)
    ## Tweak G to have mask window
    G$window <- rasta
    ##
    img <- nnmap(G, what="which")
    result <- tess(image=img)
    return(result)
  }
  ## Result is a polygonal tessellation
  Gxy <- as.matrix(as.data.frame(G))
  n <- nrow(Gxy)
  ## Hexagon centred at origin
  hex0 <- disc(npoly=6, radius=s)
  ## Form hexagons
  hexes <- vector(mode="list", length=n)
  for(i in 1:n) 
    hexes[[i]] <- shift(hex0, Gxy[i,])
  ## Determine whether tiles intersect window wholly or partly
  suspect <- rep(TRUE, n)
  GW <- G[W]
  GinW <- inside.owin(G, w=W) 
  suspect[GinW] <- (bdist.points(GW) <= s)
  ## Compute intersection of tiles with window
  trimmed <- hexes
  trimmed[suspect] <- trimmed.suspect <- 
    lapply(trimmed[suspect], intersect.owin, B=W, fatal=FALSE)
  nonempty <- rep(TRUE, n)
  nonempty[suspect] <- !unlist(lapply(trimmed.suspect, is.empty))
  if(trim) {
    ## return the tiles intersected with W
    result <- tess(tiles=trimmed[nonempty], window=W)
  } else {
    ## return the tiles that have nonempty intersection with W
    result <- tess(tiles=hexes[nonempty])
  }
  return(result)
}


  
  
  
