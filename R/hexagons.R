## hexagons.R
## $Revision: 1.6 $ $Date: 2017/02/07 07:35:32 $

hexgrid <- function(W, s, offset=c(0,0), origin=NULL, trim=TRUE) {
  W <- as.owin(W)
  check.1.real(s)
  stopifnot(s > 0)
  hstep <- 3 * s
  vstep <- sqrt(3) * s
  R <- grow.rectangle(as.rectangle(W), hstep)
  xr <- R$xrange
  yr <- R$yrange
  ## initial positions for 'odd' and 'even grids
  p0 <- as2vector(origin %orifnull% centroid.owin(R))
  p0 <- p0 + as2vector(offset)
  q0 <- p0 + c(hstep, vstep)/2
  ## 'even' points
  p0 <- c(startinrange(p0[1L], hstep, xr),
          startinrange(p0[2L], vstep, yr))
  if(!anyNA(p0)) {
    xeven <- prolongseq(p0[1L], xr, step=hstep)
    yeven <- prolongseq(p0[2L], yr, step=vstep)
    xyeven <- expand.grid(x=xeven, y=yeven)
  } else xyeven <- list(x=numeric(0), y=numeric(0))
  ## 'odd' points
  q0 <- c(startinrange(q0[1L], hstep, xr),
          startinrange(q0[2L], vstep, yr))
  if(!anyNA(q0)) {
    xodd <- prolongseq(q0[1L], xr, step=hstep)
    yodd <- prolongseq(q0[2L], yr, step=vstep)
    xyodd <- expand.grid(x=xodd, y=yodd)
  } else xyodd <- list(x=numeric(0), y=numeric(0))
  ##
  xy <- concatxy(xyeven, xyodd)
  XY <- as.ppp(xy, W=R)
  ##
  if(trim) return(XY[W])
  ok <- inside.owin(XY, w=dilation.owin(W, s))
  return(XY[ok])
}

hextess <- function(W, s, offset=c(0,0), origin=NULL, trim=TRUE) {
  W <- as.owin(W)
  G <- hexgrid(W=W, s=s, offset=offset, origin=origin, trim=FALSE)
  if(trim && is.mask(W)) {
    ## Result is a pixel image tessellation
    ## Determine pixel resolution by extending 'W' to larger domain of 'G'
    rasta <- harmonise.im(as.im(1, W), as.owin(G))[[1L]]
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


  
  
  
