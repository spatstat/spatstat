#
#   pp3.R
#
#  class of three-dimensional point patterns in rectangular boxes
#
#  $Revision: 1.22 $  $Date: 2015/05/04 05:02:31 $
#

box3 <- function(xrange=c(0,1), yrange=xrange, zrange=yrange, unitname=NULL) {
  stopifnot(is.numeric(xrange) && length(xrange) == 2 && diff(xrange) > 0)
  stopifnot(is.numeric(yrange) && length(yrange) == 2 && diff(yrange) > 0)
  stopifnot(is.numeric(zrange) && length(zrange) == 2 && diff(zrange) > 0)
  out <- list(xrange=xrange, yrange=yrange, zrange=zrange,
              units=as.units(unitname))
  class(out) <- "box3"
  return(out)
}

as.box3 <- function(...) {
  a <- list(...)
  n <- length(a)
  if(n == 0)
    stop("No arguments given")
  if(n == 1) {
    a <- a[[1]]
    if(inherits(a, "box3"))
      return(a)
    if(inherits(a, "pp3"))
      return(a$domain)
    if(inherits(a, "boxx")){
      if(ncol(a$ranges)==3)
        return(box3(a$ranges[,1], a$ranges[,2], a$ranges[,3]))
      stop("Supplied boxx object does not have dimension three")
    }
    if(inherits(a, "ppx"))
      return(as.box3(a$domain))
    if(is.numeric(a)) {
      if(length(a) == 6)
        return(box3(a[1:2], a[3:4], a[5:6]))
      stop(paste("Don't know how to interpret", length(a), "numbers as a box"))
    }
    if(!is.list(a))
      stop("Don't know how to interpret data as a box")
  }
  return(do.call("box3", a))
}

print.box3 <- function(x, ...) {
  bracket <- function(z) paste("[",
                               paste(signif(z, 5), collapse=", "),
                               "]", sep="")
  v <- paste(unlist(lapply(x[1:3], bracket)), collapse=" x ")
  s <- summary(unitname(x))
  splat("Box:", v, s$plural, s$explain)
  invisible(NULL)
}

unitname.box3 <- function(x) { x$units }

"unitname<-.box3" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

eroded.volumes <- function(x, r) { UseMethod("eroded.volumes") }

eroded.volumes.box3 <- function(x, r) {
  b <- as.box3(x)
  ax <- pmax.int(0, diff(b$xrange) - 2 * r)
  ay <- pmax.int(0, diff(b$yrange) - 2 * r)
  az <- pmax.int(0, diff(b$zrange) - 2 * r)
  ax * ay * az
}

shortside <- function(x) { UseMethod("shortside") }

shortside.box3 <- function(x) {
  min(sidelengths(x))
}

sidelengths <- function(x) { UseMethod("sidelengths") }

sidelengths.box3 <- function(x) {
  with(x, c(diff(xrange), diff(yrange), diff(zrange)))
}

bounding.box3 <- function(...) {
  wins <- list(...)
  boxes <- lapply(wins, as.box3)
  xr <- range(unlist(lapply(boxes, getElement, name="xrange")))
  yr <- range(unlist(lapply(boxes, getElement, name="yrange")))
  zr <- range(unlist(lapply(boxes, getElement, name="zrange")))
  box3(xr, yr, zr)
}

pp3 <- function(x, y, z, ...) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(z)) 
  b <- as.box3(...)
  out <- ppx(data=data.frame(x=x,y=y,z=z), domain=b)
  class(out) <- c("pp3", class(out))
  return(out)
}

domain.pp3 <- function(X, ...) { X$domain }

is.pp3 <- function(x) { inherits(x, "pp3") }

npoints.pp3 <- function(x) { nrow(x$data) }

print.pp3 <- function(x, ...) {
  splat("Three-dimensional point pattern")
  sd <- summary(x$data)
  np <- sd$ncases
  splat(np, ngettext(np, "point", "points"))
  print(x$domain)
  invisible(NULL)
}

summary.pp3 <- function(object, ...) {
  sd <- summary(object$data)
  np <- sd$ncases
  dom <- object$domain
  v <- volume.box3(dom)
  u <- summary(unitname(dom))
  intens <- np/v
  out <-  list(np=np, sumdat=sd, dom=dom, v=v, u=u, intensity=intens)
  class(out) <- "summary.pp3"
  return(out)
}

print.summary.pp3 <- function(x, ...) {
  splat("Three-dimensional point pattern")
  splat(x$np, ngettext(x$np, "point", "points"))
  print(x$dom)
  u <- x$u
  v <- x$v
  splat("Volume", v, "cubic",
        if(v == 1) u$singular else u$plural,
        u$explain)
  splat("Average intensity", x$intensity,
        "points per cubic", u$singular, u$explain)
  invisible(NULL)
}

plot.pp3 <- function(x, ..., eye=NULL, org=NULL, theta=25, phi=15,
                     type=c("p", "n", "h")) {
  xname <- short.deparse(substitute(x))
  type <- match.arg(type)
  coo <- as.matrix(coords(x))
  xlim <- x$domain$xrange
  ylim <- x$domain$yrange
  zlim <- x$domain$zrange
  if(is.null(org)) org <- c(mean(xlim), mean(ylim), mean(zlim))
  if(is.null(eye)) {
    theta <- theta * pi/180
    phi   <- phi * pi/180
    d <- 2 * diameter(x$domain)
    eye <- org + d * c(cos(phi) * c(sin(theta), -cos(theta)), sin(phi))
  }
  deefolts <- spatstat.options('par.pp3')
  ## determine default eye position and centre of view
  do.call(plot3Dpoints,
          resolve.defaults(list(xyz=coo, eye=eye, org=org, type=type),
                           list(...),
                           deefolts,
                           list(main=xname,
                                xlim=x$domain$xrange,
                                ylim=x$domain$yrange,
                                zlim=x$domain$zrange)))
}

"[.pp3" <- function(x, i, drop=FALSE, ...) {
  answer <- NextMethod("[")
  if(is.ppx(answer))
    class(answer) <- c("pp3", class(answer))
  return(answer)
}
  
unitname.pp3 <- function(x) { unitname(x$domain) }

"unitname<-.pp3" <- function(x, value) {
  d <- x$domain
  unitname(d) <- value
  x$domain <- d
  return(x)
}

diameter.box3 <- function(x) {
  stopifnot(inherits(x, "box3"))
  with(x, sqrt(diff(xrange)^2+diff(yrange)^2+diff(zrange)^2))
}

volume <- function(x) { UseMethod("volume") }

volume.box3 <- function(x) {
  stopifnot(inherits(x, "box3"))
  with(x, prod(diff(xrange), diff(yrange), diff(zrange)))
}

runifpoint3 <- function(n, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  result <- vector(mode="list", length=nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for(i in 1:nsim) {
    x <- with(dd, runif(n, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(n, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(n, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

rpoispp3 <- function(lambda, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  v <- volume(domain)
  if(!(is.numeric(lambda) && length(lambda) == 1))
    stop("lambda must be a single numeric value")
  np <- rpois(nsim, lambda * v)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    ni <- np[i]
    x <- with(dd, runif(ni, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(ni, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(ni, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}



