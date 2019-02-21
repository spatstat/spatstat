#
#	$Revision: 1.27 $	$Date: 2019/02/21 01:18:21 $
#
#	Estimates of F, G and K for three-dimensional point patterns
#
#
#  ............ user interface .............................
#

K3est <- function(X, ...,
                  rmax=NULL, nrval=128,
                  correction=c("translation", "isotropic"),
                  ratio=FALSE)
{
  stopifnot(inherits(X, "pp3"))
  correction <- pickoption("correction", correction,
                           c(translation="translation",
			     trans="translation",
                             isotropic="isotropic",
                             iso="isotropic",
                             best="isotropic"),
                           multi=TRUE)
  trap.extra.arguments(..., .Context="In K3est")
  B <- X$domain
  if(is.null(rmax))
    rmax <- diameter(B)/2
  r <- seq(from=0, to=rmax, length.out=nrval)
  np <- npoints(X)
  denom <- np * (np-1)/volume(B)
  
  # this will be the output data frame
  K <- data.frame(r=r, theo= (4/3) * pi * r^3)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- ratfv(K, NULL, denom,
             "r", quote(K[3](r)), 
             "theo", NULL, c(0,rmax/2), c("r","{%s[%s]^{pois}}(r)"), desc,
             fname=c("K", "3"),
             ratio=ratio)

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # extract coordinates
  coo <- coords(X)
  
  if(any(correction %in% "translation")) {
    u <- k3engine(coo$x, coo$y, coo$z, flatbox,
                  rmax=rmax, nrval=nrval, correction="translation")
    K <- bind.ratfv(K,
                    data.frame(trans=u$num), u$denom,
                    "{hat(%s)[%s]^{trans}}(r)",
                    "translation-corrected estimate of %s",
                    "trans",
                    ratio=ratio)
  }
  if(any(correction %in% "isotropic")) {
    u <- k3engine(coo$x, coo$y, coo$z, flatbox,
                  rmax=rmax, nrval=nrval, correction="isotropic")
    K <- bind.ratfv(K,
                    data.frame(iso=u$num), u$denom,
                    "{hat(%s)[%s]^{iso}}(r)",
                    "isotropic-corrected estimate of %s",
                    "iso",
                    ratio=ratio)
  }
  # default is to display them all
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  return(K)
}
                  
G3est <- function(X, ...,
                  rmax=NULL, nrval=128,
                  correction=c("rs", "km", "Hanisch"))
{
  stopifnot(inherits(X, "pp3"))
  correction <- pickoption("correction", correction,
                           c(rs="rs",
                             border="rs",
                             km="km",
                             KM="km",
                             Hanisch="han",
                             hanisch="han",
                             best="km"),
                           multi=TRUE)
  trap.extra.arguments(..., .Context="In G3est")
  B <- X$domain
  if(is.null(rmax))
    rmax <- diameter(B)/2
  r <- seq(from=0, to=rmax, length.out=nrval)

  coo <- coords(X)
  lambda <- nrow(coo)/volume(B)
  
  # this will be the output data frame
  G <- data.frame(r=r, theo= 1 - exp( - lambda * (4/3) * pi * r^3))
  desc <- c("distance argument r", "theoretical Poisson %s")
  G <- fv(G, "r", substitute(G3(r), NULL),
          "theo", , c(0,rmax/2), c("r","%s[pois](r)"), desc, fname="G3")

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # collect four histograms for censored data
  u <- g3Cengine(coo$x, coo$y, coo$z, flatbox,
                 rmax=rmax, nrval=nrval)

  if("rs" %in% correction) 
    G <- bind.fv(G, data.frame(rs=u$rs), "%s[rs](r)",
                  "reduced sample estimate of %s",
                  "rs")
  if("km" %in% correction)
    G <- bind.fv(G, data.frame(km=u$km), "%s[km](r)",
                  "Kaplan-Meier estimate of %s",
                  "km")
  if("han" %in% correction) 
    G <- bind.fv(G, data.frame(han=u$han), "%s[han](r)",
                  "Normalised Hanisch estimate of %s",
                  "han")
  # default is to display them all
  formula(G) <- . ~ r
  unitname(G) <- unitname(X)
  return(G)
}

F3est <- function(X, ...,
                  rmax=NULL, nrval=128, vside=NULL,
                  correction=c("rs", "km", "cs"),
                  sphere=c("fudge", "ideal", "digital"))
{
  stopifnot(inherits(X, "pp3"))
  sphere <- match.arg(sphere)
  correction <- pickoption("correction", correction,
                           c(rs="rs",
                             border="rs",
                             km="km",
                             KM="km",
                             Kaplan="km",
                             cs="cs",
                             CS="cs",
                             best="km"),
                           multi=TRUE)
  trap.extra.arguments(..., .Context="In F3est")
  B <- X$domain
  if(is.null(rmax))
    rmax <- diameter(B)/2
  r <- seq(from=0, to=rmax, length.out=nrval)

  coo <- coords(X)
  vol <- volume(B)
  lambda <- nrow(coo)/vol

  # determine voxel size
  if(missing(vside)) {
    voxvol <- vol/spatstat.options("nvoxel")
    vside <- voxvol^(1/3)
    # ensure the shortest side is a whole number of voxels
    s <- shortside(B)
    m <- ceiling(s/vside)
    vside <- s/m
  }

  # compute theoretical value
  switch(sphere,
         ideal = {
           volsph <- (4/3) * pi * r^3
           spherename <- "ideal sphere"
         },
         fudge = {
           volsph <- 0.78 * (4/3) * pi * r^3
           spherename <- "approximate sphere"
         },
         digital = {
           volsph <- digital.volume(c(0, rmax), nrval, vside)
           spherename <- "digital sphere"
         })
  theo.desc <- paste("theoretical Poisson %s using", spherename)
           
  # this will be the output data frame
  FF <- data.frame(r     = r,
                   theo  = 1 - exp( - lambda * volsph))
  desc <- c("distance argument r", theo.desc)
  labl <- c("r","%s[pois](r)")
  FF <- fv(FF, "r", substitute(F3(r), NULL),
          "theo", , c(0,rmax/2), labl, desc, fname="F3")

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # go
  u <- f3Cengine(coo$x, coo$y, coo$z, flatbox,
                 rmax=rmax, nrval=nrval, vside=vside)

  if("rs" %in% correction) 
    FF <- bind.fv(FF, data.frame(rs=u$rs), "%s[rs](r)",
                  "reduced sample estimate of %s",
                  "rs")
  if("km" %in% correction)
    FF <- bind.fv(FF, data.frame(km=u$km), "%s[km](r)",
                  "Kaplan-Meier estimate of %s",
                  "km")
  if("cs" %in% correction)
    FF <- bind.fv(FF, data.frame(cs=u$cs), "%s[cs](r)",
                  "Chiu-Stoyan estimate of %s",
                  "cs")
  # default is to display them all
  formula(FF) <- . ~ r
  unitname(FF) <- unitname(X)
  return(FF)
}

pcf3est <- function(X, ...,
                    rmax=NULL, nrval=128,
                    correction=c("translation", "isotropic"),
                    delta=NULL, adjust=1, biascorrect=TRUE)
{
  stopifnot(inherits(X, "pp3"))
  correction <- pickoption("correction", correction,
                           c(translation="translation",
                             trans="translation",
                             isotropic="isotropic",
                             iso="isotropic",
                             best="isotropic"),
                           multi=TRUE)
  trap.extra.arguments(..., .Context="In pcf3est")
  B <- X$domain
  if(is.null(rmax))
    rmax <- diameter(B)/2
  r <- seq(from=0, to=rmax, length.out=nrval)

  if(is.null(delta)) {
    lambda <- npoints(X)/volume(B)
    delta <- adjust * 0.26/lambda^(1/3)
  }
  if(biascorrect) {
    # bias correction
    rondel <- r/delta
    biasbit <- ifelseAX(rondel > 1, 1, (3/4)*(rondel + 2/3 - (1/3)*rondel^3))
  }

  # this will be the output data frame
  g <- data.frame(r=r, theo=rep.int(1, length(r)))
  desc <- c("distance argument r", "theoretical Poisson %s")
  g <- fv(g, "r", quote(g[3](r)),
          "theo", , c(0,rmax/2),
          c("r", "{%s[%s]^{pois}}(r)"),
          desc, fname=c("g", "3"))

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # extract coordinates
  coo <- coords(X)
  
  if(any(correction %in% "translation")) {
    u <- pcf3engine(coo$x, coo$y, coo$z, flatbox,
                  rmax=rmax, nrval=nrval, correction="translation", delta=delta)
    gt <- u$f
    if(biascorrect)
      gt <- gt/biasbit
    g <- bind.fv(g, data.frame(trans=gt),
                 "{hat(%s)[%s]^{trans}}(r)",
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction %in% "isotropic")) {
    u <- pcf3engine(coo$x, coo$y, coo$z, flatbox,
                  rmax=rmax, nrval=nrval, correction="isotropic", delta=delta)
    gi <- u$f
    if(biascorrect)
      gi <- gi/biasbit
    g <- bind.fv(g, data.frame(iso=gi), 
                 "{hat(%s)[%s]^{iso}}(r)",
                 "isotropic-corrected estimate of %s",
                 "iso")
  }
  # default is to display them all
  formula(g) <- . ~ r
  unitname(g) <- unitname(X)
  attr(g, "delta") <- delta
  return(g)
}

#  ............ low level code ..............................
#
k3engine <- function(x, y, z, box=c(0,1,0,1,0,1),
                     rmax=1, nrval=100, correction="translation") 
{
  code <- switch(correction, translation=0, isotropic=1)
  res <- .C("RcallK3",
            as.double(x), as.double(y), as.double(z), 
            as.integer(length(x)),
            as.double(box[1L]), as.double(box[2L]), 
            as.double(box[3L]), as.double(box[4L]), 
            as.double(box[5L]), as.double(box[6L]), 
            as.double(0), as.double(rmax), 
            as.integer(nrval),
            f = as.double(numeric(nrval)),
            num = as.double(numeric(nrval)),
            denom = as.double(numeric(nrval)),
            as.integer(code),
            PACKAGE = "spatstat")
  return(list(range = c(0,rmax),
              f = res$f, num=res$num, denom=res$denom, 
              correction=correction))
}
#
#
#
g3engine <- function(x, y, z, box=c(0,1,0,1,0,1), 
                     rmax=1, nrval=10, correction="Hanisch G3") 
{
	code <- switch(correction, "minus sampling"=1, "Hanisch G3"=3)
	res <- .C("RcallG3",
		as.double(x), as.double(y), as.double(z), 
		as.integer(length(x)),
		as.double(box[1L]), as.double(box[2L]), 
		as.double(box[3L]), as.double(box[4L]), 
		as.double(box[5L]), as.double(box[6L]), 
		as.double(0), as.double(rmax), 
		as.integer(nrval),
		f = as.double(numeric(nrval)),
		num = as.double(numeric(nrval)),
		denom = as.double(numeric(nrval)),
		as.integer(code),
	  PACKAGE = "spatstat")
	return(list(range = c(0, rmax),
                    f = res$f, num=res$num, denom=res$denom, 
                    correction=correction))
}
#
#
f3engine <- function(x, y, z, box=c(0,1,0,1,0,1), 
	vside=0.05, 
	range=c(0,1.414), nval=25, correction="minus sampling") 
	
{
#
	code <- switch(correction, "minus sampling"=1, no=0)
	res <- .C("RcallF3",
		as.double(x), as.double(y), as.double(z), 
		as.integer(length(x)),
		as.double(box[1L]), as.double(box[2L]), 
		as.double(box[3L]), as.double(box[4L]), 
		as.double(box[5L]), as.double(box[6L]), 
		as.double(vside), 
		as.double(range[1L]), as.double(range[2L]),
		m=as.integer(nval),
		num = as.integer(integer(nval)),
		denom = as.integer(integer(nval)),
		as.integer(code),
	  PACKAGE = "spatstat")
	r <- seq(from=range[1L], to=range[2L], length.out=nval)
	f <- with(res, ifelseXB(denom > 0, num/denom, 1))

	return(list(r = r, f = f, num=res$num, denom=res$denom, 
		correction=correction))
}

f3Cengine <- function(x, y, z, box=c(0,1,0,1,0,1), 
	vside=0.05, rmax=1, nrval=25)
{
#
  res <- .C("RcallF3cen",
            as.double(x), as.double(y), as.double(z), 
            as.integer(length(x)),
            as.double(box[1L]), as.double(box[2L]), 
            as.double(box[3L]), as.double(box[4L]), 
            as.double(box[5L]), as.double(box[6L]), 
            as.double(vside), 
            as.double(0), as.double(rmax),
            m=as.integer(nrval),
            obs = as.integer(integer(nrval)),
            nco = as.integer(integer(nrval)),
            cen = as.integer(integer(nrval)),
            ncc = as.integer(integer(nrval)),
            upperobs = as.integer(integer(1L)),
            uppercen = as.integer(integer(1L)),
            PACKAGE = "spatstat")
  r <- seq(from=0, to=rmax, length.out=nrval)
  #
  obs <- res$obs
  nco <- res$nco
  cen <- res$cen
  ncc <- res$ncc
  upperobs <- res$upperobs
  uppercen <- res$uppercen
  #
  breaks <- breakpts.from.r(r)
  km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
  rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
  #
  ero <- eroded.volumes(as.box3(box), r)
  H <- cumsum(nco/ero)
  cs <- H/max(H[is.finite(H)])
  #
  return(list(rs=rs, km=km$km, hazard=km$lambda, cs=cs, r=r))
}

g3Cengine <- function(x, y, z, box=c(0,1,0,1,0,1), 
	rmax=1, nrval=25)
{
#
  res <- .C("RcallG3cen",
            as.double(x), as.double(y), as.double(z), 
            as.integer(length(x)),
            as.double(box[1L]), as.double(box[2L]), 
            as.double(box[3L]), as.double(box[4L]), 
            as.double(box[5L]), as.double(box[6L]), 
            as.double(0), as.double(rmax),
            m=as.integer(nrval),
            obs = as.integer(integer(nrval)),
            nco = as.integer(integer(nrval)),
            cen = as.integer(integer(nrval)),
            ncc = as.integer(integer(nrval)),
            upperobs = as.integer(integer(1L)),
            uppercen = as.integer(integer(1L)),
            PACKAGE = "spatstat")
  r <- seq(from=0, to=rmax, length.out=nrval)
  #
  obs <- res$obs
  nco <- res$nco
  cen <- res$cen
  ncc <- res$ncc
  upperobs <- res$upperobs
  uppercen <- res$uppercen
  #
  breaks <- breakpts.from.r(r)
  km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
  rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
  #
  ero <- eroded.volumes(as.box3(box), r)
  H <- cumsum(nco/ero)
  han <- H/max(H[is.finite(H)])
  return(list(rs=rs, km=km$km, hazard=km$lambda, han=han, r=r))
}

pcf3engine <- function(x, y, z, box=c(0,1,0,1,0,1),
                       rmax=1, nrval=100, correction="translation",
                       delta=rmax/10) 
{
  code <- switch(correction, translation=0, isotropic=1)
  res <- .C("Rcallpcf3",
            as.double(x), as.double(y), as.double(z), 
            as.integer(length(x)),
            as.double(box[1L]), as.double(box[2L]), 
            as.double(box[3L]), as.double(box[4L]), 
            as.double(box[5L]), as.double(box[6L]), 
            as.double(0), as.double(rmax), 
            as.integer(nrval),
            f = as.double(numeric(nrval)),
            num = as.double(numeric(nrval)),
            denom = as.double(numeric(nrval)),
            method=as.integer(code),
            delta=as.double(delta),
            PACKAGE = "spatstat")
	return(list(range = c(0,rmax),
                    f = res$f, num=res$num, denom=res$denom, 
                    correction=correction))
}
#
# ------------------------------------------------------------
# volume of a sphere (exact and approximate)
#

sphere.volume <- function(range=c(0,1.414), nval=10) 
{
  rr <- seq(from=range[1L], to=range[2L], length.out=nval)
  return( (4/3) * pi * rr^3)
}

digital.volume <- function(range=c(0, 1.414),  nval=25, vside= 0.05) 
{
#	Calculate number of points in digital sphere 
#	by performing distance transform for a single point
#	in the middle of a suitably large box
#
#	This takes EIGHT TIMES AS LONG as the corresponding empirical F-hat !!!
#
	w <- 2 * range[2L] + 2 * vside
#
	dvol <- .C("RcallF3",
                   as.double(w/2), as.double(w/2), as.double(w/2),
                   as.integer(1L),
                   as.double(0), as.double(w), 
                   as.double(0), as.double(w), 
                   as.double(0), as.double(w), 
                   as.double(vside),
                   as.double(range[1L]), as.double(range[2L]),
                   as.integer(nval),
                   num = as.integer(integer(nval)),
                   denom = as.integer(integer(nval)),
                   as.integer(0),
	                 PACKAGE = "spatstat")$num
#	
        (vside^3) * dvol 
      }

