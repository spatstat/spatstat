#	Gmulti.S
#
#	Compute estimates of nearest neighbour distance distribution functions
#	for multitype point patterns
#
#	S functions:	
#		Gcross                G_{ij}
#		Gdot		      G_{i\bullet}
#		Gmulti	              (generic)
#
#	$Revision: 4.42 $	$Date: 2014/10/24 00:22:30 $
#
################################################################################

"Gcross" <-		
function(X, i, j, r=NULL, breaks=NULL, ..., correction=c("rs", "km", "han"))
{
#	computes G_{ij} estimates
#
#	X		marked point pattern (of class 'ppp')
#	i,j		the two mark values to be compared
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
  X <- as.ppp(X)
  if(!is.marked(X, dfok=FALSE))
    stop(paste("point pattern has no", sQuote("marks")))
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
  if(missing(j)) j <- levels(marx)[2]
#  
  I <- (marx == i)
  if(sum(I) == 0) stop("No points are of type i")
        
  if(i == j)
    result <- Gest(X[I], r=r, breaks=breaks, ...)
  else {
    J <- (marx == j)
    if(sum(J) == 0) stop("No points are of type j")
    result <- Gmulti(X, I, J, r=r, breaks=breaks, disjoint=FALSE, ...,
                     correction=correction)
  }
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result,
               substitute(G[i,j](r), list(i=iname, j=jname)),
               c("G", paste0("list(", iname, ",", jname, ")")),
               new.yexp=substitute(G[list(i,j)](r),
                                   list(i=iname,j=jname)))
  return(result)
}	

"Gdot" <- 	
function(X, i, r=NULL, breaks=NULL, ..., correction=c("km","rs","han")) {
#  Computes estimate of 
#      G_{i\bullet}(t) = 
#  P( a further point of pattern in B(0,t)| a type i point at 0 )
#	
#	X		marked point pattern (of class ppp)
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
  X <- as.ppp(X)
  if(!is.marked(X))
    stop(paste("point pattern has no", sQuote("marks")))
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
  I <- (marx == i)
  if(sum(I) == 0) stop("No points are of type i")
  J <- rep.int(TRUE, X$n)	# i.e. all points
# 
  result <- Gmulti(X, I, J, r, breaks, disjoint=FALSE, ...,
                   correction=correction)
  iname <- make.parseable(paste(i))
  result <- rebadge.fv(result,
                  substitute(G[i ~ dot](r), list(i=iname)),
                  c("G", paste(iname, "~ symbol(\"\\267\")")),
                  new.yexp=substitute(G[i ~ symbol("\267")](r), list(i=iname)))
  return(result)
}	

	
##########

"Gmulti" <- 	
function(X, I, J, r=NULL, breaks=NULL, ..., disjoint=NULL,
         correction=c("rs", "km", "han")) {
#
#  engine for computing the estimate of G_{ij} or G_{i\bullet}
#  depending on selection of I, J
#  
#	X		marked point pattern (of class ppp)
#	
#	I,J		logical vectors of length equal to the number of points
#			and identifying the two subsets of points to be
#			compared.
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
  verifyclass(X, "ppp")
  W <- X$window
  npts <- npoints(X)
  areaW <- area(W)
# check I and J
  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop("No points satisfy condition I")
  if(nJ == 0) stop("No points satisfy condition J")

  if(is.null(disjoint))
    disjoint <- !any(I & J)
# choose correction(s)
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("rs", "km", "han")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             han="han",
                             Hanisch="han",
                             best="km"),
                           multi=TRUE)
#  determine breakpoints for r values
  lamJ <- nJ/areaW
  rmaxdefault <- rmax.rule("G", W, lamJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  brks <- breaks$val
  rmax <- breaks$max
  rvals <- breaks$r
  zeroes <- numeric(length(rvals))
# initialise fv object
  df <- data.frame(r=rvals, theo=1-exp(-lamJ * pi * rvals^2))
  fname <- c("G", "list(I,J)")
  Z <- fv(df, "r", quote(G[I,J](r)), "theo", . ~ r,
          c(0,rmax),
          c("r", makefvlabel(NULL, NULL, fname, "pois")),
          c("distance argument r", "theoretical Poisson %s"),
          fname=fname,
          yexp=quote(G[list(I,J)](r)))
#  "type I to type J" nearest neighbour distances
  XI <- X[I]
  XJ <- X[J]
  if(disjoint) 
    nnd <- nncross(XI, XJ, what="dist")
  else {
    seqnp <- seq_len(npts)
    iX <- seqnp[I]
    iY <- seqnp[J]
    nnd <- nncross(XI, XJ, iX, iY, what="dist")
  }
#  distance to boundary from each type i point
  bdry <- bdist.points(XI)
#  observations
  o <- pmin.int(nnd,bdry)
#  censoring indicators
  d <- (nnd <= bdry)
#
# calculate estimates
  
  if("none" %in% correction) {
    #  UNCORRECTED e.d.f. of nearest neighbour distances: use with care
    if(npts == 0)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax],breaks=breaks$val,plot=FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    Z <- bind.fv(Z, data.frame(raw=edf),
                 makefvlabel(NULL, "hat", fname, "raw"),
                 "uncorrected estimate of %s", "raw")
  }

  if("han" %in% correction) {
    # Hanisch style estimator
    if(npts == 0)
      G <- zeroes
    else {
      #  uncensored distances
      x <- nnd[d]
      #  weights
      a <- eroded.areas(W, rvals)
      # calculate Hanisch estimator
      h <- hist(x[x <= rmax], breaks=breaks$val, plot=FALSE)$counts
      G <- cumsum(h/a)
      G <- G/max(G[is.finite(G)])
    }
    # add to fv object
    Z <- bind.fv(Z, data.frame(han=G),
                 makefvlabel(NULL, "hat", fname, "han"),
                 "Hanisch estimate of %s",
                 "han")
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[G <= 0.9])
  }
  
  if(any(correction %in% c("rs", "km"))) {
    # calculate Kaplan-Meier and border correction (Reduced Sample) estimators
    if(npts == 0)
      result <- data.frame(rs=zeroes, km=zeroes, hazard=zeroes)
    else {
      result <- km.rs(o, bdry, d, breaks)
      result <- as.data.frame(result[c("rs", "km", "hazard")])
    }
    # add to fv object
    Z <- bind.fv(Z, result,
                 c(makefvlabel(NULL, "hat", fname, "bord"),
                   makefvlabel(NULL, "hat", fname, "km"),
                   "hazard(r)"),
                 c("border corrected estimate of %s",
                   "Kaplan-Meier estimate of %s",
                   "Kaplan-Meier estimate of hazard function lambda(r)"),
                 "km")
    # modify recommended plot range
    attr(Z, "alim") <- range(rvals[result$km <= 0.9])
  }
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  unitname(Z) <- unitname(X)
  return(Z)
}	


