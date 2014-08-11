#	Jmulti.S
#
#	Usual invocations to compute multitype J function(s)
#	if F and G are not required 
#
#	$Revision: 4.36 $	$Date: 2013/04/25 06:37:43 $
#
#
#
"Jcross" <-
function(X, i, j, eps=NULL, r=NULL, breaks=NULL, ..., correction=NULL) {
#
#       multitype J function J_{ij}(r)
#  
#	X:		point pattern (an object of class 'ppp')
#       i, j:           types for which J_{i,j}(r) is calculated  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
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
  if(missing(j)) j <- levels(marx)[2]
#
  I <- (marx == i)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))
#        
  if(i == j)
    result <- Jest(X[I], eps=eps, r=r, breaks=breaks, correction=correction)
  else {
    J <- (marx == j)
    result <- Jmulti(X, I, J,
                     eps=eps, r=r, breaks=breaks, disjoint=TRUE,
                     correction=correction)
  }
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result,
               substitute(J[i,j](r),
                          list(i=iname,j=jname)),
               sprintf("J[list(%s,%s)]", iname, jname),
               new.yexp=substitute(J[list(i,j)](r),
                                      list(i=iname,j=jname)))
  return(result)
}

"Jdot" <-
function(X, i, eps=NULL, r=NULL, breaks=NULL, ..., correction=NULL) {
#  
#    multitype J function J_{i\dot}(r)
#  
#	X:		point pattern (an object of class 'ppp')
#       i:              mark i for which we calculate J_{i\cdot}(r)  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
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
#  
  I <- (marx == i)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))          
  J <- rep.int(TRUE, X$n)
#  
  result <- Jmulti(X, I, J,
                   eps=eps, r=r, breaks=breaks, disjoint=FALSE,
                   correction=correction)
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(J[i ~ dot](r), list(i=iname)),
               paste("J[", iname, "~ symbol(\"\\267\")]"),
               new.yexp=substitute(J[i ~ symbol("\267")](r), list(i=iname)))
  return(result)
}

"Jmulti" <- 	
function(X, I, J, eps=NULL, r=NULL, breaks=NULL, ..., disjoint=NULL,
         correction=NULL) {
#  
#    multitype J function (generic engine)
#  
#	X		marked point pattern (of class ppp)
#	
#	I,J		logical vectors of length equal to the number of points
#			and identifying the two subsets of points to be
#			compared.
#  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#  
#
  X <- as.ppp(X)
  W<- X$window
  rmaxdefault <- rmax.rule("J", W)
  brks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)$val
  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
  FJ <- Fest(X[J], eps, breaks=brks, correction=correction)
  GIJ <- Gmulti(X, I, J, breaks=brks, disjoint=disjoint, correction=correction)
  rvals <- FJ$r
  Fnames <- names(FJ)
  Gnames <- names(GIJ)
  bothnames <- Fnames[Fnames %in% Gnames]
  # initialise fv object
  alim <- attr(FJ, "alim")
  Z <- fv(data.frame(r=rvals, theo=1),
          "r", substitute(Jmulti(r), NULL), "theo",
          . ~ r, alim,
          c("r", "{%s^{pois}}(r)"),
          c("distance argument r", "theoretical Poisson %s"),
          fname="J[multi]")
  # add pieces manually
  ratio <- function(a, b) {
    result <- a/b
    result[ b == 0 ] <- NA
    result
  }
  if("raw" %in% bothnames) {
    Jun <- ratio(1-GIJ$raw, 1-FJ$raw)
    Z <- bind.fv(Z, data.frame(un=Jun), "hat(%s^{un})(r)",
                 "uncorrected estimate of %s", "un")
  }
  if("rs" %in% bothnames) {
    Jrs <- ratio(1-GIJ$rs, 1-FJ$rs)
    Z <- bind.fv(Z, data.frame(rs=Jrs), "hat(%s^{rs})(r)",
                 "border corrected estimate of %s", "rs")
  }
  if("han" %in% Gnames && "cs" %in% Fnames) {
    Jhan <- ratio(1-GIJ$han, 1-FJ$cs)
    Z <- bind.fv(Z, data.frame(han=Jhan), "hat(%s^{han})(r)",
                 "Hanisch-style estimate of %s", "han")
  }
  if("km" %in% bothnames) {
    Jkm <- ratio(1-GIJ$km, 1-FJ$km)
    Z <- bind.fv(Z, data.frame(km=Jkm), "hat(%s^{km})(r)",
                 "Kaplan-Meier estimate of %s", "km")
    if("hazard" %in% names(GIJ) && "hazard" %in% names(FJ)) {
      Jhaz <- GIJ$hazard - FJ$hazard
      Z <- bind.fv(Z, data.frame(hazard=Jhaz), "hazard(r)",
                   "Kaplan-Meier estimate of derivative of log(%s)")
    } 
  }
# set default plotting values and order
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
# add other info
  attr(Z, "G") <- GIJ
  attr(Z, "F") <- FJ
  unitname(Z) <- unitname(X)
  return(Z)
}
