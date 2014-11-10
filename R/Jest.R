#	Jest.S
#
#	Usual invocation to compute J function
#	if F and G are not required 
#
#	$Revision: 4.19 $	$Date: 2014/11/10 08:23:33 $
#
#
#
"Jest" <-
function(X, ..., eps=NULL, r=NULL, breaks=NULL, correction=NULL) {
  X <- as.ppp(X)
  W<- X$window
  rmaxdefault <- rmax.rule("J", W)
  brks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)$val
  # compute F and G 
  FF <- Fest(X, eps, breaks=brks, correction=correction)
  G <- Gest(X, breaks=brks, correction=correction)
  # initialise fv object
  rvals <- FF$r
  rmax  <- max(rvals)
#  Fvals <- FF[[attr(FF, "valu")]]
  Z <- fv(data.frame(r=rvals, theo=1),
          "r", substitute(J(r), NULL),
          "theo",
          . ~ r,
          c(0,rmax),
          c("r", "%s[pois](r)"),
          c("distance argument r", "theoretical Poisson %s"),
          fname="J")
  # compute J function estimates
  # this has to be done manually because of the mismatch between names
  ratio <- function(a, b) {
    result <- a/b
    result[ b == 0 ] <- NA
    result
  }
  Fnames <- names(FF)
  Gnames <- names(G)
  if("raw" %in% Gnames && "raw" %in% Fnames) {
    Jun <- ratio(1-G$raw, 1-FF$raw)
    Z <- bind.fv(Z, data.frame(un=Jun), "hat(%s)[un](r)",
                 "uncorrected estimate of %s", "un")
    attr(Z, "alim") <- range(rvals[FF$raw <= 0.9])
  }
  if("rs" %in% Gnames && "rs" %in% Fnames) {
    Jrs <- ratio(1-G$rs, 1-FF$rs)
    Z <- bind.fv(Z, data.frame(rs=Jrs), "hat(%s)[rs](r)",
                 "border corrected estimate of %s", "rs")
    attr(Z, "alim") <- range(rvals[FF$rs <= 0.9])
  }
  if("han" %in% Gnames && "cs" %in% Fnames) {
    Jhan <- ratio(1-G$han, 1-FF$cs)
    Z <- bind.fv(Z, data.frame(han=Jhan), "hat(%s)[han](r)",
                 "Hanisch-style estimate of %s", "han")
    attr(Z, "alim") <- range(rvals[FF$cs <= 0.9])
  }
  if("km" %in% Gnames && "km" %in% Fnames) {
    Jkm <- ratio(1-G$km, 1-FF$km)
    Z <- bind.fv(Z, data.frame(km=Jkm), "hat(%s)[km](r)",
                 "Kaplan-Meier estimate of %s", "km")
    attr(Z, "alim") <- range(rvals[FF$km <= 0.9])
  }
  if("hazard" %in% Gnames && "hazard" %in% Fnames) {
    Jhaz <- G$hazard - FF$hazard
    Z <- bind.fv(Z, data.frame(hazard=Jhaz), "hazard(r)",
                 "Kaplan-Meier estimate of derivative of log(%s)")
  }
# set default plotting values and order
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  
# add more info        
  attr(Z, "F") <- FF
  attr(Z, "G") <- G

  unitname(Z) <- unitname(X)
  return(Z)
}

