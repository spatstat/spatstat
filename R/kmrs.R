#
#	kmrs.S
#
#	S code for Kaplan-Meier, Reduced Sample and Hanisch
#	estimates of a distribution function
#	from _histograms_ of censored data.
#
#	kaplan.meier()
#	reduced.sample()
#       km.rs()
#
#	$Revision: 3.26 $	$Date: 2013/06/27 08:59:16 $
#
#	The functions in this file produce vectors `km' and `rs'
#	where km[k] and rs[k] are estimates of F(breaks[k+1]),
#	i.e. an estimate of the c.d.f. at the RIGHT endpoint of the interval.
#

"kaplan.meier" <-
function(obs, nco, breaks, upperobs=0) {
#	obs: histogram of all observations : min(T_i,C_i)
#	nco: histogram of noncensored observations : T_i such that T_i <= C_i
# 	breaks: breakpoints (vector or 'breakpts' object, see breaks.S)
#       upperobs: number of observations beyond rightmost breakpoint
#  
        breaks <- as.breakpts(breaks)

	n <- length(obs)
	if(n != length(nco)) 
		stop("lengths of histograms do not match")
	check.hist.lengths(nco, breaks)
#
#	
#   reverse cumulative histogram of observations
	d <- revcumsum(obs) + upperobs
#
#  product integrand
	s <- ifelseXB(d > 0, 1 - nco/d, 1)
#
	km <- 1 - cumprod(s)
#  km has length n;  km[i] is an estimate of F(r) for r=breaks[i+1]
#	
	widths <- diff(breaks$val)
        lambda <- numeric(n)
        pos <- (s > 0)
        lambda[pos] <- -log(s[pos])/widths[pos]
#  lambda has length n; lambda[i] is an estimate of
#  the average of \lambda(r) over the interval (breaks[i],breaks[i+1]).
#	
	return(list(km=km, lambda=lambda))
}

"reduced.sample" <-
function(nco, cen, ncc, show=FALSE, uppercen=0)
#	nco: histogram of noncensored observations: T_i such that T_i <= C_i
#	cen: histogram of all censoring times: C_i
#	ncc: histogram of censoring times for noncensored obs:
#		C_i such that T_i <= C_i
#
#	Then nco[k] = #{i: T_i <= C_i, T_i \in I_k}
#	     cen[k] = #{i: C_i \in I_k}
#	     ncc[k] = #{i: T_i <= C_i, C_i \in I_k}.
#
#       The intervals I_k must span an interval [0,R] beginning at 0.
#       If this interval did not include all censoring times,
#       then `uppercen' must be the number of censoring times
#       that were not counted in 'cen'.
{
	n <- length(nco)
	if(n != length(cen) || n != length(ncc))
		stop("histogram lengths do not match")
#
#	denominator: reverse cumulative histogram of censoring times
#		denom(r) = #{i : C_i >= r}
#	We compute 
#		cc[k] = #{i: C_i > breaks[k]}	
#	except that > becomes >= for k=0.
#
	cc <- revcumsum(cen) + uppercen
#
#
#	numerator
#	#{i: T_i <= r <= C_i }
#	= #{i: T_i <= r, T_i <= C_i} - #{i: C_i < r, T_i <= C_i}
#	We compute
#		u[k] = #{i: T_i <= C_i, T_i <= breaks[k+1]}
#			- #{i: T_i <= C_i, C_i <= breaks[k]}
#		     = #{i: T_i <= C_i, C_i > breaks[k], T_i <= breaks[k+1]}
#	this ensures that numerator and denominator are 
#	comparable, u[k] <= cc[k] always.
#
	u <- cumsum(nco) - c(0,cumsum(ncc)[1:(n-1)])
	rs <- u/cc
#
#	Hence rs[k] = u[k]/cc[k] is an estimator of F(r) 
#	for r = breaks[k+1], i.e. for the right hand end of the interval.
#
        if(!show)
          return(rs)
        else
          return(list(rs=rs, numerator=u, denominator=cc))
}

"km.rs" <-
function(o, cc, d, breaks) {
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
# compile histograms (breakpoints may not span data)
  obs <- whist( o,     breaks=bval)
  nco <- whist( o[d],  breaks=bval)
  cen <- whist( cc,    breaks=bval)
  ncc <- whist( cc[d], breaks=bval)
# number of observations exceeding largest breakpoint
  upperobs <- attr(obs, "high")
  uppercen <- attr(cen, "high")
# go
  km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
  rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
#
  return(list(rs=rs, km=km$km, hazard=km$lambda,
              r=breaks$r, breaks=bval))
}

"km.rs.opt" <-
function(o, cc, d, breaks, KM=TRUE, RS=TRUE) {
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
  out <- list(r=breaks$r, breaks=bval)
  if(KM || RS)
    nco <- whist( o[d],  breaks=bval)
  if(KM) {
    obs <- whist( o,     breaks=bval)
    upperobs <- attr(obs, "high")
    km <- kaplan.meier(obs, nco, breaks, upperobs=upperobs)
    out <- append(list(km=km$km, hazard=km$lambda), out)
  }
  if(RS) {
    cen <- whist( cc,    breaks=bval)
    ncc <- whist( cc[d], breaks=bval)
    uppercen <- attr(cen, "high")
    rs <- reduced.sample(nco, cen, ncc, uppercen=uppercen)
    out <- append(list(rs=rs), out)
  }
  return(out)
}


censtimeCDFest <- function(o, cc, d, breaks, ...,
                           KM=TRUE, RS=TRUE, HAN=TRUE, RAW=TRUE,
                           han.denom=NULL, tt=NULL, pmax=0.9) {
# Histogram-based estimation of cumulative distribution function
# of lifetimes subject to censoring.
#	o: censored lifetimes min(T_i,C_i)
#	cc: censoring times C_i
#	d: censoring indicators 1(T_i <= C_i)
#	breaks: histogram breakpoints (vector or 'breakpts' object)
#       han.denom: denominator (eroded area) for each value of r
#       tt: uncensored lifetimes T_i, if known  
  breaks <- as.breakpts(breaks)
  bval <- breaks$val
  rval <- breaks$r
  rmax <- breaks$max
  # Kaplan-Meier and/or Reduced Sample
  out <- km.rs.opt(o, cc, d, breaks, KM=KM, RS=RS)
  # convert to data frame
  out$breaks <- NULL
  df <- as.data.frame(out)
  # Raw ecdf of observed lifetimes if available
  if(RAW && !is.null(tt)) {
    h <- whist(tt[tt <= rmax], breaks=bval)
    df <- cbind(df, data.frame(raw=cumsum(h)/length(tt)))
  }
  # Hanisch
  if(HAN) {
    if(is.null(han.denom))
      stop("Internal error: missing denominator for Hanisch estimator")
    if(length(han.denom) != length(rval))
      stop(paste("Internal error:",
                 "length(han.denom) =", length(han.denom),
                 "!=", length(rval), "= length(rvals)"))
    #  uncensored distances
    x <- o[d]
    # calculate Hanisch estimator
    h <- whist(x[x <= rmax], breaks=bval)
    H <- cumsum(h/han.denom)
    df <- cbind(df, data.frame(han=H/max(H[is.finite(H)])))
  }
  # determine appropriate plotting range
  bestest <- if(KM) "km" else if(HAN) "han" else if(RS) "rs" else "raw"
  alim <- range(df$r[df[[bestest]] <= pmax])
  # convert to fv object
  nama <-  c("r",  "km", "hazard", "han", "rs", "raw")
  avail <- c(TRUE,  KM,  KM,       HAN,   RS,   RAW)
  iscdf <- c(FALSE, TRUE, FALSE,   TRUE,  TRUE, TRUE)
  labl <- c("r", "hat(%s)[km](r)", "lambda(r)", "hat(%s)[han](r)",
            "hat(%s)[bord](r)", "hat(%s)[raw](r)")[avail]
  desc <- c("distance argument r",
            "Kaplan-Meier estimate of %s",
            "Kaplan-Meier estimate of hazard function lambda(r)",
            "Hanisch estimate of %s",
            "border corrected estimate of %s",
            "uncorrected estimate of %s")[avail]
  df <- df[, nama[avail]]
  Z <- fv(df, "r", substitute(CDF(r), NULL), bestest, . ~ r, alim, labl, desc,
          fname="CDF")
  fvnames(Z, ".") <- nama[iscdf & avail]
  return(Z)
}

# simple interface for students and code development

compileCDF <- function(D, B, r, ..., han.denom=NULL, check=TRUE) {
  han <- !is.null(han.denom)
  breaks <- breakpts.from.r(r)
  if(check) {
    stopifnot(length(D) == length(B) && all(D >= 0) && all(B >= 0))
    if(han)
      stopifnot(length(han.denom) == length(r))
  }
  D <- as.vector(D)
  B <- as.vector(B)
  # observed (censored) lifetimes
  o <- pmin.int(D, B)
  # censoring indicators
  d <- (D <= B)
  # go
  result <- censtimeCDFest(o, B, d, breaks,
                           HAN=han, 
                           han.denom=han.denom,
                           RAW=TRUE, tt=D)
  result <- rebadge.fv(result, new.fname="compileCDF")
}
