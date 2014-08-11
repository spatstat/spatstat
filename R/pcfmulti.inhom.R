#
#   pcfmulti.inhom.R
#
#   $Revision: 1.12 $   $Date: 2014/04/16 05:31:22 $
#
#   inhomogeneous multitype pair correlation functions
#
#

pcfcross.inhom <- 
  function(X, i, j, lambdaI=NULL, lambdaJ=NULL, ...,
         r=NULL, breaks=NULL,
         kernel="epanechnikov", bw=NULL, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  g <- pcfmulti.inhom(X, I, J, lambdaI, lambdaJ, ...,
                      r=r,breaks=breaks,
                      kernel=kernel, bw=bw, stoyan=stoyan,
                      correction=correction,
                      sigma=sigma, varcov=varcov,
                      Iname=Iname, Jname=Jname)
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(g,
               substitute(g[inhom,i,j](r),
                          list(i=iname,j=jname)),
               c("g", paste0("list", paren(paste("inhom", i, j, sep=",")))),
               new.yexp=substitute(g[list(inhom,i,j)](r),
                                   list(i=iname,j=jname)))
  attr(result, "dangerous") <- attr(g, "dangerous")
  return(result)
}

pcfdot.inhom <- 
function(X, i, lambdaI=NULL, lambdadot=NULL, ...,
         r=NULL, breaks=NULL,
         kernel="epanechnikov", bw=NULL, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         sigma=NULL, varcov=NULL)
{
  verifyclass(X, "ppp")
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points")
	
  g <- pcfmulti.inhom(X, I, J, lambdaI, lambdadot, ...,
                      r=r,breaks=breaks,
                      kernel=kernel, bw=bw, stoyan=stoyan,
                      correction=correction,
                      sigma=sigma, varcov=varcov,
                      Iname=Iname, Jname=Jname)
  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(g,
               substitute(g[inhom, i ~ dot](r), list(i=iname)),
               c("g", paste0("list(inhom,", iname, "~symbol(\"\\267\"))")),
               new.yexp=substitute(g[list(inhom, i ~ symbol("\267"))](r),
                 list(i=iname)))
  if(!is.null(dang <- attr(g, "dangerous"))) {
    dang[dang == "lambdaJ"] <- "lambdadot"
    dang[dang == "lambdaIJ"] <- "lambdaIdot"
    attr(result, "dangerous") <- dang
  }
  return(result)
}


pcfmulti.inhom <- function(X, I, J, lambdaI=NULL, lambdaJ=NULL, ...,
                           r=NULL, breaks=NULL, 
                           kernel="epanechnikov", bw=NULL, stoyan=0.15,
                           correction=c("translate", "Ripley"),
                           sigma=NULL, varcov=NULL,
                           Iname="points satisfying condition I",
                           Jname="points satisfying condition J")
{
  verifyclass(X, "ppp")
  r.override <- !is.null(r)

  win <- X$window
  area <- area.owin(win)
  npts <- npoints(X)
  
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("translate", "Ripley")
  correction <- pickoption("correction", correction,
                           c(isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, win$type, correction.given)
  
  # bandwidth  
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(npts/area)
    hmax <- h
    # conversion to standard deviation
    bw <- h/sqrt(5)
  } else if(is.numeric(bw)) {
    # standard deviation of kernel specified
    # upper bound on half-width
    hmax <- 3 * bw
  } else {
    # data-dependent bandwidth selection: guess upper bound on half-width
    hmax <- 2 * stoyan /sqrt(npts/area)
  }

  ##########  indices I and J  ########################
  
  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != npts || length(J) != npts)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  XI <- X[I]
  XJ <- X[J]
  
  ########## intensity values #########################

  dangerous <- c("lambdaI", "lambdaJ")
  dangerI <- dangerJ <- TRUE
  
  if(is.null(lambdaI)) {
      # Estimate density by leave-one-out kernel smoothing
    dangerI <- FALSE
    lambdaI <- density(XI, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
  } else {
    # lambda values provided
    if(is.vector(lambdaI)) 
      check.nvector(lambdaI, nI)
    else if(is.im(lambdaI)) 
      lambdaI <- safelookup(lambdaI, XI)
    else if(is.function(lambdaI)) 
      lambdaI <- lambdaI(XI$x, XI$y)
    else stop(paste(sQuote("lambdaI"),
                    "should be a vector, a pixel image, or a function"))
  }

  if(is.null(lambdaJ)) {
      # Estimate density by leave-one-out kernel smoothing
    dangerJ <- FALSE
    lambdaJ <- density(XJ, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=TRUE)
  } else {
    # lambda values provided
    if(is.vector(lambdaJ)) 
      check.nvector(lambdaJ, nJ)
    else if(is.im(lambdaJ)) 
      lambdaJ <- safelookup(lambdaJ, XJ)
    else if(is.function(lambdaJ)) 
      lambdaJ <- lambdaJ(XJ$x, XJ$y)
    else stop(paste(sQuote("lambdaJ"),
                    "should be a vector, a pixel image, or a function"))
  }

  danger <- dangerI || dangerJ
  
  ########## r values ############################
  # handle arguments r and breaks 

  rmaxdefault <- rmax.rule("K", win, npts/area)        
  breaks <- handle.r.b.args(r, breaks, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  fname <- c("g", "list(inhom,I,J)")
  out <- fv(df, "r",
            quote(g[inhom,I,J](r)), "theo", ,
            alim,
            c("r", makefvlabel(NULL, NULL, fname, "pois")),            
            c("distance argument r", "theoretical Poisson %s"),
            fname=fname,
            yexp=quote(g[list(inhom,I,J)](r)))
  
  ########## smoothing parameters for pcf ############################  
  # arguments for 'density'

  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=length(r), from=0, to=rmax))
  
  #################################################
  
  # compute pairwise distances
  
# identify close pairs of points
  close <- crosspairs(XI, XJ, rmax+hmax)
# map (i,j) to original serial numbers in X
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
# eliminate any identical pairs
  if(any(I & J)) {
    ok <- (iX != jX)
    if(!all(ok)) {
      close$i  <- close$i[ok]
      close$j  <- close$j[ok]
      close$xi <- close$xi[ok]
      close$yi <- close$yi[ok]
      close$xj <- close$xj[ok]
      close$yj <- close$yj[ok]
      close$dx <- close$dx[ok]
      close$dy <- close$dy[ok]
      close$d  <- close$d[ok]
    }
  }
# extract information for these pairs (relative to orderings of XI, XJ)
  dclose <- close$d
  icloseI  <- close$i
  jcloseJ  <- close$j

# Form weight for each pair
  weight <- 1/(lambdaI[icloseI] * lambdaJ[jcloseJ])

  ###### compute #######

  if(any(correction=="translate")) {
    # translation correction
    edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
    gT <- sewpcf(dclose, edgewt * weight, denargs, area)$g
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   makefvlabel(NULL, "hat", fname, "Trans"),
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
    gR <- sewpcf(dclose, edgewt * weight, denargs, area)$g
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   makefvlabel(NULL, "hat", fname, "Ripley"),
                   "isotropic-corrected estimate of %s",
                   "iso")
  }
  
  # sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # which corrections have been computed?
  nama2 <- names(out)
  corrxns <- rev(nama2[nama2 != "r"])

  # default is to display them all
  formula(out) <- deparse(as.formula(paste(
                       "cbind(",
                        paste(corrxns, collapse=","),
                        ") ~ r")))
  unitname(out) <- unitname(X)

  if(danger)
    attr(out, "dangerous") <- dangerous
  return(out)
}

