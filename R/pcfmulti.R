#
#   pcfmulti.R
#
#   $Revision: 1.8 $   $Date: 2016/09/21 07:28:58 $
#
#   multitype pair correlation functions
#

pcfcross <- 
  function(X, i, j, ...,
         r=NULL, kernel="epanechnikov", bw=NULL, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         divisor=c("r","d"))
{
  verifyclass(X, "ppp")
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL
  divisor <- match.arg(divisor)
  ##
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  ##
  result <- pcfmulti(X, I, J, ...,
                     r=r, 
                     kernel=kernel, bw=bw, stoyan=stoyan,
                     correction=correction,
                     divisor=divisor,
                     Iname=Iname, Jname=Jname)
  ##
  iname <- make.parseable(paste(i))
  jname <- make.parseable(paste(j))
  result <-
    rebadge.fv(result,
               substitute(g[i,j](r),
                          list(i=iname,j=jname)),
               c("g", paste0("list", paren(paste(iname, jname, sep=",")))),
               new.yexp=substitute(g[list(i,j)](r),
                                   list(i=iname,j=jname)))
  return(result)
}

pcfdot <- 
function(X, i, ...,
         r=NULL, kernel="epanechnikov", bw=NULL, stoyan=0.15,
         correction = c("isotropic", "Ripley", "translate"),
         divisor=c("r", "d"))
{
  verifyclass(X, "ppp")
  stopifnot(is.multitype(X))
  if(missing(correction))
    correction <- NULL
  divisor <- match.arg(divisor)

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- "points"
	
  result <- pcfmulti(X, I, J, ...,
                     r=r, kernel=kernel, bw=bw, stoyan=stoyan,
                     correction=correction,
                     divisor=divisor,
                     Iname=Iname, Jname=Jname)

  iname <- make.parseable(paste(i))
  result <-
    rebadge.fv(result,
               substitute(g[i ~ dot](r), list(i=iname)),
               c("g", paste0(iname, "~symbol(\"\\267\")")),
               new.yexp=substitute(g[i ~ symbol("\267")](r),
                 list(i=iname)))
  return(result)
}


pcfmulti <- function(X, I, J, ...,
                     r=NULL, 
                     kernel="epanechnikov", bw=NULL, stoyan=0.15,
                     correction=c("translate", "Ripley"),
                     divisor=c("r","d"),
                     Iname="points satisfying condition I",
                     Jname="points satisfying condition J")
{
  verifyclass(X, "ppp")
#  r.override <- !is.null(r)
  divisor <- match.arg(divisor)

  win <- X$window
  areaW <- area(win)
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
  
  ## .......... indices I and J .............................
  
  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if(is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")

  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  XI <- X[I]
  XJ <- X[J]

#  lambdaI <- nI/areaW
  lambdaJ <- nJ/areaW
  nIJ <- sum(I & J)
  lambdaIJarea <- (nI * nJ - nIJ)/areaW
  
  ## ...........  kernel bandwidth and support .........................
  
  if(is.null(bw) && kernel=="epanechnikov") {
    # Stoyan & Stoyan 1995, eq (15.16), page 285
    h <- stoyan /sqrt(lambdaJ)
    hmax <- h
    # conversion to standard deviation
    bw <- h/sqrt(5)
  } else if(is.numeric(bw)) {
    # standard deviation of kernel specified
    # upper bound on half-width
    hmax <- 3 * bw
  } else {
    # data-dependent bandwidth selection: guess upper bound on half-width
    hmax <- 2 * stoyan /sqrt(lambdaJ)
  }


########## r values ############################
  # handle argument r 

  rmaxdefault <- rmax.rule("K", win, lambdaJ)
  breaks <- handle.r.b.args(r, NULL, win, rmaxdefault=rmaxdefault)
  if(!(breaks$even))
    stop("r values must be evenly spaced")
  # extract r values
  r <- breaks$r
  rmax <- breaks$max
  # recommended range of r values for plotting
  alim <- c(0, min(rmax, rmaxdefault))

  # initialise fv object
  
  df <- data.frame(r=r, theo=rep.int(1,length(r)))
  fname <- c("g", "list(I,J)")
  yexp <- quote(g[list(I,J)](r))
  out <- fv(df, "r",
            quote(g[I,J](r)), "theo", ,
            alim,
            c("r", makefvlabel(NULL, NULL, fname, "Pois")),
            c("distance argument r", "theoretical Poisson %s"),
            fname=fname,
            yexp=yexp)
  
  ########## smoothing parameters for pcf ############################  
  # arguments for 'density'

  denargs <- resolve.defaults(list(kernel=kernel, bw=bw),
                              list(...),
                              list(n=length(r), from=0, to=rmax))
  
  #################################################
  
  ## compute pairwise distances
  
  ## identify close pairs of points
  what <- if(any(correction == "translate")) "all" else "ijd"
  close <- crosspairs(XI, XJ, rmax+hmax, what=what)
  ## map (i,j) to original serial numbers in X
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
  ## eliminate any identical pairs
  if(nIJ > 0) {
    ok <- (iX != jX)
    if(!all(ok))
      close <- as.list(as.data.frame(close)[ok, , drop=FALSE])
  }
  ## extract information for these pairs (relative to orderings of XI, XJ)
  dclose <- close$d
  icloseI  <- close$i
#  jcloseJ  <- close$j

  ###### compute #######

  if(any(correction=="translate")) {
    # translation correction
    edgewt <- edge.Trans(dx=close$dx, dy=close$dy, W=win, paired=TRUE)
    gT <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)$g
    out <- bind.fv(out,
                   data.frame(trans=gT),
                   makefvlabel(NULL, "hat", fname, "Trans"),
                   "translation-corrected estimate of %s",
                   "trans")
  }
  if(any(correction=="isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
    gR <- sewpcf(dclose, edgewt, denargs, lambdaIJarea, divisor)$g
    out <- bind.fv(out,
                   data.frame(iso=gR),
                   makefvlabel(NULL, "hat", fname, "Ripley"),
                   "isotropic-corrected estimate of %s",
                   "iso")
  }
  
  ## sanity check
  if(is.null(out)) {
    warning("Nothing computed - no edge corrections chosen")
    return(NULL)
  }
  
  # which corrections have been computed?
  corrxns <- rev(setdiff(names(out), "r"))

  # default is to display them all
  formula(out) <- . ~ r
  fvnames(out, ".") <- corrxns

  # 
  unitname(out) <- unitname(X)
  return(out)
}

