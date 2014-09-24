#
#	breakpts.S
#
#	A simple class definition for the specification
#       of histogram breakpoints in the special form we need them.
#
#	even.breaks()
#
#	$Revision: 1.14 $	$Date: 2014/09/24 05:23:36 $
#
#
#       Other functions in this directory use the standard Splus function
#	hist() to compute histograms of distance values.
#       One argument of hist() is the vector 'breaks'
#	of breakpoints for the histogram cells. 
#
#       The breakpoints must
#            (a) span the range of the data
#            (b) be given in increasing order
#            (c) satisfy breaks[2] = 0,
#
#	The function make.even.breaks() will create suitable breakpoints.
#
#       Condition (c) means that the first histogram cell has
#       *right* endpoint equal to 0.
#
#       Since all our distance values are nonnegative, the effect of (c) is
#       that the first histogram cell counts the distance values which are
#       exactly equal to 0. Hence F(0), the probability P{X = 0},
#       is estimated without a discretisation bias.
#
#	We assume the histograms have followed the default counting rule
#	in hist(), which is such that the k-th entry of the histogram
#	counts the number of data values in 
#		I_k = ( breaks[k],breaks[k+1] ]	for k > 1
#		I_1 = [ breaks[1],breaks[2]   ]
#
#	The implementations of estimators of c.d.f's in this directory
#       produce vectors of length = length(breaks)-1
#       with value[k] = estimate of F(breaks[k+1]),
#       i.e. value[k] is an estimate of the c.d.f. at the RIGHT endpoint
#       of the kth histogram cell.
#
#       An object of class 'breakpts' contains:
#
#              $val     the actual breakpoints
#              $max     the maximum value (= last breakpoint)
#              $ncells  total number of histogram cells
#              $r       right endpoints, r = val[-1]
#              $even    logical = TRUE if cells known to be evenly spaced
#              $npos    number of histogram cells on the positive halfline
#                        = length(val) - 2,
#                       or NULL if cells not evenly spaced
#              $step    histogram cell width
#                       or NULL if cells not evenly spaced
#       
# --------------------------------------------------------------------
breakpts <- function(val, maxi, even=FALSE, npos=NULL, step=NULL) {
  out <- list(val=val, max=maxi, ncells=length(val)-1, r = val[-1],
              even=even, npos=npos, step=step)
  class(out) <- "breakpts"
  out
}

scalardilate.breakpts <- function(X, f, ...) {
  out <- with(X,
              list(val    = f*val,
                   max    = f*max,
                   ncells = ncells,
                   r      = f*r,
                   even   = even,
                   npos   = npos,
                   step   = f*step))
  class(out) <- "breakpts"
  out
}  
                            
"make.even.breaks" <- 
function(bmax, npos, bstep) {
  if(bmax <= 0)
    stop("bmax must be positive")
  if(missing(bstep) && missing(npos))
    stop(paste("Must specify either", sQuote("bstep"),
               "or", sQuote("npos")))
  if(!missing(npos)) {
    bstep <- bmax/npos
    val <- seq(from=0, to=bmax, length.out=npos+1)
    val <- c(-bstep,val)
    right <- bmax
  } else {
    npos <- ceiling(bmax/bstep)
    right <- bstep * npos
    val <- seq(from=0, to=right, length.out=npos+1)
    val <- c(-bstep,val)
  }
  breakpts(val, right, TRUE, npos, bstep)
}

"as.breakpts" <- function(...) {

  XL <- list(...)

  if(length(XL) == 1) {
    # single argument
    X <- XL[[1]]

    if(!is.null(class(X)) && class(X) == "breakpts")
    # X already in correct form
      return(X)
  
    if(is.vector(X) && length(X) > 2) {
    # it's a vector
      if(X[2] != 0)
        stop("breakpoints do not satisfy breaks[2] = 0")
      # The following test for equal spacing is used in hist.default
      steps <- diff(X)
      if(diff(range(steps)) < 1e-07 * mean(steps))
        # equally spaced
        return(breakpts(X, max(X), TRUE, length(X)-2, steps[1]))
      else
        # unknown spacing
        return(breakpts(X, max(X), FALSE))
    }
  } else {

    # There are multiple arguments.
  
    # exactly two arguments - interpret as even.breaks()
    if(length(XL) == 2)
      return(make.even.breaks(XL[[1]], XL[[2]]))

    # two arguments 'max' and 'npos'
  
    if(!is.null(XL$max) && !is.null(XL$npos))
      return(make.even.breaks(XL$max, XL$npos))

    # otherwise
    stop("Don't know how to convert these data to breakpoints")
  }
  # never reached
}


check.hist.lengths <- function(hist, breaks) {
  verifyclass(breaks, "breakpts")
  nh <- length(hist)
  nb <- breaks$ncells
  if(nh != nb)
    stop(paste("Length of histogram =", nh,
               "not equal to number of histogram cells =", nb))
}

breakpts.from.r <- function(r) {
  if(!is.numeric(r) && !is.vector(r))
    stop("r must be a numeric vector")
  if(length(r) < 2)
    stop(paste("r has length", length(r), "- must be at least 2"))
  if(r[1] != 0)
    stop("First r value must be 0")
  if(any(diff(r) <= 0))
    stop("successive values of r must be increasing")
  dr <- r[2] - r[1]
  b <- c(-dr, r)
  return(as.breakpts(b))
}

handle.r.b.args <- function(r=NULL, breaks=NULL, window, pixeps=NULL,
                            rmaxdefault=NULL) {

        if(!is.null(r) && !is.null(breaks))
          stop(paste("Do not specify both",
                     sQuote("r"), "and", sQuote("breaks")))
  
        if(!is.null(breaks)) {
          breaks <- as.breakpts(breaks)
        } else if(!is.null(r)) {
          breaks <- breakpts.from.r(r)
	} else {
          rmax <- rmaxdefault %orifnull% diameter(Frame(window))
          if(is.null(pixeps)) {
            pixeps <- if(is.mask(window))
                      min(window$xstep, window$ystep) else rmax/128
          }
          rstep <- pixeps/4
          breaks <- make.even.breaks(rmax, bstep=rstep)
        }

        return(breaks)
}
