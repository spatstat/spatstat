#
#
#    pairwise.family.S
#
#    $Revision: 1.40 $	$Date: 2013/05/23 07:26:11 $
#
#    The pairwise interaction family of point process models
#
#    pairwise.family:      object of class 'isf' defining pairwise interaction
#	
#
# -------------------------------------------------------------------
#	

pairwise.family <-
  list(
       name  = "pairwise",
       print = function(self) {
         cat("Pairwise interaction family\n")
       },
       plot = function(fint, ..., d=NULL, plotit=TRUE) {
         verifyclass(fint, "fii")
         inter <- fint$interaction
         if(is.null(inter) || is.null(inter$family)
            || inter$family$name != "pairwise")
           stop("Tried to plot the wrong kind of interaction")
         # get fitted coefficients of interaction terms
         # and set coefficients of offset terms to 1
         Vnames <- fint$Vnames
         IsOffset <- fint$IsOffset
         coeff <- rep.int(1, length(Vnames))
         names(coeff) <- Vnames
         coeff[!IsOffset] <- fint$coefs[Vnames[!IsOffset]]
         # 
         pairpot <- inter$pot
         potpars <- inter$par
         rmax <- reach(fint, epsilon=1e-3)
         xlim <- list(...)$xlim
         if(is.infinite(rmax)) {
           if(!is.null(xlim))
             rmax <- max(xlim)
           else {
             warning("Reach of interaction is infinite; need xlim to plot it")
             return(invisible(NULL))
           }
         }
         if(is.null(d)) {
           dmax <- 1.25 * rmax
           d <- seq(from=0, to=dmax, length.out=1024)
         } else {
           stopifnot(is.numeric(d) &&
                     all(is.finite(d)) &&
                     all(diff(d) > 0))
           dmax <- max(d)
         }
         if(is.null(xlim))
           xlim <- c(0, dmax)
         types <- potpars$types
         if(is.null(types)) {
           # compute potential function as `fv' object
           dd <- matrix(d, ncol=1)
           p <- pairpot(dd, potpars)
           if(length(dim(p))==2)
             p <- array(p, dim=c(dim(p),1), dimnames=NULL)
           if(dim(p)[3] != length(coeff))
             stop("Dimensions of potential do not match coefficient vector")
           for(k in seq_len(dim(p)[3])) 
             p[,,k] <- multiply.only.finite.entries( p[,,k] , coeff[k] )
           y <- exp(apply(p, c(1,2), sum))
           ylim <- range(0, 1.1, y, finite=TRUE)
           fun <- fv(data.frame(r=d, h=y, one=1),
                     "r", substitute(h(r), NULL), "h", cbind(h,one) ~ r,
                     xlim, c("r", "h(r)", "1"),
                     c("distance argument r",
                       "pairwise interaction term h(r)",
                       "reference value 1"))
           if(plotit)
             do.call("plot.fv",
                     resolve.defaults(list(fun),
                                      list(...),
                                      list(ylab="Pairwise interaction",
                                           xlab="Distance",
                                           ylim=ylim)))
           return(invisible(fun))
         } else{
           # compute each potential and store in `fasp' object
           if(!is.factor(types))
             types <- factor(types, levels=types)
           m <- length(types)
           nd <- length(d)
           dd <- matrix(rep.int(d, m), nrow=nd * m, ncol=m)
           tx <- rep.int(types, rep.int(nd, m))
           ty <- types
           p <- pairpot(dd, tx, ty, potpars)
           if(length(dim(p))==2)
             p <- array(p, dim=c(dim(p),1), dimnames=NULL)
           if(dim(p)[3] != length(coeff))
             stop("Dimensions of potential do not match coefficient vector")
           for(k in seq_len(dim(p)[3]))
             p[,,k] <- multiply.only.finite.entries( p[,,k] , coeff[k] )
           y <- exp(apply(p, c(1,2), sum))
           ylim <- range(0, 1.1, y, finite=TRUE)
           fns <- vector(m^2, mode="list")
           which <- matrix(, m, m)
           for(i in seq_len(m)) {
             for(j in seq_len(m)) {
               # relevant position in matrix
               ijpos <- i + (j-1) * m
               which[i,j] <- ijpos
               # extract values of potential
               yy <- y[tx == types[i], j]
               # make fv object
               fns[[ijpos]] <- fv(data.frame(r=d, h=yy, one=1),
                     "r", substitute(h(r), NULL), "h", cbind(h,one) ~ r,
                     xlim, c("r", "h(r)", "1"),
                     c("distance argument r",
                       "pairwise interaction term h(r)",
                       "reference value 1"))
               #
             }
           }
           funz <- fasp(fns, which=which,
                        formulae=list(cbind(h, one) ~ r),
                        title="Fitted pairwise interactions",
                        rowNames=paste(types), colNames=paste(types))
           if(plotit)
             do.call("plot.fasp",
                     resolve.defaults(list(funz),
                                      list(...),
                                      list(ylim=ylim,
                                           ylab="Pairwise interaction",
                                           xlab="Distance")))
           return(invisible(funz))
         }
       },
       # end of function `plot'
       # ----------------------------------------------------
       eval  = function(X,U,EqualPairs,pairpot,potpars,correction,
           ..., precomputed=NULL, savecomputed=FALSE, pot.only=FALSE) {
  #
  # This is the eval function for the `pairwise' family.
  # 
  # This internal function is not meant to be called by the user.
  # It is called by mpl.prepare() during execution of ppm().
  #         
  # The eval functions perform all the manipulations that are common to
  # a given class of interactions. 
  #
  # For the `pairwise' family of pairwise-interaction processes,
  # this eval function computes the distances between points,
  # invokes 'pairpot' to evaluate the potential between each pair of points,
  # applies edge corrections, and then sums the pair potential terms.
  #
  # ARGUMENTS:
  #   All 'eval' functions have the following arguments 
  #   which are called in sequence (without formal names)
  #   by mpl.prepare():
  #       
  #   X           data point pattern                      'ppp' object
  #   U           points at which to evaluate potential   list(x,y) suffices
  #   EqualPairs  two-column matrix of indices i, j such that X[i] == U[j]
  #               (or NULL, meaning all comparisons are FALSE)
  #   pot         potential function 
  #   potpars     auxiliary parameters for pot            list(......)
  #   correction  edge correction type                    (string)
  #
  # VALUE:
  #    All `eval' functions must return a        
  #    matrix of values of the total potential
  #    induced by the pattern X at each location given in U.
  #    The rows of this matrix correspond to the rows of U (the sample points);
  #    the k columns are the coordinates of the k-dimensional potential.
  #
  ##########################################################################

  # POTENTIAL:
  #
  # The pair potential function 'pairpot' should be either
  #    pairpot(d, par)            [for potentials that don't depend on marks]
  # or
  #    pairpot(d, tx, tu, par)    [for potentials that do depend on mark]
  # where d is a matrix of interpoint distances,
  # tx is the vector of types for the data points,
  # tu is the vector of types for all quadrature points          
  # and
  #  par is a list of parameters for the potential.
  #         
  # It must return a matrix with the same dimensions as d
  # or an array with its first two dimensions the same as the dimensions of d.

fop <- names(formals(pairpot))
if(identical(all.equal(fop, c("d", "par")), TRUE))
  marx <- FALSE
else if(identical(all.equal(fop, c("d", "tx", "tu", "par")), TRUE))
  marx <- TRUE
else 
  stop("Formal arguments of pair potential function are not understood")

## edge correction argument

if(length(correction) > 1)
  stop("Only one edge correction allowed at a time!")

if(!any(correction == c("periodic", "border", "translate", "translation", "isotropic", "Ripley", "none")))
  stop(paste("Unrecognised edge correction", sQuote(correction)))

#### Compute basic data

if(!is.null(precomputed)) {
  # precomputed
  X <- precomputed$X
  U <- precomputed$U
  EqualPairs <- precomputed$E
  M <- precomputed$M
} else {
  U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window
  # Form the matrix of distances
  M <- crossdist(X, U, periodic=(correction=="periodic"))
}

# Evaluate the pairwise potential 

if(!marx) 
  POT <- pairpot(M, potpars)
else
  POT <- pairpot(M, marks(X), marks(U), potpars)

# Determine whether each column of potential is an offset

  IsOffset <- attr(POT, "IsOffset")

# Check errors and special cases

if(!is.matrix(POT) && !is.array(POT)) {
  if(length(POT) == 0 && X$n ==  0) # empty pattern
    POT <- array(POT, dim=c(dim(M),1))
  else
    stop("Pair potential did not return a matrix or array")
}

if(length(dim(POT)) == 1 || any(dim(POT)[1:2] != dim(M))) {
        whinge <- paste(
           "The pair potential function ",short.deparse(substitute(pairpot)),
           "must produce a matrix or array with its first two dimensions\n",
           "the same as the dimensions of its input.\n", sep="")
	stop(whinge)
}

# make it a 3D array
if(length(dim(POT))==2)
        POT <- array(POT, dim=c(dim(POT),1), dimnames=NULL)
                          
if(correction == "translate" || correction == "translation") {
        edgewt <- edge.Trans(X, U)
        # sanity check ("everybody knows there ain't no...")
        if(!is.matrix(edgewt))
          stop("internal error: edge.Trans() did not yield a matrix")
        if(nrow(edgewt) != X$n || ncol(edgewt) != length(U$x))
          stop("internal error: edge weights matrix returned by edge.Trans() has wrong dimensions")
        POT <- c(edgewt) * POT
} else if(correction == "isotropic" || correction == "Ripley") {
        # weights are required for contributions from QUADRATURE points
        edgewt <- t(edge.Ripley(U, t(M), X$window))
        if(!is.matrix(edgewt))
          stop("internal error: edge.Ripley() did not return a matrix")
        if(nrow(edgewt) != X$n || ncol(edgewt) != length(U$x))
          stop("internal error: edge weights matrix returned by edge.Ripley() has wrong dimensions")
        POT <- c(edgewt) * POT
}

# No pair potential term between a point and itself
if(length(EqualPairs) > 0) {
  nplanes <- dim(POT)[3]
  for(k in 1:nplanes)
    POT[cbind(EqualPairs, k)] <- 0
}

# Return just the pair potential?
if(pot.only)
  return(POT)

# Sum the pairwise potentials 

V <- apply(POT, c(2,3), sum)

# attach the original pair potentials
attr(V, "POT") <- POT

# attach the offset identifier
attr(V, "IsOffset") <- IsOffset

# pass computed information out the back door
if(savecomputed)
  attr(V, "computed") <- list(E=EqualPairs, M=M)
return(V)

},
######### end of function $eval
       suffstat = function(model, X=NULL, callstring="pairwise.family$suffstat") {
# for pairwise models only  (possibly nonstationary)
  verifyclass(model, "ppm")
  if(!identical(model$interaction$family$name,"pairwise"))
    stop("Model is not a pairwise interaction process")

  if(is.null(X)) {
    X <- data.ppm(model)
    modelX <- model
  } else {
    verifyclass(X, "ppp")
    modelX <- update(model, X, method="mpl")
  }

  # find data points which do not contribute to pseudolikelihood
  mplsubset <- getglmdata(modelX)$.mpl.SUBSET
  mpldata   <- is.data(quad.ppm(modelX))
  contribute <- mplsubset[mpldata]

  Xin  <- X[contribute]
  Xout <- X[!contribute]
  
  # partial model matrix arising from ordered pairs of data points
  # which both contribute to the pseudolikelihood
  Empty <- X[numeric(0)]
  momINxIN <- partialModelMatrix(Xin, Empty, model, "suffstat")

  # partial model matrix arising from ordered pairs of data points
  # the second of which does not contribute to the pseudolikelihood
  mom <- partialModelMatrix(Xout, Xin, model, "suffstat")
  indx <- Xout$n + (1:(Xin$n))
  momINxOUT <- mom[indx, , drop=FALSE]

  # parameters
  order2  <- names(coef(model)) %in% model$internal$Vnames
  order1  <- !order2

  result <- 0 * coef(model)
  
  if(any(order1)) {
    # first order contributions can be determined from INxIN
    o1terms  <- momINxIN[ , order1, drop=FALSE]
    o1sum   <- apply(o1terms, 2, sum)
    result[order1] <- o1sum
  }
  if(any(order2)) {
    # adjust for double counting of ordered pairs in INxIN but not INxOUT
    o2termsINxIN  <- momINxIN[, order2, drop=FALSE]
    o2termsINxOUT <- momINxOUT[, order2, drop=FALSE]
    o2sum   <- apply(o2termsINxIN, 2, sum)/2 + apply(o2termsINxOUT, 2, sum)
    result[order2] <- o2sum
  }

  return(result)
  },
######### end of function $suffstat
  delta2 = function(X, inte, correction, ...) {
    # Sufficient statistic for second order conditional intensity
    # Evaluate \Delta_{x_i} \Delta_{x_j} S(x) for data points x_i, x_j
    # i.e.  h(X[i]|X) - h(X[i]|X[-j]) where h is first order cif statistic
    nX <- npoints(X)
    if(!(correction %in% c("border", "none"))) {
      # calculation involves edge correction.
      # Use generic evaluator
      E <- cbind(1:nX, 1:nX)
      result <- pairwise.family$eval(X,X,E,inte$pot,inte$par,correction,
                                     pot.only=TRUE)
    } else {
      # No edge correction weights
      R <- reach(inte)
      # identify close pairs (without repetition)
      cl <- closepairs(X, R, ordered=FALSE)
      I <- cl$i
      J <- cl$j
      D <- matrix(cl$d, ncol=1)
      # evaluate potential for these pairs
      POT <- inte$pot(D, inte$par)
      # ensure 3D array
      nd <- length(dim(POT))
      if(nd == 0) {
        POT <- array(POT, dim=c(length(POT), 1, 1))
      } else if(nd == 2) {
        POT <- array(POT, dim=c(dim(POT), 1))
      }
      p <- dim(POT)[3]
      # create result array
      result <- array(0, dim=c(nX, nX, p))
      # insert results
      II <- rep(I, p)
      JJ <- rep(J, p)
      KK <- rep(1:p, each=length(I))
      result[cbind(II,JJ,KK)] <- POT
      result[cbind(JJ,II,KK)] <- POT
    }
    #
    return(result)
  }
######### end of function $delta2
)
######### end of list

class(pairwise.family) <- "isf"


