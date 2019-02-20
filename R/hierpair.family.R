#
#
#    hierpair.family.R
#
#    $Revision: 1.11 $	$Date: 2019/02/20 03:34:50 $
#
#    The family of hierarchical pairwise interactions
#
#
# -------------------------------------------------------------------
#	

hierpair.family <-
  list(
       name  = "hierpair",
       print = function(self) {
         splat("Hierarchical pairwise interaction family")
       },
       plot = function(fint, ..., d=NULL, plotit=TRUE) {
         verifyclass(fint, "fii")
         inter <- fint$interaction
         if(is.null(inter) || is.null(inter$family)
            || inter$family$name != "hierpair")
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
         if(is.null(types))
           stop("Unable to determine types of points")
         if(!is.factor(types))
           types <- factor(types, levels=types)
         ## compute each potential and store in `fasp' object
         m <- length(types)
         nd <- length(d)
         dd <- matrix(rep.int(d, m), nrow=nd * m, ncol=m)
         tx <- rep.int(types, rep.int(nd, m))
         ty <- types
         p <- pairpot(dd, tx, ty, potpars)
         if(length(dim(p))==2)
           p <- array(p, dim=c(dim(p),1), dimnames=NULL)
         if(dim(p)[3L] != length(coeff))
           stop("Dimensions of potential do not match coefficient vector")
         for(k in seq_len(dim(p)[3L]))
           p[,,k] <- multiply.only.finite.entries( p[,,k] , coeff[k] )
         y <- exp(apply(p, c(1,2), sum))
         ylim <- range(0, 1.1, y, finite=TRUE)
         fns <- vector(m^2, mode="list")
         which <- matrix(, m, m)
         for(i in seq_len(m)) {
           for(j in seq_len(m)) {
             ## relevant position in matrix
             ijpos <- i + (j-1L) * m
             which[i,j] <- ijpos
             ## extract values of potential
             yy <- y[tx == types[i], j]
             ## make fv object
             fns[[ijpos]] <- fv(data.frame(r=d, h=yy, one=1),
                                "r", quote(h(r)), "h", cbind(h,one) ~ r,
                                xlim, c("r", "h(r)", "1"),
                                c("distance argument r",
                                  "pairwise interaction term h(r)",
                                  "reference value 1"))
           }
         }
         funz <- fasp(fns, which=which,
                      formulae=list(cbind(h, one) ~ r),
                      title="Fitted pairwise interactions",
                      rowNames=paste(types), colNames=paste(types))
         if(plotit)
           do.call(plot.fasp,
                   resolve.defaults(list(funz),
                                    list(...),
                                    list(ylim=ylim,
                                         ylab="Pairwise interaction",
                                         xlab="Distance")))
         return(invisible(funz))
       },
       # end of function `plot'
       # ----------------------------------------------------
       eval  = function(X,U,EqualPairs,pairpot,potpars,correction,
           ..., Reach=NULL, precomputed=NULL, savecomputed=FALSE,
           pot.only=FALSE) {
         ##
         ## This is the eval function for the `hierpair' family.
         ## 

fop <- names(formals(pairpot))
if(isTRUE(all.equal(fop, c("d", "par"))))
  marx <- FALSE
else if(isTRUE(all.equal(fop, c("d", "tx", "tu", "par"))))
  marx <- TRUE
else 
  stop("Formal arguments of pair potential function are not understood")

## edge correction argument

if(length(correction) > 1)
  stop("Only one edge correction allowed at a time!")

if(!any(correction == c("periodic", "border", "translate", "translation", "isotropic", "Ripley", "none")))
  stop(paste("Unrecognised edge correction", sQuote(correction)))

 no.correction <- 

#### Compute basic data

   # Decide whether to apply faster algorithm using 'closepairs'
   use.closepairs <- FALSE &&
     (correction %in% c("none", "border", "translate", "translation")) &&
     !is.null(Reach) && is.finite(Reach) &&
     is.null(precomputed) && !savecomputed 

if(!is.null(precomputed)) {
  # precomputed
  X <- precomputed$X
  U <- precomputed$U
  EqualPairs <- precomputed$E
  M <- precomputed$M
} else {
  U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window
  if(!use.closepairs) 
    # Form the matrix of distances
    M <- crossdist(X, U, periodic=(correction=="periodic"))
}

nX <- npoints(X)
nU <- npoints(U)
dimM <- c(nX, nU)

# Evaluate the pairwise potential without edge correction

if(use.closepairs)
  POT <- evalPairPotential(X,U,EqualPairs,pairpot,potpars,Reach)
else if(!marx) 
  POT <- pairpot(M, potpars)
else
  POT <- pairpot(M, marks(X), marks(U), potpars)

# Determine whether each column of potential is an offset

  IsOffset <- attr(POT, "IsOffset")

# Check errors and special cases

if(!is.matrix(POT) && !is.array(POT)) {
  if(length(POT) == 0 && X$n ==  0) # empty pattern
    POT <- array(POT, dim=c(dimM,1))
  else
    stop("Pair potential did not return a matrix or array")
}

if(length(dim(POT)) == 1 || any(dim(POT)[1:2] != dimM)) {
        whinge <- paste0(
           "The pair potential function ",short.deparse(substitute(pairpot)),
           " must produce a matrix or array with its first two dimensions\n",
           "the same as the dimensions of its input.\n")
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
  nplanes <- dim(POT)[3L]
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
       suffstat = function(model, X=NULL, callstring="hierpair.family$suffstat") {
# for hierarchical pairwise models only  (possibly nonstationary)
  verifyclass(model, "ppm")
  if(!identical(model$interaction$family$name,"hierpair"))
    stop("Model is not a hierarchical pairwise interaction process")

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
  Empty <- X[integer(0)]
  momINxIN <- partialModelMatrix(Xin, Empty, model, "suffstat")

  # partial model matrix at data points which contribute to the pseudolikelihood
  momIN <-
    partialModelMatrix(X, Empty, model, "suffstat")[contribute, , drop=FALSE]
  
  # partial model matrix arising from ordered pairs of data points
  # the second of which does not contribute to the pseudolikelihood
  mom <- partialModelMatrix(Xout, Xin, model, "suffstat")
  indx <- Xout$n + seq_len(Xin$n)
  momINxOUT <- mom[indx, , drop=FALSE]

  ## determine which canonical covariates are true second-order terms
  ## eg 'mark1x1' 
  typ <- levels(marks(X))
  vn <- paste0("mark", typ, "x", typ)
  order2  <- names(coef(model)) %in% vn
  order1  <- !order2

  result <- 0 * coef(model)
  
  if(any(order1)) {
    # first order contributions (including 'mark1x2' etc)
    o1terms  <- momIN[ , order1, drop=FALSE]
    o1sum   <- colSums(o1terms)
    result[order1] <- o1sum
  }
  if(any(order2)) {
    # adjust for double counting of ordered pairs in INxIN but not INxOUT
    o2termsINxIN  <- momINxIN[, order2, drop=FALSE]
    o2termsINxOUT <- momINxOUT[, order2, drop=FALSE]
    o2sum   <- colSums(o2termsINxIN)/2 + colSums(o2termsINxOUT)
    result[order2] <- o2sum
  }

  return(result)
  },
######### end of function $suffstat
  delta2 = function(X, inte, correction, ...) {
  # Sufficient statistic for second order conditional intensity
  # for hierarchical pairwise interaction processes
  # Equivalent to evaluating pair potential.
    if(is.ppp(X)) {
      seqX <- seq_len(npoints(X))
      E <- cbind(seqX, seqX)
      R <- reach(inte)
      POT <- hierpair.family$eval(X,X,E,
                                  inte$pot,inte$par,
                                  correction,
                                  pot.only=TRUE,
                                  Reach=R, splitInf=TRUE)
      result <- aperm(POT, c(2,1,3))
      M <- attr(POT, "IsNegInf")
      if(!is.null(M)) {
        #' validate
        if(length(dim(M)) != 3)
          stop("Internal error: IsNegInf is not a 3D array")
        M <- aperm(M, c(2,1,3))
        #' collapse vector-valued potential, yielding a matrix
        M <- apply(M, c(1,2), any)
        if(!is.matrix(M)) M <- matrix(M, nrow=nX)
        #' count conflicts
        hits <- colSums(M)
        #'  hits[j] == 1 implies that X[j] violates hard core with only one X[i]
        #'  and therefore changes status if X[i] is deleted.
        deltaInf <- M
        deltaInf[, hits != 1] <- FALSE
        attr(result, "deltaInf") <- deltaInf
      }
    } else if(is.quad(X)) {
      U <- union.quad(X)
      izdat <- is.data(X)
      nU <- npoints(U)
      nX <- npoints(X$data)
      seqU <- seq_len(nU)
      E <- cbind(seqU, seqU)
      R <- reach(inte)
      POT <- hierpair.family$eval(U,U,E,
                                  inte$pot,inte$par,
                                  correction,
                                  pot.only=TRUE,
                                  Reach=R, splitInf=TRUE)
      result <- aperm(POT, c(2,1,3))
      M <- attr(POT, "IsNegInf")
      if(!is.null(M)) {
        #' validate
        if(length(dim(M)) != 3)
          stop("Internal error: IsNegInf is not a 3D array")
        M <- aperm(M, c(2,1,3))
        #' consider conflicts with data points
        MXU <- M[izdat, , , drop=FALSE]
        #' collapse vector-valued potential, yielding a matrix
        MXU <- apply(MXU, c(1,2), any)
        if(!is.matrix(MXU)) MXU <- matrix(MXU, nrow=nX)
        #' count data points conflicting with each quadrature point
        nhitdata <- colSums(MXU)
        #' for a conflicting pair U[i], U[j],
        #' status of U[j] will change when U[i] is added/deleted
        #' iff EITHER
        #'     U[i] = X[i] is a data point and
        #'     U[j] is only in conflict with X[i],
        deltaInf <- apply(M, c(1,2), any)
        deltaInf[izdat, nhitdata != 1] <- FALSE
        #' OR
        #'     U[i] is a dummy point,
        #'     U[j] has no conflicts with X.
        deltaInf[!izdat, nhitdata != 0] <- FALSE
        attr(result, "deltaInf") <- deltaInf
      }
    }
    return(result)
  }
######### end of function $delta2
)
######### end of list

class(hierpair.family) <- "isf"

