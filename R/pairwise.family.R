#
#
#    pairwise.family.S
#
#    $Revision: 1.71 $	$Date: 2018/04/06 08:55:03 $
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
         unitz <- unitname(fint)
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
                       "reference value 1"),
                     unitname=unitz)
           if(plotit)
             do.call(plot.fv,
                     resolve.defaults(list(fun),
                                      list(...),
                                      list(ylim=ylim)))
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
               fns[[ijpos]] <-
                   fv(data.frame(r=d, h=yy, one=1),
                      "r", substitute(h(r), NULL), "h", cbind(h,one) ~ r,
                      xlim, c("r", "h(r)", "1"),
                      c("distance argument r",
                        "pairwise interaction term h(r)",
                        "reference value 1"),
                      unitname=unitz)
               #
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
                                      list(ylim=ylim)))
           return(invisible(funz))
         }
       },
       # end of function `plot'
       # ----------------------------------------------------
    eval  = function(X,U,EqualPairs,pairpot,potpars,correction,
                     splitInf=FALSE,
                     ..., Reach=NULL,
                     precomputed=NULL, savecomputed=FALSE,
                     pot.only=FALSE) {
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
  # The additional argument 'splitInf' is also permitted.
  #         
  # It must return a matrix with the same dimensions as d
  # or an array with its first two dimensions the same as the dimensions of d.

pt <- PairPotentialType(pairpot) # includes validation of pair potential

## edge correction argument

if(length(correction) > 1)
  stop("Only one edge correction allowed at a time!")

if(!any(correction == c("periodic", "border", "translate", "translation", "isotropic", "Ripley", "none")))
  stop(paste("Unrecognised edge correction", sQuote(correction)))

 no.correction <- 

#### Compute basic data

   # Decide whether to apply faster algorithm using 'closepairs'
   use.closepairs <-
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

if(use.closepairs) {
  POT <- evalPairPotential(X,U,EqualPairs,pairpot,potpars,Reach)
} else {
  POT <- do.call.matched(pairpot,
                         list(d=M,
                              tx=marks(X),
                              tu=marks(U),
                              par=potpars))
}

# Determine whether each component of potential is an offset

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

#' positive case

  if(splitInf) {
    IsNegInf <- (POT == -Inf)
    POT[IsNegInf] <- 0
  }
      

# handle corrections      

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
  for(k in 1:nplanes) {
    POT[cbind(EqualPairs, k)] <- 0
    if(splitInf) IsNegInf[cbind(EqualPairs, k)] <- FALSE
  }
}

# reattach the negative infinity for re-use by special code 
if(splitInf) attr(POT, "IsNegInf") <- IsNegInf
      
# Return just the pair potential?
if(pot.only) 
  return(POT)

# Sum the pairwise potentials over data points for each quadrature point

V <- apply(POT, c(2,3), sum)

# Handle positive case      

if(splitInf) 
  attr(V, "-Inf") <- apply(IsNegInf, 2, any)
      
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
  indx <- Xout$n + seq_len(Xin$n)
  momINxOUT <- mom[indx, , drop=FALSE]

  # parameters
  order2  <- names(coef(model)) %in% model$internal$Vnames
  order1  <- !order2

  result <- 0 * coef(model)
  
  if(any(order1)) {
    # first order contributions can be determined from INxIN
    o1terms  <- momINxIN[ , order1, drop=FALSE]
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
  delta2 = function(X, inte, correction, ..., sparseOK=FALSE) {
    #' Sufficient statistic for second order conditional intensity
    #' for pairwise interaction processes
    #' Equivalent to evaluating pair potential.
    if(is.ppp(X)) {
      seqX <- seq_len(npoints(X))
      E <- cbind(seqX, seqX)
      R <- reach(inte)
      result <- pairwise.family$eval(X,X,E,
                                     inte$pot,inte$par,
                                     correction,
                                     pot.only=TRUE,
                                     Reach=R, splitInf=TRUE)
      M <- attr(result, "IsNegInf")
      if(sparseOK)
        result <- as.sparse3Darray(result)
      if(!is.null(M)) {
        #' validate
        if(length(dim(M)) != 3)
          stop("Internal error: IsNegInf is not a 3D array")
        #' collapse vector-valued potential, yielding a matrix
        M <- apply(M, c(1,2), any)
        if(!is.matrix(M)) M <- matrix(M, nrow=nX)
        #' count conflicts
        hits <- colSums(M)
        #'  hits[j] == 1 implies that X[j] violates hard core with only one X[i]
        #'  and therefore changes status if X[i] is deleted.
        deltaInf <- M
        deltaInf[, hits != 1] <- FALSE
        if(sparseOK)
          deltaInf <- as(deltaInf, "sparseMatrix")
        #' 
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
      result <- pairwise.family$eval(U,U,E,
                                     inte$pot,inte$par,
                                     correction,
                                     pot.only=TRUE,
                                     Reach=R, splitInf=TRUE)
      M <- attr(result, "IsNegInf")
      if(sparseOK)
        result <- as.sparse3Darray(result)
      if(!is.null(M)) {
        #' validate
        if(length(dim(M)) != 3)
          stop("Internal error: IsNegInf is not a 3D array")
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
        #'
        if(sparseOK)
          deltaInf <- as(deltaInf, "sparseMatrix")
        #'
        attr(result, "deltaInf") <- deltaInf
      }
    }
    return(result)
  }
######### end of function $delta2
)
######### end of list

class(pairwise.family) <- "isf"


# externally visible

PairPotentialType <- function(pairpot) {
  stopifnot(is.function(pairpot))
  fop <- names(formals(pairpot))
  v <- match(list(fop),
             list(c("d", "par"),
                  c("d", "tx", "tu", "par")))
  if(is.na(v))
    stop("Formal arguments of pair potential function are not understood",
         call.=FALSE)
  marked <- (v == 2)
  return(list(marked=marked))
}

evalPairPotential <- function(X, P, E, pairpot, potpars, R) {
  # Evaluate pair potential without edge correction weights
  nX <- npoints(X)
  nP <- npoints(P)
  pt <- PairPotentialType(pairpot) # includes validation
  # determine dimension of potential, etc
  fakePOT <- do.call.matched(pairpot,
                             list(d=matrix(, 0, 0),
                                  tx=marks(X)[integer(0)],
                                  tu=marks(P)[integer(0)],
                                  par=potpars))
  IsOffset <- attr(fakePOT, "IsOffset")
  fakePOT <- ensure3Darray(fakePOT)
  Vnames <- dimnames(fakePOT)[[3]]
  p <- dim(fakePOT)[3]
  # Identify close pairs X[i], P[j]
  cl <- crosspairs(X, P, R, what="ijd")
  I <- cl$i
  J <- cl$j
  D <- matrix(cl$d, ncol=1)
  # deal with empty cases
  if(nX == 0 || nP == 0 || length(I) == 0) {
    di <- c(nX, nP, p)
    dn <- list(NULL, NULL, Vnames)
    result <- array(0, dim=di, dimnames=dn)
    attr(result, "IsOffset") <- IsOffset
    return(result)
  }
  # evaluate potential for close pairs
  # POT is a 1-column matrix or array, with rows corresponding to close pairs
  if(!pt$marked) {
    # unmarked
    POT <- do.call.matched(pairpot,
                           list(d=D,
                                par=potpars))
    IsOffset <- attr(POT, "IsOffset")
  } else {
    # marked
    marX <- marks(X)
    marP <- marks(P)
    if(!identical(levels(marX), levels(marP)))
      stop("Internal error: marks of X and P have different levels")
    types <- levels(marX)
    mI <- marX[I]
    mJ <- marP[J]
    POT <- NULL
    # split data by type of P[j]
    for(k in types) {
      relevant <- which(mJ == k)
      if(length(relevant) > 0) {
        fk <- factor(k, levels=types)
        POTk <- do.call.matched(pairpot,
                                list(d=D[relevant,  , drop=FALSE],
                                     tx=mI[relevant],
                                     tu=fk,
                                     par=potpars))
        POTk <- ensure3Darray(POTk)
        if(is.null(POT)) {
          #' use first result of 'pairpot' to determine dimension
          POT <- array(, dim=c(length(I), 1, dim(POTk)[3]))
          #' capture information about offsets, and names of interaction terms
          IsOffset <- attr(POTk, "IsOffset")
          Vnames <- dimnames(POTk)[[3]]
        }
        # insert values just computed
        POT[relevant, , ] <- POTk
      }
    }
  }
  POT <- ensure3Darray(POT)
  p <- dim(POT)[3]
  # create result array
  result <- array(0, dim=c(npoints(X), npoints(P), p),
                  dimnames=list(NULL, NULL, Vnames))
  # insert results
  II <- rep(I, p)
  JJ <- rep(J, p)
  KK <- rep(1:p, each=length(I))
  IJK <- cbind(II, JJ, KK)
  result[IJK] <- POT
  # finally identify identical pairs and set value to 0
  if(length(E) > 0) {
    E.rep <- apply(E, 2, rep, times=p)
    p.rep <- rep(1:p, each=nrow(E))
    result[cbind(E.rep, p.rep)] <- 0
  }
  attr(result, "IsOffset") <- IsOffset
  return(result)
}
