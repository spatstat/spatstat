#
#
#    pairsat.family.S
#
#    $Revision: 1.38 $	$Date: 2013/04/25 06:37:43 $
#
#    The saturated pairwise interaction family of point process models
#
#    (an extension of Geyer's saturation process to all pairwise interactions)
#
#    pairsat.family:         object of class 'isf'
#                     defining saturated pairwise interaction
#	
#
# -------------------------------------------------------------------
#	

pairsat.family <-
  list(
         name  = "saturated pairwise",
         print = function(self) {
                      cat("Saturated pairwise interaction family\n")
         },
         eval  = function(X,U,EqualPairs,pairpot,potpars,correction,
                          ..., precomputed=NULL, savecomputed=FALSE,
                               halfway=FALSE) {
  #
  # This is the eval function for the `pairsat' family.
  # 
  # This internal function is not meant to be called by the user.
  # It is called by mpl.prepare() during execution of ppm().
  #         
  # The eval functions perform all the manipulations that are common to
  # a given class of interactions. 
  #
  # For the `pairsat' family of pairwise-interaction processes,
  # this eval function computes the distances between points,
  # invokes 'pairpot' to evaluate the potential between each pair of points,
  # applies edge corrections, and then sums the pair potential terms
  # applying the saturation threshold.
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
  ########################################################################
  #
  # POTENTIAL:
  # The pair potential function 'pairpot' will be called as
  #    pairpot(M, potpars)   where M is a matrix of interpoint distances.
  # It must return a matrix with the same dimensions as M
  # or an array with its first two dimensions the same as the dimensions of M.
  #           
  # NOTE:
  #   Note the Geyer saturation threshold must be given in 'potpars$sat'

  ##########################################################################

           
# coercion should be unnecessary, but this is useful for debugging
X <- as.ppp(X)
U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window

# saturation parameter(s)
saturate <- potpars$sat

if(is.null(saturate)) {
  # pairwise interaction 
  V <- pairwise.family$eval(X, U, EqualPairs,
                            pairpot, potpars, correction, ...,
                            precomputed=precomputed,
                            savecomputed=savecomputed)
  return(V)
}

# first ensure all data points are included in the quadrature points
nX <- npoints(X)
nU <- npoints(U)
Xseq  <- seq_len(nX)
if(length(EqualPairs) == 0) {
  # no data points currently included 
  missingdata <- rep.int(TRUE, nX)
} else {
  Xused <- EqualPairs[,1]
  missingdata <- !(Xseq %in% Xused)
}
somemissing <- any(missingdata)
if(somemissing) {
  # add the missing data points
  originalrows <- seq_len(nU)
  nmiss <- sum(missingdata)
  U <- superimpose(U, X[missingdata], W=X$window)
  # correspondingly augment the list of equal pairs
  newXindex <- Xseq[missingdata]
  newUindex <- nU + seq_len(nmiss)
  EqualPairs <- rbind(EqualPairs, cbind(newXindex, newUindex))
  nU <- nU + nmiss
}

# compute the pair potentials POT and the unsaturated potential sums V

V <- pairwise.family$eval(X, U, EqualPairs, pairpot, potpars, correction, ...)
POT <- attr(V, "POT")

computed <- attr(V, "computed")   # could be NULL

#
# V is a matrix with rows = quadrature points,
#                    columns = coordinates of potential
# POT is an array with rows = data points
#                      columns = quadrature points
#                      planes = coordinates of potential

#################################################################
################## saturation part ##############################
#################################################################

# check dimensions and ensure 'saturate' is a vector
ns <- length(saturate)
np <- ncol(V)
if(ns == 1 && np > 1)
  saturate <- rep.int(saturate, np)
else if(ns != np)
  stop("Length of vector of saturation parameters is incompatible with the pair potential", call.=FALSE)

# replicate as a matrix and as an array
saturate2 <- array(saturate[slice.index(V, 2)], dim=dim(V))
saturate3 <- array(saturate[slice.index(POT, 3)], dim=dim(POT))
#
# (a) compute SATURATED potential sums
V.sat <- pmin(V, saturate2)

if(halfway)
  return(V.sat)
#
# (b) compute effect of addition/deletion of dummy/data point j
# on the UNSATURATED potential sum of each data point i
#
# Identify data points
is.data <- seq_len(npoints(U)) %in% EqualPairs[,2] # logical vector corresp. to rows of V

# Extract potential sums for data points only
V.data <- V[is.data, , drop=FALSE]

# replicate them so that V.dat.rep[i,j,k] = V.data[i, k]
V.dat.rep <- aperm(array(V.data, dim=c(dim(V.data), U$n)), c(1,3,2))

# make a logical array   col.is.data[i,j,k] = is.data[j]
col.is.data <- array(is.data[slice.index(POT, 2)], dim=dim(POT))

# compute value of unsaturated potential sum for each data point i
# obtained after addition/deletion of each dummy/data point j
                                  
V.after <- V.dat.rep + ifelse(col.is.data, -POT, POT)
#
#
# (c) difference of SATURATED potential sums for each data point i
# before & after increment/decrement of each dummy/data point j
#
# saturated values after increment/decrement
V.after.sat <- array(pmin(saturate3, V.after), dim=dim(V.after))
# saturated values before
V.dat.rep.sat <- array(pmin(saturate3, V.dat.rep), dim=dim(V.dat.rep))
# difference
V.delta <- V.after.sat - V.dat.rep.sat
V.delta <- ifelse(col.is.data, -V.delta, V.delta)
#
# (d) Sum (c) over all data points i
V.delta.sum <- apply(V.delta, c(2,3), sum)
#
# (e) Result
V <- V.sat + V.delta.sum

##########################################
# remove rows corresponding to supplementary points
if(somemissing)
      V <- V[originalrows, , drop=FALSE]

### tack on the saved computations from pairwise.family$eval
if(savecomputed)
  attr(V, "computed") <- computed

return(V)

},     ######### end of function $eval                            
suffstat = function(model, X=NULL, callstring="pairsat.family$suffstat") {

# for saturated pairwise models only  (possibly nonstationary)
  verifyclass(model, "ppm")
  if(!identical(model$interaction$family$name,"saturated pairwise"))
    stop("Model is not a saturated pairwise interaction process") 

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
  
  Empty <- X[integer(0)]
  mom <- partialModelMatrix(X, Empty, model, "suffstat", halfway=TRUE)
  # halfway=TRUE is passed to pairsat.family$eval
  # and yields matrix of saturated potential sums 

  # take only those terms that contribute to the pseudolikelihood
  mom <- mom[contribute, , drop=FALSE]
  
  result <- apply(mom, 2, sum)
  return(result)
         

} ######### end of function $suffstat
)     ######### end of list

class(pairsat.family) <- "isf"
