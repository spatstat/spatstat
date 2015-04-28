#
#
#    ord.family.S
#
#    $Revision: 1.16 $	$Date: 2013/04/25 06:37:43 $
#
#    The Ord model (family of point process models)
#
#    ord.family:      object of class 'isf' defining Ord model structure
#	
#
# -------------------------------------------------------------------
#	

ord.family <-
  list(
         name  = "ord",
         print = function(self) {
                      cat("Ord model family\n")
         },
         eval  = function(X, U, EqualPairs, pot, pars, ...) {
  #
  # This auxiliary function is not meant to be called by the user.
  # It computes the distances between points,
  # evaluates the pair potential and applies edge corrections.
  #
  # Arguments:
  #   X           data point pattern                      'ppp' object
  #   U           points at which to evaluate potential   list(x,y) suffices
  #   EqualPairs  two-column matrix of indices i, j such that X[i] == U[j]
  #               (or NULL, meaning all comparisons are FALSE)
  #   pot         potential function                      function(d, p)
  #   pars        auxiliary parameters for pot            list(......)
  #   ...         IGNORED                             
  #
  # Value:
  #    matrix of values of the potential
  #    induced by the pattern X at each location given in U.
  #    The rows of this matrix correspond to the rows of U (the sample points);
  #    the k columns are the coordinates of the k-dimensional potential.
  #
  # Note:
  # The potential function 'pot' will be called as
  #    pot(M, pars)   where M is a vector of tile areas.
  # It must return a vector of the same length as M
  # or a matrix with number of rows equal to the length of M
  ##########################################################################

nX <- npoints(X)
nU <- length(U$x)       # number of data + dummy points

seqX <- seq_len(nX)
seqU <- seq_len(nU)

# determine which points in the combined list are data points
if(length(EqualPairs) > 0)           
  is.data <- seqU %in% EqualPairs[,2] 
else
  is.data <- rep.int(FALSE, nU)

#############################################################################
# First compute Dirichlet tessellation of data
# and its total potential (which could be vector-valued)
#############################################################################

marks(X) <- NULL
Wdata <- dirichletWeights(X)   # sic - these are the tile areas.
Pdata <- pot(Wdata, pars)
summa <- function(P) {
  if(is.matrix(P))
    matrowsum(P)
  else if(is.vector(P) || length(dim(P))==1 )
    sum(P)
  else
    stop("Don't know how to take row sums of this object")
}
total.data.potential <- summa(Pdata)

# Initialise V

dimpot <- dim(Pdata)[-1]  # dimension of each value of the potential function
                          # (= numeric(0) if potential is a scalar)

dimV <- c(nU, dimpot)
if(length(dimV) == 1)
  dimV <- c(dimV, 1)

V <- array(0, dim=dimV)

rowV <- array(seqU, dim=dimV)

#################### Next, evaluate V for the data points.  ###############
# For each data point, compute Dirichlet tessellation
# of the data with this point removed.
# Compute difference of total potential.
#############################################################################


for(j in seq_len(nX)) {
        #  Dirichlet tessellation of data without point j
  Wminus <- dirichletWeights(X[-j])
        #  regressor is the difference in total potential
  V[rowV == j] <- total.data.potential - summa(pot(Wminus, pars))
}


#################### Next, evaluate V for the dummy points   ################
# For each dummy point, compute Dirichlet tessellation
# of (data points together with this dummy point) only. 
# Take difference of total potential.
#############################################################################

for(j in seqU[!is.data]) {
  Xplus <- superimpose(X, list(x=U$x[j], y=U$y[j]), W=X$window)
  #  compute Dirichlet tessellation (of these points only!)
  Wplus <- dirichletWeights(Xplus)
  #  regressor is difference in total potential
  V[rowV == j] <- summa(pot(Wplus, pars)) - total.data.potential
}

cat("dim(V) = \n")
print(dim(V))

return(V)

} ######### end of function $eval                            

) ######### end of list

class(ord.family) <- "isf"
