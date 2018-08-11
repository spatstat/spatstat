#
#
#    triplet.family.R
#
#    $Revision: 1.1 $	$Date: 2011/11/05 07:18:51 $
#
#    Family of `third-order' point process models
#
#    triplet.family:      object of class 'isf' 
#	
#
# -------------------------------------------------------------------
#	

triplet.family <-
  list(
       name  = "triplet",
       print = function(self) {
         cat("Family of third-order interactions\n")
       },
       plot = NULL,
       # ----------------------------------------------------
       eval  = function(X,U,EqualPairs,pot,pars,correction, ...) {
  #
  # This is the eval function for the `triplet' family.
  # 
  # This internal function is not meant to be called by the user.
  # It is called by mpl.prepare() during execution of ppm().
  #         
  # The eval functions perform all the manipulations that are common to
  # a given class of interactions. 
  #
  # This function is currently modelled on 'inforder.family'.
  # It simply invokes the potential 'pot' directly
  # and expects 'pot' to return the values of the sufficient statistic S(u,X).
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
  #   potpars     auxiliary parameters for pairpot        list(......)
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
  # In this case the potential function 'pot' should have arguments
  #    pot(X, U, EqualPairs, pars, correction, ...)
  #         
  # It must return a vector with length equal to the number of points in U,
  # or a matrix with as many rows as there are points in U.

         if(!is.ppp(U))
           U <- ppp(U$x, U$y, window=X$window)
         
         POT <- pot(X, U, EqualPairs, pars, correction, ...)

         if(is.matrix(POT)) {
           if(nrow(POT) != U$n)
             stop("Internal error: the potential returned a matrix with the wrong number of rows")
         } else if(is.array(POT) && length(dim(POT)) > 2)
           stop("Internal error: the potential returned an array with more than 2 dimensions")
         else if(is.vector(POT)) {
           if(length(POT) != U$n)
             stop("Internal error: the potential returned a vector with the wrong length")
           POT <- matrix(POT, ncol=1)
         } else
         stop("Internal error: the return value from the potential is not understood")

         return(POT)
       },
######### end of function $eval
       suffstat = NULL
######### end of function $suffstat
)
######### end of list

class(triplet.family) <- "isf"


