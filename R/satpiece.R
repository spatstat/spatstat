#
#
#    satpiece.S
#
#    $Revision: 1.14 $	$Date: 2012/01/18 11:04:54 $
#
#    Saturated pairwise interaction process with piecewise constant potential
#
#    SatPiece()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

SatPiece <- local({

  # ..... auxiliary functions ......

  delSP <- function(i, r, sat) {
    r   <- r[-i]
    sat <- sat[-i]
    nr <- length(r)
    if(nr == 0) return(Poisson())
    if(nr == 1) return(Geyer(r, sat))
    return(SatPiece(r, sat))
  }

  # ....... template object ..........
  
  BlankSatPiece <- 
    list(
         name     = "piecewise constant Saturated pairwise interaction process",
         creator  = "SatPiece",
         family   = "pairsat.family", # evaluated later
         pot      = function(d, par) {
                       r <- par$r
                       nr <- length(r)
                       out <- array(FALSE, dim=c(dim(d), nr))
                       out[,,1] <- (d < r[1])
                       if(nr > 1) {
                         for(i in 2:nr) 
                           out[,,i] <- (d >= r[i-1]) & (d < r[i])
                       }
                       out
                    },
         par      = list(r = NULL, sat=NULL), # filled in later
         parnames = c("interaction thresholds", "saturation parameters"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || !all(r > 0))
                        stop("interaction thresholds r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction thresholds r must be strictly increasing")
                      if(!is.numeric(sat) || any(sat < 0))
                        stop("saturation parameters must be nonnegative numbers")
                      if(any(ceiling(sat) != floor(sat)))
                        warning("saturation parameter has a non-integer value")
                      if(length(sat) != length(r) && length(sat) != 1)
                        stop("vectors r and sat must have equal length")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           r <- self$par$r
           npiece <- length(r)
           # extract coefficients
           gammas <- exp(as.numeric(coeffs))
           # name them
           gn <- gammas
           names(gn) <- paste("[", c(0,r[-npiece]),",", r, ")", sep="")
           #
           return(list(param=list(gammas=gammas),
                       inames="interaction parameters gamma_i",
                       printable=round(gn,4)))
         },
        valid = function(coeffs, self) {
           # interaction parameters gamma must be
           #   non-NA 
           #   finite, if sat > 0
           #   less than 1, if sat = Inf
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           sat <- self$par$sat
           if(any(is.na(gamma)))
             return(FALSE)
           return(all((is.finite(gamma) | sat == 0)
                      & (gamma <= 1 | sat != Inf)))
        },
        project = function(coeffs, self){
          loggammas <- as.numeric(coeffs)
          sat <- self$par$sat
          r   <- self$par$r
          ok <- is.finite(loggammas) & (is.finite(sat) | loggammas <= 0)
          if(all(ok))
            return(NULL)
          if(!any(ok))
            return(Poisson())
          bad <- !ok
          if(spatstat.options("project.fast") || sum(bad) == 1) {
            # remove smallest threshold with an unidentifiable parameter
            firstbad <- min(which(bad))
            return(delSP(firstbad, r, sat))
          } else {
            # consider all candidate submodels
            subs <- lapply(which(bad), delSP, r=r, sat=sat)
            return(subs)
          }
        },
        irange = function(self, coeffs=NA, epsilon=0, ...) {
          r <- self$par$r
          sat <- self$par$sat
          if(all(is.na(coeffs)))
            return(2 * max(r))
          gamma <- (self$interpret)(coeffs, self)$param$gammas
          gamma[is.na(gamma)] <- 1
          active <- (abs(log(gamma)) > epsilon) & (sat > 0)
          if(!any(active))
            return(0)
          else return(2 * max(r[active]))
        },
       version=NULL # added later
  )
  class(BlankSatPiece) <- "interact"

  SatPiece <- function(r, sat) {
    instantiate.interact(BlankSatPiece, list(r=r, sat=sat))
  }

  SatPiece
})


                  
