#
#
#    pairpiece.S
#
#    $Revision: 1.23 $	$Date: 2018/03/15 07:37:41 $
#
#    A pairwise interaction process with piecewise constant potential
#
#    PairPiece()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

PairPiece <- local({

  # .... auxiliary functions ........
  delP <- function(i, r) {
    r <- r[-i]
    nr <- length(r)
    if(nr == 0) return(Poisson())
    if(nr == 1) return(Strauss(r))
    return(PairPiece(r))
  }

  # ..... template ..........

  BlankPairPiece <- 
  list(
         name     = "Piecewise constant pairwise interaction process",
         creator  = "PairPiece",
         family   = "pairwise.family", # evaluated later
         pot      = function(d, par) {
                       r <- par$r
                       nr <- length(r)
                       out <- array(FALSE, dim=c(dim(d), nr))
                       out[,,1] <-  (d < r[1])
                       if(nr > 1) {
                         for(i in 2:nr) 
                           out[,,i] <- (d >= r[i-1]) & (d < r[i])
                       }
                       out
                     },
         par      = list(r = NULL), # filled in later
         parnames = "interaction thresholds",
         hasInf = FALSE,
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || !all(r > 0))
                       stop("interaction thresholds r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction thresholds r must be strictly increasing")
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
                       printable=dround(gn)))
         },
        valid = function(coeffs, self) {
           # interaction parameters gamma
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           if(!all(is.finite(gamma))) return(FALSE)
           return(all(gamma <= 1) || gamma[1] == 0)
        },
        project = function(coeffs, self){
           # interaction parameters gamma
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # interaction thresholds r[i]
           r <- self$par$r
           # check for NA or Inf
           bad <- !is.finite(gamma)
           # gamma > 1 forbidden unless hard core
           ishard <- is.finite(gamma[1]) && (gamma[1] == 0)
           if(!ishard)
             bad <- bad | (gamma > 1)
           if(!any(bad))
             return(NULL)
           if(spatstat.options("project.fast") || sum(bad) == 1) {
             # remove smallest threshold with an unidentifiable parameter
             firstbad <- min(which(bad))
             return(delP(firstbad, r))
           } else {
             # consider all candidate submodels
             subs <- lapply(which(bad), delP, r=r)
             return(subs)
           }
        },
        irange = function(self, coeffs=NA, epsilon=0, ...) {
          r <- self$par$r
          if(all(is.na(coeffs)))
            return(max(r))
          gamma <- (self$interpret)(coeffs, self)$param$gammas
          gamma[is.na(gamma)] <- 1
          active <- (abs(log(gamma)) > epsilon)
          if(!any(active))
            return(0)
          else return(max(r[active]))
        },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         r     <- self$par$r
         gamma <- (self$interpret)(coeffs, self)$param$gammas
         # areas of annuli between r[i-1], r[i]
         areas <- pi * diff(c(0,r)^2)
         return(sum(areas * (1-gamma)))
       },
       version=NULL # filled in later
       )
  class(BlankPairPiece) <- "interact"

  PairPiece <- function(r) {
    instantiate.interact(BlankPairPiece, list(r=r))
  }

  PairPiece <- intermaker(PairPiece, BlankPairPiece)
  
  PairPiece
})

                   
