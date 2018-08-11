#
#
#    badgey.S
#
#    $Revision: 1.17 $	$Date: 2018/03/15 07:37:41 $
#
#    Hybrid Geyer process
#
#    BadGey()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

BadGey <- local({

  # ........... auxiliary functions ..............
  delBG <- function(i, r, sat) {
    r   <- r[-i]
    if(length(r) == length(sat)) {
      r   <- r[-i]
      sat <- sat[-i]
    } else if(length(sat) == 1) {
      r <- r[-i]
    } else stop("Mismatch in dimensions of arguments r and sat")
    nr <- length(r)
    if(nr == 0) return(Poisson())
    if(nr == 1) return(Geyer(r, sat))
    return(BadGey(r, sat))
  }

  # .............. template ....................
  
  BlankBG <- 
  list(
         name     = "hybrid Geyer process",
         creator  = "BadGey",
         family   = "pairsat.family",  # will be evaluated later
         pot      = function(d, par) {
                       r <- par$r
                       nr <- length(r)
                       out <- array(FALSE, dim=c(dim(d), nr))
                       for(i in 1:nr) 
                         out[,,i] <- (d <= r[i])
                       out
                    },
         par      = list(r = NULL, sat=NULL), # to fill in later
         parnames = c("interaction radii", "saturation parameters"),
         hasInf   = FALSE,
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || !all(r > 0))
                        stop("interaction radii r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction radii r must be strictly increasing")
                      if(!is.numeric(sat) || any(sat < 0))
                        stop("saturation parameters must be nonnegative numbers")
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
           names(gn) <- paste("[0,", r, ")", sep="")
           #
           return(list(param=list(gammas=gammas),
                       inames="interaction parameters gamma_i",
                       printable=dround(gn)))
         },
        valid = function(coeffs, self) {
           # interaction parameters gamma must be
           #   non-NA 
           #   finite, if sat > 0
           #   less than 1, if sat = Inf
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           sat <- self$par$sat
           if(anyNA(gamma))
             return(FALSE)
           return(all((is.finite(gamma) | sat == 0)
                      & (gamma <= 1 | sat != Inf)))
        },
        project = function(coeffs, self){
          loggammas <- as.numeric(coeffs)
          sat <- self$par$sat
          r   <- self$par$r
          good <- is.finite(loggammas) & (is.finite(sat) | loggammas <= 0)
          if(all(good))
            return(NULL)
          if(!any(good))
            return(Poisson())
          bad <- !good
          if(spatstat.options("project.fast") || sum(bad) == 1) {
            # remove smallest threshold with an unidentifiable parameter
            firstbad <- min(which(bad))
            return(delBG(firstbad, r, sat))
          } else {
            # consider all candidate submodels
            subs <- lapply(which(bad), delBG, r=r, sat=sat)
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
       version=NULL, # to be added later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         ..., halfway=FALSE) {
         # fast evaluator for BadGey interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for BadGey")
         r   <- potpars$r
         sat <- potpars$sat
         # ensure r and sat have equal length
         if(length(r) != length(sat)) {
           if(length(r) == 1)
             r <- rep.int(r, length(sat))
           else if(length(sat) == 1)
             sat <- rep.int(sat, length(r))
           else stop("lengths of r and sat do not match")
         }
         # first ensure all data points are in U
         nX <- npoints(X)
         nU <- npoints(U)
         Xseq  <- seq_len(nX)
         if(length(EqualPairs) == 0) {
           # no data points currently included 
           missingdata <- rep.int(TRUE, nX)
         } else {
           Xused <- EqualPairs[,1L]
           missingdata <- !(Xseq %in% Xused)
         }
         somemissing <- any(missingdata)
         if(somemissing) {
           # add the missing data points
           nmiss <- sum(missingdata)
           U <- superimpose(U, X[missingdata], W=X$window)
           # correspondingly augment the list of equal pairs
           originalrows <- seq_len(nU)
           newXindex <- Xseq[missingdata]
           newUindex <- nU + seq_len(nmiss)
           EqualPairs <- rbind(EqualPairs, cbind(newXindex, newUindex))
           nU <- nU + nmiss
         }
         nterms <- length(r)
         answer <- matrix(, nrow=nU, ncol=nterms)
         for(k in 1:nterms) {
           # first determine saturated pair counts
           counts <- strausscounts(U, X, r[k], EqualPairs) 
           satcounts <- pmin.int(sat[k], counts)
           # trapdoor used by suffstat() 
           if(halfway) 
             answer[,k] <- satcounts
           else if(sat[k] == Inf)
             answer[,k] <- 2 * satcounts
           else {
             # extract counts for data points
             Uindex <- EqualPairs[,2L]
             Xindex <- EqualPairs[,1L]
             Xcounts <- integer(npoints(X))
             Xcounts[Xindex] <- counts[Uindex]
             # evaluate change in saturated counts of other data points
             change <- geyercounts(U, X, r[k], sat[k], Xcounts, EqualPairs)
             answer[,k] <- satcounts + change
           }
         }
         if(somemissing)
           answer <- answer[originalrows, , drop=FALSE]
         return(answer)
       }
  )
  class(BlankBG) <- "interact"

  BadGey <- function(r, sat) {
    instantiate.interact(BlankBG, list(r=r, sat=sat))
  }

  BadGey <- intermaker(BadGey, BlankBG)
  
  BadGey

})


