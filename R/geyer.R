#
#
#    geyer.S
#
#    $Revision: 2.22 $	$Date: 2012/11/06 08:15:58 $
#
#    Geyer's saturation process
#
#    Geyer()    create an instance of Geyer's saturation process
#                 [an object of class 'interact']
#
#	

Geyer <- local({

  # .......... template ..........

  BlankGeyer <- 
  list(
         name     = "Geyer saturation process",
         creator  = "Geyer",
         family   = "pairsat.family",  # evaluated later
         pot      = function(d, par) {
                         (d <= par$r)  # same as for Strauss
                    },
         par      = list(r = NULL, sat=NULL),  # filled in later
         parnames = c("interaction distance","saturation parameter"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                      if(!is.numeric(sat) || length(sat) != 1 || sat < 1)
                       stop("saturation parameter sat must be a number >= 1")
                      if(ceiling(sat) != floor(sat))
                        warning(paste("saturation parameter sat",
                                      "has a non-integer value"))
                    },
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         valid = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           sat <- self$par$sat
           return(is.finite(loggamma) && (is.finite(sat) || loggamma <= 0))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(!is.na(coeffs))) {
             loggamma <- coeffs[1]
             if(!is.na(loggamma) && (abs(loggamma) <= epsilon))
               return(0)
           }
           return(2 * r)
         },
       version=NULL, # evaluated later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         ..., halfway=FALSE, check=TRUE) {
         # fast evaluator for Geyer interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Geyer")
         r   <- potpars$r
         sat <- potpars$sat
         # first ensure all data points are in U
         nX <- npoints(X)
         nU <- npoints(U)
         Xseq  <- seq_len(nX)
         if(length(EqualPairs) == 0) {
           # no data points currently included 
           missingdata <- rep(TRUE, nX)
         } else {
           Xused <- EqualPairs[,1]
           missingdata <- !(Xseq %in% Xused)
         }
         somemissing <- any(missingdata)
         if(somemissing) {
           # add the missing data points
           nmiss <- sum(missingdata)
           U <- superimpose(U, X[missingdata], W=X$window, check=check)
           # correspondingly augment the list of equal pairs
           originalrows <- seq_len(nU)
           newXindex <- Xseq[missingdata]
           newUindex <- nU + seq_len(nmiss)
           EqualPairs <- rbind(EqualPairs, cbind(newXindex, newUindex))
           nU <- nU + nmiss
         }
         # determine saturated pair counts
         counts <- strausscounts(U, X, r, EqualPairs) 
         satcounts <- pmin(sat, counts)
         satcounts <- matrix(satcounts, ncol=1)
         if(halfway) {
           # trapdoor used by suffstat()
           answer <- satcounts
         } else if(sat == Inf) {
           # no saturation: fast code
           answer <- 2 * satcounts
         } else {
           # extract counts for data points
           Uindex <- EqualPairs[,2]
           Xindex <- EqualPairs[,1]
           Xcounts <- integer(npoints(X))
           Xcounts[Xindex] <- counts[Uindex]
           # evaluate change in saturated counts of other data points
           change <- geyercounts(U, X, r, sat, Xcounts, EqualPairs)
           answer <- satcounts + change
           answer <- matrix(answer, ncol=1)
         }
         if(somemissing)
           answer <- answer[originalrows, , drop=FALSE]
         return(answer)
       }
  )
  class(BlankGeyer) <- "interact"
  
  Geyer <- function(r, sat) {
    instantiate.interact(BlankGeyer, list(r = r, sat=sat))
  }

  Geyer
})

  # ........... externally visible auxiliary functions .........
  
  geyercounts <- function(U, X, r, sat, Xcounts, EqualPairs) {
    # evaluate effect of adding dummy point or deleting data point
    # on saturated counts of other data points
    stopifnot(is.numeric(r))
    stopifnot(is.numeric(sat))
    # for C calls we need finite numbers
    stopifnot(is.finite(r))
    stopifnot(is.finite(sat))
    # sort in increasing order of x coordinate
    oX <- fave.order(X$x)
    oU <- fave.order(U$x)
    Xsort <- X[oX]
    Usort <- U[oU]
    nX <- npoints(X)
    nU <- npoints(U)
    Xcountsort <- Xcounts[oX]
    # inverse: data point i has sorted position i' = rankX[i]
    rankX <- integer(nX)
    rankX[oX] <- seq_len(nX)
    rankU <- integer(nU)
    rankU[oU] <- seq_len(nU)
    # map from quadrature points to data points
    Uindex <- EqualPairs[,2]
    Xindex <- EqualPairs[,1]
    Xsortindex <- rankX[Xindex]
    Usortindex <- rankU[Uindex]
    Cmap <- rep(-1, nU)
    Cmap[Usortindex] <- Xsortindex - 1
    # call C routine
    zz <- .C("Egeyer",
             nnquad = as.integer(nU),
             xquad  = as.double(Usort$x),
             yquad  = as.double(Usort$y),
             quadtodata = as.integer(Cmap),
             nndata = as.integer(nX),
             xdata  = as.double(Xsort$x),
             ydata  = as.double(Xsort$y),
             tdata  = as.integer(Xcountsort),
             rrmax  = as.double(r),
             ssat   = as.double(sat),
             result = as.double(numeric(nU)),
             PACKAGE="spatstat")
    result <- zz$result[rankU]
    return(result)
  }

