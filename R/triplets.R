#
#
#    triplets.R
#
#    $Revision: 1.15 $	$Date: 2015/10/21 09:06:57 $
#
#    The triplets interaction
#
#    Triplets()    create an instance of the triplets process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#

Triplets <- local({

  DebugTriplets <- FALSE
  
  # define triplet potential
  TripletPotential <- function(X,U,EqualPairs,pars,correction, ...) {
    if(!all(ok <- correction %in% c("border", "none"))) {
      nbad <- sum(bad <- !ok)
      warning(paste(ngettext(nbad, "Correction", "Corrections"),
                    commasep(sQuote(correction[bad])),
                    ngettext(nbad,
                             "is unavailable and was ignored",
                             "are unavailable and were ignored")))
    }
    # check that all points of X are included in U
    nX <- npoints(X)
    nU <- npoints(U)
    XinU <- if(length(EqualPairs) == 0) integer(0) else EqualPairs[,1]
    missX <- which(table(factor(XinU, levels=1:nX)) == 0)
    if((nmiss <- length(missX)) > 0) {
      # add missing points to (the end of) U
      U <- superimpose(U, X[missX], W=as.owin(X), check=FALSE)
      EqualPairs <- rbind(EqualPairs, cbind(missX, nU + 1:nmiss))
      nU <- nU + nmiss
    }
    iXX <- EqualPairs[,1]
    iXU <- EqualPairs[,2]
    # construct map from X index to U index 
    mapXU <- integer(nX)
    mapXU[iXX] <- iXU
    # construct map from U index to X index 
    mapUX <- rep.int(NA_integer_, nU)
    mapUX[iXU] <- iXX
    # logical vector identifying which quadrature points are in X
    isdata <- rep.int(FALSE, nU)
    isdata[iXU] <- TRUE
    # identify all close pairs u, x
    r <- pars$r
    cp <- crosspairs(U, X, r, what="indices")
    if(DebugTriplets)
      cat(paste("crosspairs at distance", r, "yields", length(cp$i), "pairs\n"))
    IU <- cp$i
    J <- cp$j
    # map X index to U index
    JU <- mapXU[J]
    # Each (Xi, Xj) pair will appear twice - eliminate duplicates
    dupX <- isdata[IU] & isdata[JU] & (IU > JU)
    retain <- !dupX
    IU <- IU[retain]
    JU <- JU[retain]
    if(DebugTriplets)
      cat(paste(sum(dupX), "duplicate pairs removed\n"))
    # find all triangles
    tri <- edges2triangles(IU, JU, nU, friendly=isdata)
    if(DebugTriplets)
      cat(paste(nrow(tri), "triangles identified\n"))
    if(nrow(tri) == 0) {
      # there are no triangles; return vector of zeroes
      return(rep.int(0, nU-nmiss))
    }
    # count triangles containing a given quadrature point
    tcount <- apply(tri, 2,
                    function(x, n) { table(factor(x, levels=1:n)) }, n=nU)
    tcount <- .rowSums(tcount, nrow(tcount), ncol(tcount))
    # select triangles consisting only of data points
    triX <- matrix(mapUX[tri], nrow=nrow(tri))
    isX <- apply(!is.na(triX), 1, all)
    triX <- triX[isX, , drop=FALSE]
    #
    if(nrow(triX) > 0) {
      # count triangles of data points containing each given data point
      tXcount <- apply(triX, 2,
                       function(x, n) { table(factor(x, levels=1:n)) }, n=nX)
      tXcount <- .rowSums(tXcount, nrow(tXcount), ncol(tXcount))
    } else {
      # there are no triangles of data points
      tXcount <- rep.int(0, nX)
    }
    #
    answer <- tcount
    answer[iXU] <- tXcount[iXX]
    if(DebugTriplets)
      cat(paste("Max suff stat: data ", max(tXcount),
                ", dummy ", max(tcount[isdata]), "\n", sep=""))
    # truncate to original size
    if(nmiss > 0)
      answer <- answer[-((nU-nmiss+1):nU)]
    return(answer)
  }
  # set up basic 'triplets' object except for family and parameters
  BlankTripletsObject <- 
    list(
         name     = "Triplets process",
         creator  = "Triplets",
         family   = "triplet.family", # evaluated later
         pot      = TripletPotential,
         par      = list(r=NULL), # filled in later
         parnames = "interaction distance",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=dround(gamma)))
         },
         valid = function(coeffs, self) {
           gamma <- ((self$interpret)(coeffs, self))$param$gamma
           return(is.finite(gamma) && (gamma <= 1))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(is.na(coeffs)))
             return(r)
           loggamma <- coeffs[1]
           if(abs(loggamma) <= epsilon)
             return(0)
           else
             return(r)
         },
         version=NULL # to be added
         )
  class(BlankTripletsObject) <- "interact"
  # define Triplets function
  Triplets <- function(r) {
    instantiate.interact(BlankTripletsObject, list(r=r))
  }
  Triplets <- intermaker(Triplets, BlankTripletsObject)
  
  Triplets
})

