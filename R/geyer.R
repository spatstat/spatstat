#
#
#    geyer.S
#
#    $Revision: 2.42 $	$Date: 2018/03/15 07:37:41 $
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
         hasInf   = FALSE,
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                      if(!is.numeric(sat) || length(sat) != 1 || sat < 0)
                       stop("saturation parameter sat must be a positive number")
                    },
         update = NULL, # default OK
         print = NULL,    # default OK
         plot = function(fint, ..., d=NULL, plotit=TRUE) {
           verifyclass(fint, "fii")
           inter <- fint$interaction
           unitz <- unitname(fint)
           if(!identical(inter$name, "Geyer saturation process"))
             stop("Tried to plot the wrong kind of interaction")
           #' fitted interaction coefficient
           theta <- fint$coefs[fint$Vnames]
           #' interaction radius
           r <- inter$par$r
           sat <- inter$par$sat
           xlim <- resolve.1.default(list(xlim=c(0, 1.25 * r)), list(...)) 
           rmax <- max(xlim, d)
           if(is.null(d)) {
             d <- seq(from=0, to=rmax, length.out=1024)
           } else {
             stopifnot(is.numeric(d) &&
                       all(is.finite(d)) &&
                       all(diff(d) > 0))
           }
           #' compute interaction between two points at distance d
           y <- exp(theta * sat * (d <= r))
           #' compute `fv' object
           fun <- fv(data.frame(r=d, h=y, one=1),
                     "r", substitute(h(r), NULL), "h", cbind(h,one) ~ r,
                     xlim, c("r", "h(r)", "1"),
                     c("distance argument r",
                       "maximal interaction h(r)",
                       "reference value 1"),
                     unitname=unitz)
           if(plotit)
             do.call(plot.fv,
                     resolve.defaults(list(fun),
                                      list(...),
                                      list(ylim=range(0,1,y))))
           return(invisible(fun))
         },
         #' end of function 'plot'
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1L])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=dround(gamma)))
         },
         valid = function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1L])
           sat <- self$par$sat
           return(is.finite(loggamma) && (is.finite(sat) || loggamma <= 0))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(!is.na(coeffs))) {
             loggamma <- coeffs[1L]
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
                         splitInf=FALSE,
                         ..., halfway=FALSE, check=TRUE) {
         #' fast evaluator for Geyer interaction
         dont.complain.about(splitInf)
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
           missingdata <- rep.int(TRUE, nX)
         } else {
           Xused <- EqualPairs[,1L]
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
         satcounts <- pmin.int(sat, counts)
         satcounts <- matrix(satcounts, ncol=1)
         if(halfway) {
           # trapdoor used by suffstat()
           answer <- satcounts
         } else if(sat == Inf) {
           # no saturation: fast code
           answer <- 2 * satcounts
         } else {
           # extract counts for data points
           Uindex <- EqualPairs[,2L]
           Xindex <- EqualPairs[,1L]
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
       },
       delta2 = function(X,inte,correction, ..., sparseOK=FALSE) {
         # Sufficient statistic for second order conditional intensity
         # h(X[i] | X) - h(X[i] | X[-j])
         # Geyer interaction
         if(correction == "periodic") return(NULL)
         r   <- inte$par$r
         sat <- inte$par$sat
         result <- geyerdelta2(X, r, sat,
                               correction=correction, sparseOK=sparseOK)
         return(result)
       }
  )
  class(BlankGeyer) <- "interact"
  
  Geyer <- function(r, sat) {
    instantiate.interact(BlankGeyer, list(r = r, sat=sat))
  }

  Geyer <- intermaker(Geyer, BlankGeyer)
  
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
  Uindex <- EqualPairs[,2L]
  Xindex <- EqualPairs[,1L]
  Xsortindex <- rankX[Xindex]
  Usortindex <- rankU[Uindex]
  Cmap <- rep.int(-1L, nU)
  Cmap[Usortindex] <- Xsortindex - 1L
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
           PACKAGE = "spatstat")
  result <- zz$result[rankU]
  return(result)
}

geyerdelta2 <- local({

  geyerdelta2 <- function(X, r, sat, ..., sparseOK=FALSE, correction="none") {
    ## Sufficient statistic for second order conditional intensity
    ## Geyer model
    stopifnot(is.numeric(sat) && length(sat) == 1 && sat > 0)
    ## X could be a ppp or quad.
    if(is.ppp(X)) {
      # evaluate \Delta_{x_i} \Delta_{x_j} S(x) for data points x_i, x_j
      # i.e.  h(X[i]|X) - h(X[i]|X[-j]) where h is first order cif statistic
      return(geydelppp(X, r, sat, correction, sparseOK))
    } else if(is.quad(X)) {
      # evaluate \Delta_{u_i} \Delta_{u_j} S(x) for quadrature points u_i, u_j
      return(geydelquad(X, r, sat, correction, sparseOK))
    } else stop("Internal error: X should be a ppp or quad object")
  }

  geydelppp <- function(X, r, sat, correction, sparseOK) {
    # initialise
    nX <- npoints(X)
    result <- if(!sparseOK) matrix(0, nX, nX) else
              sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                           dims=c(nX, nX))
    ## identify all r-close pairs (ordered pairs i ~ j)
    unweighted <- correction %in% c("border", "none")
    if(unweighted) {
      a <- closepairs(X, r, what="indices")
      I <- a$i
      J <- a$j
      IJ <- cbind(I,J)
      ## count number of r-neighbours for each point
      tvals <- table(factor(I, levels=1:nX))
    } else {
      a <- weightedclosepairs(X, r, correction=correction, what="indices")
      I <- a$i
      J <- a$j
      IJ <- cbind(I,J)
      wIJ <- a$weight
      wmatrix <- sparseMatrix(i=I,j=J,x=wIJ, dims=c(nX,nX))
      wJI <- wmatrix[cbind(J,I)]
      ## total edge-correction weight over r-neighbours for each point
      tvals <- tapplysum(wIJ, list(factor(I, levels=1:nX)))
    }
    # Compute direct part
    # (arising when i~j) 
    tI <- tvals[I]
    tJ <- tvals[J]
    if(unweighted) {
      result[IJ] <-
        pmin(sat, tI) - pmin(sat, tI - 1L) + pmin(sat, tJ) - pmin(sat, tJ - 1L)
    } else {
      result[IJ] <-
        pmin(sat, tI) - pmin(sat, tI - wIJ) +
        pmin(sat, tJ) - pmin(sat, tJ - wJI)
    }
    # Compute indirect part
    # (arising when i~k and j~k for another point k)
    # First find all such triples 
    ord <- (I < J)
    vees <- edges2vees(I[ord], J[ord], nX)
    # evaluate contribution of (k, i, j)
    II <- vees$j
    JJ <- vees$k
    KK <- vees$i
    tKK <- tvals[KK]
    if(unweighted) {
      contribKK <- pmin(sat, tKK) - 2 * pmin(sat, tKK-1L) + pmin(sat, tKK-2L)
    } else {
      wKKII <- wmatrix[cbind(KK, II)]
      wKKJJ <- wmatrix[cbind(KK, JJ)]
      contribKK <- (
        pmin(sat, tKK)
        - pmin(sat, tKK-wKKII)
        - pmin(sat, tKK-wKKJJ)
        + pmin(sat, tKK-wKKII-wKKJJ)
      )
    }
    # for each (i, j), sum the contributions over k 
    if(!sparseOK) {
      II <- factor(II, levels=1:nX)
      JJ <- factor(JJ, levels=1:nX)
      # was:
      # delta3 <- tapply(contribKK, list(I=II, J=JJ), sum)
      # delta3[is.na(delta3)] <- 0
      delta3 <- tapplysum(contribKK, list(I=II, J=JJ))
    } else {
      delta3 <- sparseMatrix(i=II, j=JJ, x=contribKK, dims=c(nX, nX))
    }
    # symmetrise and combine
    result <- result + delta3 + t(delta3)
    return(result)
  }

  geydelquad <- function(Q, r, sat, correction, sparseOK) {
    Z <- is.data(Q)
    U <- union.quad(Q)
    nU <- npoints(U)
    nX <- npoints(Q$data)
    result <- if(!sparseOK) matrix(0, nU, nU) else
              sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                           dims=c(nU, nU))
    ## identify all r-close pairs U[i], U[j]
    unweighted <- correction %in% c("none", "border")
    if(unweighted) {
      a <- closepairs(U, r, what="indices")
      I <- a$i
      J <- a$j
      IJ <- cbind(I, J)
    } else {
      a <- weightedclosepairs(U, r, correction=correction, what="indices")
      I <- a$i
      J <- a$j
      IJ <- cbind(I, J)
      wIJ <- a$weight
      wmatrix <- sparseMatrix(i=I, j=J, x=wIJ, dims=c(nU,nU))
      wJI <- wmatrix[cbind(J,I)]
    }
    ## tag which ones are data points
    zI <- Z[I]
    zJ <- Z[J]
    # count t(U[i], X)
    IzJ <- I[zJ]
    JzJ <- J[zJ]
    if(unweighted) {
      tvals <- table(factor(IzJ, levels=1:nU))
    } else {
      wzJ <- wIJ[zJ]
      tvals <- tapplysum(wzJ, list(factor(IzJ, levels=1:nU)))
    }
    ## Compute direct part
    ## (arising when U[i]~U[j]) 
    tI <- tvals[I]
    tJ <- tvals[J]
    if(unweighted) {
      tIJ <- tI - zJ
      tJI <- tJ - zI
      result[IJ] <- pmin(sat, tIJ + 1L) - pmin(sat, tIJ) +
                    pmin(sat, tJI + 1L) - pmin(sat, tJI)
    } else {
      tIJ <- tI - zJ * wIJ
      tJI <- tJ - zI * wJI
      result[IJ] <- pmin(sat, tIJ + wIJ) - pmin(sat, tIJ) +
                    pmin(sat, tJI + wJI) - pmin(sat, tJI)
    }
    # Compute indirect part
    # (arising when U[i]~X[k] and U[j]~X[k] for another point X[k])
    # First find all such triples
    # Group close pairs X[k] ~ U[j] by index k
    spl <- split(IzJ, factor(JzJ, levels=1:nX))
    grlen <- lengths(spl)
    # Assemble list of triples U[i], X[k], U[j]
    # by expanding each pair U[i], X[k]
    JJ <- unlist(spl[JzJ])
    II <- rep(IzJ, grlen[JzJ])
    KK <- rep(JzJ, grlen[JzJ])
    # remove identical pairs i = j
    ok <- II != JJ
    II <- II[ok]
    JJ <- JJ[ok]
    KK <- KK[ok]
    # evaluate contribution of each triple
    tKK <- tvals[KK]
    zII <- Z[II]
    zJJ <- Z[JJ]
    if(unweighted) {
      tKIJ <- tKK - zII - zJJ 
      contribKK <-
        pmin(sat, tKIJ + 2L) - 2 * pmin(sat, tKIJ + 1L) + pmin(sat, tKIJ)
    } else {
      wKKII <- wmatrix[cbind(KK, II)]
      wKKJJ <- wmatrix[cbind(KK, JJ)]
      tKIJ <- tKK - zII * wKKII - zJJ * wKKJJ
      contribKK <- (pmin(sat, tKIJ + wKKII + wKKJJ)
        - pmin(sat, tKIJ + wKKII)
        - pmin(sat, tKIJ + wKKJJ)
        + pmin(sat, tKIJ))
    }
    # for each (i, j), sum the contributions over k
    if(!sparseOK) {
      II <- factor(II, levels=1:nU)
      JJ <- factor(JJ, levels=1:nU)
      # was:
      # delta4 <- tapply(contribKK, list(I=II, J=JJ), sum)
      # delta4[is.na(delta4)] <- 0
      delta4 <- tapplysum(contribKK, list(I=II, J=JJ))
    } else {
      delta4 <- sparseMatrix(i=II, j=JJ, x=contribKK, dims=c(nU, nU))
    }
    # combine
    result <- result + delta4
    return(result)
  }

  geyerdelta2
})
