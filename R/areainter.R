#
#
#    areainter.R
#
#    $Revision: 1.25 $	$Date: 2013/05/26 13:41:33 $
#
#    The area interaction
#
#    AreaInter()    create an instance of the area-interaction process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#

AreaInter <- local({

  # area-interaction potential function
  areapot <- 
    function(X,U,EqualPairs,pars,correction, ...) {
      uhoh <- !(correction %in% c("border", "none"))
      if(any(uhoh)) {
        nuh <- sum(uhoh)
        warning(paste(ngettext(nuh, "Correction", "Corrections"),
                      commasep(sQuote(correction[uhoh])),
                      ngettext(nuh,
                               "is not supported and was ignored",
                               "are not supported and were ignored")))
      }
      n <- U$n
      answer <- numeric(n)
      r <- pars$r
      if(is.null(r)) stop("internal error: r parameter not found")
      dummies <- !(seq_len(n) %in% EqualPairs[,2])
      if(sum(dummies) > 0)
        answer[dummies] <- -areaGain(U[dummies], X, r)
      ii <- EqualPairs[,1]
      jj <- EqualPairs[,2]
      answer[jj] <- -areaLoss(X, r, subset=ii)
#    for(k in seq_len(nrow(EqualPairs))) {
#      i <- EqualPairs[k,1]
#      j <- EqualPairs[k,2]
#      answer[j] <- -areaGain(U[j], X[-i], r)
#    }
      return(1 + answer/(pi * r^2))
    }

  # template object without family, par, version
  BlankAI <- 
  list(
         name     = "Area-interaction process",
         creator  = "AreaInter",
         family   = "inforder.family", # evaluated later
         pot      = areapot,
         par      = list(r = NULL), # to be filled in
         parnames = "disc radius",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("disc radius r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           logeta <- as.numeric(coeffs[1])
           eta <- exp(logeta)
           return(list(param=list(eta=eta),
                       inames="interaction parameter eta",
                       printable=round(eta,4)))
         },
         valid = function(coeffs, self) {
           eta <- ((self$interpret)(coeffs, self))$param$eta
           return(is.finite(eta))
         },
         project = function(coeffs, self) {
           if((self$valid)(coeffs, self))
             return(NULL)
           return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(is.na(coeffs)))
             return(2 * r)
           logeta <- coeffs[1]
           if(abs(logeta) <= epsilon)
             return(0)
           else
             return(2 * r)
         },
         delta2 = function(X, inte, correction) {
           # Sufficient statistic for second order conditional intensity
           # Area-interaction model 
           # Evaluate \Delta_{x_i} \Delta_{x_j} S(x) for data points x_i, x_j
           # i.e.  h(X[i]|X) - h(X[i]|X[-j])
           #       where h is first order cif statistic
           if(!(correction %in% c("border", "none")))
             return(NULL)
           r <- inte$par$r
           areadelta2(X, r)
         },
       version=NULL # to be added
  )
  class(BlankAI) <- "interact"

  AreaInter <- function(r) {
    instantiate.interact(BlankAI, list(r=r))
  }

  AreaInter
})


areadelta2 <- function(X, r, algorithm=c("C", "nncross", "nnmap")) {
  # Sufficient statistic for second order conditional intensity
  # Area-interaction model 
  # Evaluate \Delta_{x_i} \Delta_{x_j} S(x) for data points x_i, x_j
  # i.e.  h(X[i]|X) - h(X[i]|X[-j])
  #       where h is first order cif statistic
  nX <- npoints(X)
  result <- matrix(0, nX, nX)
  if(nX < 2)
    return(result)
  algorithm <- match.arg(algorithm)
  if(algorithm == "C") {
    # use special purpose C routine
    # called once for each interacting pair of points
    xx <- X$x
    yy <- X$y
    cl <- closepairs(X, 2 * r, what="indices", ordered=FALSE)
    I <- cl$i
    J <- cl$j
    eps <- r/spatstat.options("ngrid.disc")
    for(k in seq_along(I)) {
      i <- I[k]
      j <- J[k]
      # all neighbours of i
      Ki <- union(J[I==i], I[J==i])
      # all neighbours of j
      Kj <- union(J[I==j], I[J==j])
      # relevant neighbours
      K <- setdiff(union(Ki, Kj), c(i,j))
      # call C code
      DUP <- spatstat.options("dupC")
      z <- .C("delta2area",
              xa = as.double(xx[i]),
              ya = as.double(yy[i]),
              xb = as.double(xx[j]),
              yb = as.double(yy[j]),
              nother = as.integer(length(K)),
              xother = as.double(xx[K]),
              yother = as.double(yy[K]),
              radius = as.double(r),
              epsilon = as.double(eps),
              pixcount = as.integer(integer(1)),
              DUP = DUP,
              PACKAGE = "spatstat")
      result[i,j] <- result[j,i] <- z$pixcount
    }
    # normalise
    result <- result * (eps^2)/(pi * r^2)
    return(result)
  }
  # remove any non-interacting points
  relevant <- (nndist(X) <= 2 * r)
  if(!all(relevant)) {
    answer <- matrix(0, nX, nX)
    if(any(relevant)) {
             # call self on subset
      Dok <- areadelta2(X[relevant], r, algorithm)
      answer[relevant,relevant] <- Dok
    }
    return(answer)
  }

  # sort pattern in increasing order of x
  sortX <- (algorithm == "nnmap")
  if(sortX) {
    oX <- fave.order(X$x)
    X <- X[oX]
  }
  
  # area calculation may be restricted to window W for efficiency
  W <- as.owin(X)
  U <- as.rectangle(W)

  # decide pixel resolution
  eps <- r/spatstat.options("ngrid.disc")
  npix <- prod(ceiling(sidelengths(U)/eps))
  if(npix <= 2^20) {
    # do it all in one go
    tile <- list(NULL)
  } else {
    # divide into rectangular tiles
    B <- as.rectangle(W)
    ntile0 <- ceiling(npix/(2^20))
    tile0area <- area.owin(B)/ntile0
    tile0side <- sqrt(tile0area)
    nx <- ceiling(sidelengths(B)[1]/tile0side)
    ny <- ceiling(sidelengths(B)[2]/tile0side)
    tile <- tiles(quadrats(B, nx, ny))
  }
           
  result <- matrix(0, nX, nX)
  for(i in seq_len(length(tile))) {
    # form pixel grid
    Ti <- tile[[i]]
    Wi <- if(is.null(Ti)) W else intersect.owin(W, Ti)
    if(algorithm == "nncross") {
      # Trusted, slow algorithm using nncross
      Z <- as.mask(Wi, eps=eps)
      G <- as.ppp(raster.xy(Z), U, check=FALSE)
      # compute 3 nearest neighbours in X of each grid point
      v <- nncross(G, X, k=1:3)
      # select pixels which have exactly 2 neighbours within distance r
      ok <- with(v, dist.3 > r & dist.2 <= r)
      if(any(ok)) {
        v <- v[ok, , drop=FALSE]
        # accumulate pixel counts -> areas
        counts <- with(v, table(i=factor(which.1, levels=1:nX),
                                j=factor(which.2, levels=1:nX)))
        pixarea <- with(Z, xstep * ystep)
        result <- result + pixarea * (counts + t(counts))
      }
    } else {
      # Faster algorithm using nnmap
      # compute 3 nearest neighbours in X of each grid point
      stuff <- nnmap(X, k=1:3, W=Wi, eps=eps,
                     is.sorted.X=TRUE, sortby="x",
                     outputarray=TRUE)
      dist.2 <- stuff$dist[2,,]
      dist.3 <- stuff$dist[3,,]
      which.1 <- stuff$which[1,,]
      which.2 <- stuff$which[2,,]
      ok <- (dist.3 > r & dist.2 <= r)
      if(any(ok)) {
        which.1 <- as.vector(which.1[ok])
        which.2 <- as.vector(which.2[ok])
        counts <- table(i=factor(which.1, levels=1:nX),
                        j=factor(which.2, levels=1:nX))
        pixarea <- attr(stuff, "pixarea")
        result <- result + pixarea * (counts + t(counts))
      }
    }
  }
  if(sortX) {
    # map back to original ordering
    result[oX, oX] <- result
  }
  # normalise
  result <- result/(pi * r^2)
  return(result)
}

             
