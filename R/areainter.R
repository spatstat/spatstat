#
#
#    areainter.R
#
#    $Revision: 1.30 $	$Date: 2014/10/24 00:22:30 $
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

  # area-interaction conditional intensity potential
  #     corresponds to potential -C(x) = n(x) - A(x)/\pi r^2
  areapot <- 
    function(X,U,EqualPairs,pars,correction, ..., W=as.owin(X)) {
      uhoh <- !(correction %in% c("border", "none"))
      if(any(uhoh)) {
        nuh <- sum(uhoh)
        warning(paste(ngettext(nuh, "Correction", "Corrections"),
                      commasep(sQuote(correction[uhoh])),
                      ngettext(nuh,
                               "is not supported and was ignored",
                               "are not supported and were ignored")))
      }
      r <- pars$r
      if(is.null(r)) stop("internal error: r parameter not found")
      n <- U$n
      areas <- numeric(n)
      dummies <- !(seq_len(n) %in% EqualPairs[,2])
      if(sum(dummies) > 0)
        areas[dummies] <- areaGain(U[dummies], X, r, W=W)
      ii <- EqualPairs[,1]
      jj <- EqualPairs[,2]
      areas[jj] <- areaLoss(X, r, subset=ii, W=W)
      return(1 - areas/(pi * r^2))
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


areadelta2 <- local({

  areadelta2 <- function(X, r, ...) {
    # Sufficient statistic for second order conditional intensity
    # Area-interaction model 
    if(is.ppp(X)) return(areadelppp(X, r, ...)) else
    if(inherits(X, "quad")) return(areadelquad(X, r)) else
    stop("internal error: X should be a ppp or quad object")
  }

  areadelppp <- function(X, r, algorithm=c("C", "nncross", "nnmap")) {
    # Evaluate \Delta_{x_i} \Delta_{x_j} S(x) for data points x_i, x_j
    # i.e.  h(X[i]|X) - h(X[i]|X[-j])
    #       where h is first order cif statistic
    algorithm <- match.arg(algorithm)
    nX <- npoints(X)
    result <- matrix(0, nX, nX)
    if(nX < 2)
      return(result)
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
                pixcount = as.integer(integer(1)))
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
        Dok <- areadelppp(X[relevant], r, algorithm)
        answer[relevant,relevant] <- Dok
      }
      return(answer)
    }

    # .............. algorithm using interpreted code ...........
    
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
      tile0area <- area(B)/ntile0
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
        G <- as.ppp(rasterxy.mask(Z), U, check=FALSE)
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

  areadelquad <- function(Q, D, r) {
    # Sufficient statistic for second order conditional intensity
    # Area-interaction model 
    # Evaluate \Delta_{u_j} \Delta_{u_i} S(x) for quadrature points 
    # answer is area(b(u[i],r) \cap b(u[j],r)\setminus \bigcup_k b(x[k],r))
    # where k ranges over all indices that are not equivalent to u[i,j]
    U <- union.quad(Q)
    Z <- is.data(Q)
    nU <- npoints(U)
    xx <- U$x
    yy <- U$y
    # identify all close pairs of quadrature points
    cl <- closepairs(U, 2 * r, what="indices")
    I <- b$i
    J <- b$j
    # find neighbours in X of each quadrature point
    zJ <- Z[J]
    neigh <- split(J[zJ], factor(I[zJ], levels=1:nU))
    # 
    result <- matrix(0, nU, nU)
    eps <- r/spatstat.options("ngrid.disc")
    #
    for(k in seq_along(I)) {
      i <- I[k]
      j <- J[k]
      # all points of X close to U[i]
      Ki <- neigh[[i]]
      # all points of X close to U[j]
      Kj <- neigh[[j]]
      # relevant neighbours
      K <- setdiff(union(Ki, Kj), c(i,j))
      # call C code
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
            pixcount = as.integer(integer(1)))
      result[i,j] <- z$pixcount
    }
    # normalise
    result <- result * (eps^2)/(pi * r^2)
    return(result)
  }

  areadelta2
})
