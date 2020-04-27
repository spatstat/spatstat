#'
#'    density.lpp.R
#'
#'    Method for 'density' for lpp objects
#'
#'    Copyright (C) 2017-2020 Greg McSwiggan and Adrian Baddeley
#'

density.lpp <- function(x, sigma=NULL, ...,
                        weights=NULL,
                        distance=c("path", "euclidean"),
                        continuous=TRUE,
                        kernel="gaussian") {
  stopifnot(inherits(x, "lpp"))
  distance <- match.arg(distance)

  if(distance == "euclidean") {
    #' Euclidean 2D kernel
    return(densityQuick.lpp(x, sigma, ..., kernel=kernel, weights=weights))
  }

  #' kernel is 1-D
  kernel <- match.kernel(kernel)
  if(continuous && (kernel == "gaussian")) {
    #' equal-split continuous with Gaussian kernel: use heat equation
    return(densityHeat(x, sigma, ..., weights=weights))
  }

  ##' Okabe-Sugihara equal-split method
  return(densityEqualSplit(x, sigma, ..., kernel=kernel, weights=weights))
}

density.splitppx <- function(x, sigma=NULL, ...) {
  if(!all(sapply(x, is.lpp)))
    stop("Only implemented for patterns on a linear network")
  solapply(x, density.lpp, sigma=sigma, ...)
}

densityEqualSplit <- function(x, sigma=NULL, ...,
                              weights=NULL,
                              kernel="epanechnikov",
                              continuous=TRUE,
                              epsilon=1e-6,
                              verbose=TRUE, debug=FALSE, savehistory=TRUE) {
  ## Based on original code by Adrian Baddeley and Greg McSwiggan 2014-2016
  check.1.real(sigma)
  if(bandwidth.is.infinite(sigma)) {
    out <- as.linim(flatdensityfunlpp(x, weights=weights))
    attr(out, "sigma") <- sigma
    return(out)
  }
  L <- as.linnet(x)
  # weights
  np <- npoints(x)
  if(is.null(weights)) {
    weights <- rep(1, np)
  } else {
    stopifnot(is.numeric(weights))
    check.nvector(weights, np, oneok=TRUE)
    if(length(weights) == 1L) weights <- rep(weights, np) 
  }
  # pixellate linear network
  Llines <- as.psp(L)
  linemask <- as.mask.psp(Llines, ...)
  lineimage <- as.im(linemask, value=0)
  # extract pixel centres
  xx <- raster.x(linemask)
  yy <- raster.y(linemask)
  mm <- linemask$m
  xx <- as.vector(xx[mm])
  yy <- as.vector(yy[mm])
  pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=FALSE)
  pixdf <- data.frame(xc=xx, yc=yy)
  # project pixel centres onto lines
  p2s <- project2segment(pixelcentres, Llines)
  projloc <- as.data.frame(p2s$Xproj)
  projmap <- as.data.frame(p2s[c("mapXY", "tp")])
  projdata <- cbind(pixdf, projloc, projmap)
  # initialise pixel values
  values <- rep(0, nrow(pixdf))
  # Extract local coordinates of data
  n <- npoints(x)
  coo <- coords(x)
  seg <- coo$seg
  tp  <- coo$tp
  # lengths of network segments
  Llengths <- lengths_psp(Llines)
  # initialise stack
  stack <- data.frame(seg=integer(0), from=logical(0), 
                      distance=numeric(0), weight=numeric(0),
                      generation=integer(0))
  # process each data point
  for(i in seq_len(n)) {
    segi <- seg[i]
    tpi  <- tp[i]
    len <- Llengths[segi]
    # evaluate kernel on segment containing x[i]
    relevant <- (projmap$mapXY == segi)
    values[relevant] <- values[relevant] +
      dkernel(len * (projmap$tp[relevant] - tpi),
              kernel=kernel, sd=sigma)
    # push the two tails onto the stack
    stack <- rbind(data.frame(seg = c(segi, segi),
                              from  = c(TRUE, FALSE), 
                              distance = len * c(tpi, 1-tpi),
                              weight = rep(weights[i], 2L),
                              generation = rep(1L, 2)),
                   stack)
  }
  Lfrom <- L$from
  Lto   <- L$to
  if(verbose)
    niter <- 0
  if(savehistory)
    history <- data.frame(iter=integer(0), qlen=integer(0),
                          totmass=numeric(0), maxmass=numeric(0))

  lastgen <- resolve.1.default(list(lastgen=Inf), list(...))
  sortgen <- resolve.1.default(list(sortgen=FALSE), list(...))
  sortgen <- sortgen || is.finite(lastgen) 

  ## process the stack
  while(nrow(stack) > 0) {
    if(debug) print(stack)
    masses <- with(stack, abs(weight) * pkernel(distance,
                                                kernel=kernel,
                                                sd=sigma,
                                                lower.tail=FALSE))
    totmass <- sum(masses)
    maxmass <- max(masses)
    if(savehistory)
      history <- rbind(history,
                       data.frame(iter=nrow(history)+1L,
                                  qlen=nrow(stack),
                                  totmass=totmass,
                                  maxmass=maxmass))
    if(verbose) {
      niter <- niter + 1L
      cat(paste("Iteration", niter, "\tStack length", nrow(stack), "\n"))
      cat(paste("Total stack mass", totmass, "\tMaximum", maxmass, "\n"))
    }
    # trim
    tiny <- (masses < epsilon)
    if(any(tiny)) {
      if(verbose) {
        ntiny <- sum(tiny)
        cat(paste("Removing", ntiny,
                  "tiny", ngettext(ntiny, "tail", "tails"), "\n"))
      }
      stack <- stack[!tiny, ]
    }
    if(nrow(stack) == 0)
      break;
    # pop the top of the stack
    H  <- stack[1L, , drop=FALSE]
    stack <- stack[-1L, , drop=FALSE]
    # segment and vertex
    Hseg <- H$seg
    Hvert <- if(H$from) Lfrom[Hseg] else Lto[Hseg]
    Hdist <- H$distance
    Hgen <- H$generation
    ## finished processing?
    if(Hgen > lastgen)
      break;
    # find all segments incident to this vertex
    incident <- which((Lfrom == Hvert) | (Lto == Hvert))
    degree <- length(incident)
    # exclude reflecting paths?
    if(!continuous)
      incident <- setdiff(incident, Hseg)
    for(J in incident) {
      lenJ <- Llengths[J]
      # determine whether Hvert is the 'to' or 'from' endpoint of segment J
      H.is.from <- (Lfrom[J] == Hvert)
      # update weight
      if(continuous) {
        Jweight <- H$weight * (2/degree - (J == Hseg))
      } else {
        Jweight <- H$weight/(degree-1)
      }
      # increment density on segment
      relevant <- (projmap$mapXY == J)
      tp.rel <- projmap$tp[relevant]
      d.rel <- lenJ * (if(H.is.from) tp.rel else (1 - tp.rel))
      values[relevant] <- values[relevant] +
        Jweight * dkernel(d.rel + Hdist, kernel=kernel, sd=sigma)
      # push other end of segment onto stack
      stack <- rbind(data.frame(seg = J,
                                from  = !(H.is.from),
                                distance = lenJ + Hdist,
                                weight = Jweight,
                                generation = Hgen + 1L),
                     stack)
      if(sortgen)
        stack <- stack[order(stack$generation), , drop=FALSE]
      print(stack)
    }
  }
  # attach values to nearest pixels
  Z <- lineimage
  Z[pixelcentres] <- values
  # attach exact line position data
  df <- cbind(projdata, values)
  out <- linim(L, Z, df=df)
  attr(out, "sigma") <- sigma
  if(savehistory)
    attr(out, "history") <- history
  return(out)
}

densityHeat <- function(x, sigma, ...,
                        at=c("pixels", "points"),
                        leaveoneout=TRUE, weights=NULL, 
                        dx=NULL, dt=NULL, iterMax=1e6, verbose=FALSE) {
  stopifnot(is.lpp(x))
  check.1.real(sigma)
  at <- match.arg(at)
  if(!is.null(weights)) 
    check.nvector(weights, npoints(x))

  if(bandwidth.is.infinite(sigma)) {
    out <- switch(at,
                  pixels = as.linim(flatdensityfunlpp(x, weights=weights)),
                  points = flatdensityatpointslpp(x, weights=weights,
                                                  leaveoneout=leaveoneout))
    attr(out, "sigma") <- sigma
    return(out)
  }
  if(at == "points") {
    out <- densitypointsLPP(x, sigma, ..., leaveoneout=leaveoneout, 
                            weights=weights, nsigma=1,
                            dx=dx, dt=dt, iterMax=iterMax,
                            verbose=verbose)
    attr(out, "sigma") <- sigma
    return(out)
  }
  ## internal arguments
  fun         <- resolve.1.default(list(fun=FALSE), list(...))
  finespacing <- resolve.1.default(list(finespacing=FALSE), list(...))
  ## 
  ## determine algorithm parameters
  L <- as.linnet(x)
  p <- resolve.heat.steps(sigma, dx=dx, dt=dt, iterMax=iterMax, L=L, ...,
                          verbose=verbose)
  ## go
  a <- FDMKERNEL(lppobj=x, dtx=p$dx, dtt=p$dt, M=p$niter,
                 weights=weights, stepnames=list(time="dt", space="dx"),
                 verbose=verbose)
  f <- a$kernel_fun
  if(fun) {
    result <- f
  } else if(!finespacing) {
    result <- as.linim(f, ...)
  } else {
    Z <- as.im(as.linim(f, ...))
    df <- a$df
    colnames(df)[colnames(df) == "seg"] <- "mapXY"
    ij <- nearest.valid.pixel(df$x, df$y, Z)
    xy <- data.frame(xc = Z$xcol[ij$col],
                     yc = Z$yrow[ij$row])
    df <- cbind(xy, df)
    result <- linim(domain(f), Z, restrict=FALSE, df=df)
  }
  attr(result, "sigma") <- sigma
  attr(result, "dx") <- a$deltax
  attr(result, "dt") <- a$deltat
  return(result)
}

FDMKERNEL <- function(lppobj, dtt, dtx, M, nsave=1,
                      weights=NULL, stepnames=list(time="dtt", space="dtx"),
                      setuponly=FALSE, verbose=FALSE) {
  ## Copyright (c) Greg McSwiggan and Adrian Baddeley 2016-2020
  ## Based on original code by Greg McSwiggan 2015-2016
  ## Internal code: parameters are now assumed to be valid.
  ## Validation code is now in 'resolve.heat.steps()'
  net2 <- as.linnet(lppobj)
  npts <- npoints(lppobj)
  if(verbose) cat("Subdividing network ...")
  lenfs <- lengths_psp(as.psp(net2))
  seg_in_lengths <- pmax(1, round(lenfs/dtx))
  new_lpp <- lixellate(lppobj, nsplit=seg_in_lengths)
  net_nodes <- as.linnet(new_lpp)
  nvert <- nvertices(net_nodes)
  if(verbose) {
    cat("Done.", fill=TRUE)
    splat("New network:")
    print(net_nodes)
    cat("Constructing update matrix A ..")
  }
  alpha <- dtt/(dtx^2)
  A <- net_nodes$m * alpha
  diag(A) <- 1 - colSums(A)
  if(verbose) {
    cat("Done.", fill=TRUE)
    splat("alpha = ", alpha)
    cat("Building initial state ..")
  }
  if(npts == 0) {
    ff <- factor(integer(0), levels=seq_len(nvert))
    ww <- numeric(0)
    U0 <- numeric(nvert)
  } else {
    tp1 <- as.numeric(new_lpp$data$tp)
    tp2 <- as.vector(rbind(1 - tp1, tp1))
    newseg <- as.integer(new_lpp$data$seg)
    vert_init_events1 <- as.vector(rbind(net_nodes$from[newseg],
                                         net_nodes$to[newseg]))
    ff <- factor(vert_init_events1, levels=seq_len(nvert))
    ww <- if(is.null(weights)) tp2 else (rep(weights, each=2) * tp2)
    ww <- ww/dtx
    U0 <- tapplysum(ww, list(ff))
  }
  if(verbose) cat("Done.", fill=TRUE)
  if(setuponly) {
    out <- list(linnet_obj   = net_nodes,
                lixelmap     = ff,   
                lixelweight  = ww,   
                Amatrix      = A,
                U0           = U0,
                deltax       = dtx,
                deltat       = dtt)
    return(out)
  }

  if(nsave == 1) {
    blockstart <- 1
    blockend <- M
  } else {
    blocksize <- ceiling(M/nsave)
    blockend <- pmin(blocksize * seq_len(nsave), M)
    blockstart <- c(1L, blockend[-nsave])
  }
  
  blocklength <- blockend - blockstart + 1L
  elapsedtime <- blockend * dtt

  if(verbose) cat("Running iterative solver ..")
  
  U <- matrix(0, nvert, nsave)

  if(npts > 0) {
    currentU <- U0
    for(i in 1:nsave) {
      v <- currentU
      nit <- blocklength[i]
      for(j in 1:nit)
        v <- A %*% v
      U[,i] <- currentU <- as.numeric(v)
    }
  }

  finalU <- U[,ncol(U)]

  if(verbose) {
    cat("Done.", fill=TRUE)
    cat("Mapping results to spatial location ..")
  }
  vert_new <- as.data.frame(vertices(net_nodes))[,c("x","y","segcoarse","tpcoarse")]
  colnames(vert_new) <- c("x", "y", "seg", "tp")
  Nodes <- lpp(vert_new, net2, check=FALSE)
  nodemap <- nnfun(Nodes)
  interpUxyst <- function(x, y, seg, tp) {
    finalU[nodemap(x,y,seg,tp)]
  }
  interpU <- linfun(interpUxyst, net2)
  df <- cbind(vert_new, data.frame(values=finalU))

  if(nsave > 1) {
    interpUxystK <- function(x, y, seg, tp, k) {
      nono <- nodemap(x,y,seg,tp)
      if(missing(k)) U[nono, ] else U[nono, k]
    }
    interpUK <- linfun(interpUxystK, net2)
  } else interpUK <- NULL

  if(verbose)
    cat("Done.", fill=TRUE)
  
  out <- list(kernel_fun  = interpU,
              elapsedtime = elapsedtime,
              tau         = sqrt(2 * elapsedtime),
              df          = df, 
              deltax      = dtx,
              deltat      = dtt,
              progressfun = interpUK)

  return(out)
}


resolve.heat.steps <-
  function(sigma, ...,
           ## main parameters (all are optional)
           ##  A=adjustable by code, F=fixed, A*=adjustable only if allow.adjust=TRUE
           dx=NULL, # spacing of sample points (A)
           dt=NULL, # time step                (A)
           niter=NULL,  # number of iterations (A*)
           iterMax=1e6, # maximum number of iterations (can be Inf) (F)
           nsave=1, # number of time points for which data should be saved (F)
                    # nsave = Inf means save all iterations, nsave = niter
           ## network information
           seglengths=NULL, # lengths of network edges
           maxdegree=NULL, # maximum vertex degree
           AMbound=NULL, # Anderson-Morley bound
           L=NULL, # optional linear network from which to extract data
           ## rules 
           finespacing=TRUE, # if FALSE, use spacing implied by pixel resolution
                             # if TRUE,  use finer spacing 
           fineNsplit=30, # finespacing rule average number of pieces per edge
           fineNlixels=100, # finespacing rule total number of pieces
           W=NULL, eps=NULL, dimyx=NULL, xy=NULL, # pixel resolution
           allow.adjust=TRUE, # 'niter' can be changed
           warn.adjust=verbose,
           verbose=TRUE,
           stepnames=list(time="dt", space="dx"))
{
  ## Based on original code by Greg McSwiggan 2015-2016

  check.1.real(sigma)  # infinite sigma is allowed
  check.1.real(nsave)  # infinite 'nsave' is allowed (will be reset to niter)
  if(is.finite(nsave)) check.1.integer(nsave)
  stopifnot(nsave >= 1)
  dx.given    <- !is.null(dx) && check.1.real(dx)
  dt.given    <- !is.null(dt) && check.1.real(dt)
  niter.given <- !is.null(niter) && check.1.integer(niter)
  nsave.given <- (nsave > 1)
  obey.nsave  <- nsave.given && is.finite(nsave) 
  save.all    <- is.infinite(nsave)
  
  one <- 1 + .Machine$double.eps # tolerance for comparisons

  if(verbose) {
    if(dx.given) splat("Given: dx =", dx)
    if(dt.given) splat("Given: dt =", dx)
    if(niter.given) splat("Given: niter =", niter)
    if(nsave.given) splat("Given: nsave =", nsave)
  }

  ## ---------- CHARACTERISTICS OF NETWORK ------------------
  if(is.null(L)) {
    check.1.integer(maxdegree)
    check.1.integer(AMbound)
  } else {
    L <- as.linnet(L)
    if(is.null(seglengths)) seglengths <- lengths_psp(as.psp(L))
    if(is.null(maxdegree) || is.null(AMbound)) {
      verdeg <- vertexdegree(L)
      maxdegree <- max(verdeg)
      AMbound <- max(verdeg[L$from] + verdeg[L$to])
    }
    if(is.null(W)) W <- Frame(L)
  }
  ## segment lengths
  nseg <- length(seglengths)
  lmin <- min(seglengths[seglengths > 0])
  lbar <- mean(seglengths[seglengths > 0])
  ltot <- sum(seglengths)
  if(verbose) {
    splat(" Network:")
    splat("    total length =", ltot)
    splat("    number of edges =", nseg)
    splat("    average nonzero edge length = ", lbar)
    splat("    shortest nonzero edge length = ", lmin)
  }

  ## ----------- NUMBER OF ITERATIONS ---------------------------------
  if(niter.given) {
    if(verbose) splat(" Validating niter =", niter)
    stopifnot(niter >= 10)
    stopifnot(niter <= iterMax)
    if(save.all) {
      nsave <- niter
    } else if(obey.nsave && ( (niter < nsave) || (niter %% nsave != 0) )) {
      if(!allow.adjust)
        stop(paste("niter =", niter, "is not a multiple of nsave =", nsave),
             call.=FALSE)
      niterOLD <- niter
      niter <- nsave * max(1L, floor(as.double(niter)/nsave))
      if(warn.adjust || verbose) {
        comment <- paste("niter was adjusted from", niterOLD, "to", niter,
                         "to ensure it is a multiple of nsave =", nsave)
        if(warn.adjust) warning(comment, call.=FALSE)
        if(verbose) splat(comment)
      }
    }
  }

  ## ----------- TIME STEP dt ------  (if given) -----------------------
  if(niter.given) {
    if(verbose) splat(" Determining dt from niter")
    dtOLD <- dt
    dt <- sigma^2/(2  * niter)
    if(dt.given) {
      if(!allow.adjust)
        stop("Only one of the arguments dt and niter should be given", call.=FALSE)
      if(warn.adjust || verbose) {
        quibble <- paste("Time step dt was adjusted from", dtOLD,
                         "to sigma^2/(2 * niter) =",
                         sigma^2, "/", 2 * niter, "=", dt)
        if(warn.adjust) warning(quibble, call.=FALSE)
        if(verbose) splat(quibble)
      }
    } else if(verbose) splat(" dt =", dt)
  } else if(dt.given) {
    if(verbose) splat(" Determining niter from dt",
                      if(obey.nsave) "and nsave" else NULL)
    stepratio <- sigma^2/(2 * dt)
    niter <- if(save.all) max(1L, round(stepratio)) else 
             nsave * max(1L, round(stepratio/nsave))
    if(niter > iterMax) {
      problem <- paste("Time step dt =", dt,
                       "implies number of iterations =", niter,
                       "exceeds maximum iterMax =", iterMax)
      if(!allow.adjust)
        stop(paste0(problem, "; increase dt or increase iterMax"),
             call.=FALSE)
      niter <- iterMax
      if(obey.nsave)
        niter <- nsave * max(1L, floor(as.double(niter)/nsave))
      if(save.all)
        nsave <- niter
      dt <- sigma^2/(2 * niter)
      if(warn.adjust || verbose) {
        comment <- paste0(problem,
                          "; niter reduced to iterMax and dt increased to ", dt)
        if(warn.adjust) warning(comment, call.=FALSE)
        if(verbose) splat(comment)
      }
    } 
    if(verbose) {
      splat(" niter =", niter)
      splat(" nsave =", nsave)
    }
  }

  ## check dt satisfies basic constraint
  if((dt.known <- dt.given || niter.given)) {
    if(verbose) splat(" Validating dt")
    dxmax <- lmin/3
    dtmax <- min(0.95 * (dxmax^2)/AMbound, sigma^2/(2 * 10), sigma * dxmax/6)
    niterMin <- max(1L, round(sigma^2/(2 * dtmax)))
    if(niterMin > iterMax) 
      stop(paste("Minimum number of iterations required is", niterMin,
                 "which exceeds iterMax =", iterMax,
                 "; increase iterMax or reduce sigma"),
           call.=FALSE)
    if(dt > dtmax) {
      #' allow rounding error
      really <- (dt > dtmax * one)
      dtOLD <- dt
      dt <- dtmax
      if(really) {
        gripe <- paste("Time step dt =", dtOLD,
                       if(allow.adjust) "reduced to" else "exceeds",
                       "maximum permitted value =", dtmax)
        if(!allow.adjust) stop(gripe, call.=FALSE)
        if(warn.adjust) warning(gripe, call.=FALSE)
        if(verbose) splat(gripe)
        if(niter.given) {
          niter <- max(1L, round(sigma^2/(2 * dt)))
          if(obey.nsave)
            niter <- nsave * max(1L, floor(as.double(niter)/nsave))
          comment <- paste("niter adjusted to", niter)
          if(warn.adjust) warning(comment, call.=FALSE)
          if(verbose) splat(comment)
          if(save.all) {
            nsave <- niter
            if(verbose) splat("  nsave = niter =", nsave)
          }
        }
      }
    }
  }
  
  #' ------------- SPACING OF SAMPLE POINTS, dx ---------------
  if(dx.given) {
    if(verbose) splat(" Validating dx =", dx)
    check.finite(dx)
    stopifnot(dx > 0)
    if(dx > lmin/3)
      stop(paste("dx must not exceed (shortest nonzero edge length)/3 =",
                 lmin/3),
           call.=FALSE)
  } else if(dt.known) {
    ## determine dx from dt
    if(verbose) splat(" Determine dx from dt")
    dx <- max(6 * dt/sigma^2, sqrt(dt * AMbound/0.95))
    if(verbose) splat(" dx =", dx)
  } else {
    #' default rule
    if(verbose) splat(" Determine dx by default rule")    
    dx <- min(lbar/fineNsplit, ltot/fineNlixels, lmin/3)
    if(verbose) {
      splat(" Mean Nonzero Edge Length/", fineNsplit, "=", lbar/fineNsplit)
      splat(" Total Network Length/", fineNlixels, "=", ltot/fineNlixels)
      splat(" Min Nonzero Edge Length/3 = ", lmin/3)
      splat(" dx = minimum of the above =", dx)
    }
    if(!finespacing && is.owin(W)) {
      W <- Frame(W)
      #' allow coarser spacing, determined by pixel size
      eps <- if(!is.null(eps)) min(eps)
             else if(!is.null(dimyx)) min(sidelengths(W)/rev(dimyx))
             else if(!is.null(xy)) with(as.mask(W, xy=xy), min(xstep, ystep))
             else min(sidelengths(W)/spatstat.options("npixel"))
      dx <- max(dx, eps/1.4)
      if(verbose) {
        splat(" Coarse spacing rule")
        splat(" Pixel size/1.4 =", eps/1.4)
      }
    } else if(verbose) splat("Fine spacing rule") 
    if(verbose) splat(" dx = ", dx)
    nlixels <- ceiling(ltot/dx)
    nlixels <- min(nlixels, .Machine$integer.max)
    dx <- ltot/nlixels
    if(verbose) {
      splat(" Rounded total number of lixels =", nlixels)
      splat(" dx =", dx)
    }
  }

  #' ------------- TIME STEP dt ----------------------------------
  dtmax <- min(0.95 * (dx^2)/AMbound, sigma^2/(2 * 10), sigma * dx/6)
  if(verbose) splat(" Applying full set of constraints")
  if(!dt.known) {
    dt <- dtmax
    if(verbose) 
      splat(" dt (determined by all constraints) = ", dt)
  } else if(dt > dtmax) {
    really <- (dt > dtmax * one)
    dtOLD <- dt
    dt <- dtmax
    if(really) {
      gripe <- paste("Time step dt =", dtOLD,
                     if(allow.adjust) "reduced to" else "exceeds",
                     "maximum permitted value =", dtmax)
      if(!allow.adjust) stop(gripe, call.=FALSE)
      if(warn.adjust) warning(gripe, call.=FALSE)
      if(verbose) splat(gripe)
      if(niter.given) {
        niter <- max(1L, round(sigma^2/(2 * dt)))
        if(obey.nsave)
          niter <- nsave * max(1L, floor(as.double(niter)/nsave))
        comment <- paste("niter adjusted to", niter)
        if(warn.adjust) warning(comment, call.=FALSE)
        if(verbose) splat(comment)
        if(save.all) {
          nsave <- niter
          if(verbose) splat(" nsave = niter =", nsave)
        }
      }
    }
  }

  #' finally determine the number of iterations, if not already done.

  if(is.null(niter)) {
    niter <- if(save.all) max(1L, round(sigma^2/(2 * dt)))
             else nsave * max(1L, round(sigma^2/(nsave * 2 * dt)))
    dt <- sigma^2/(2 * niter)
    if(verbose) {
      splat(" Number of iterations",
            paren(paste0("determined from dt",
                        if(obey.nsave) " and nsave" else NULL)),
            "=", niter)
      splat(" Updated dt =", dt)
    }
    if(save.all) {
      nsave <- niter
      if(verbose) splat(" nsave = niter =", nsave)
    }
  }

  if(niter > iterMax)
    stop(paste("Required number of iterations =", niter,
               "exceeds iterMax =", iterMax,
               "; either increase iterMax, dx, dt or reduce sigma"),
         call.=FALSE)

  alpha <- dt/dx^2
  if(verbose) splat(" alpha =", alpha)
  
  if(1 - maxdegree * alpha < 0)
    stop(paste0("Algorithm is unstable: alpha = ",
                stepnames[["time"]], "/", stepnames[["space"]], "^2 = ", alpha,
                " does not satisfy (maxdegree * alpha <= 1)",
                " where maxdegree = highest vertex degree = ", maxdegree,
                "; decrease time step ", stepnames[["time"]],
                ", or increase spacing ", stepnames[["space"]]),
         call.=FALSE)

  if(verbose) {
    splat(" Final values:")
    splat("   Time step                    dt = ", dt)
    splat("   Sample spacing               dx = ", dx)
    splat("   Number of iterations      niter = ", niter)
    splat("   Number of states saved    nsave = ", nsave)
  }
  return(list(dt=dt, dx=dx, niter=niter, nsave=nsave))
}

flatdensityfunlpp <- function(X, ..., disconnect=TRUE, weights=NULL,
                              what=c("estimate", "var", "se")) {
  stopifnot(is.lpp(X))
  trap.extra.arguments(...)
  what <- match.arg(what)
  L <- domain(X)
  nX <- npoints(X)
  if(is.null(weights)) { weights <- rep(1, nX) } else check.nvector(weights, nX)
  if(!disconnect) {
    #' constant intensity across entire network
    num <- sum(weights)
    vol <- volume(L)
    value <- switch(what,
                    estimate = num/vol,
                    var      = num/vol^2,
                    se       = sqrt(num)/vol)
    fff <- function(x, y, seg, tp) { rep(value, length(x)) }
  } else {
    #' divide L into connected components and assign each vertex to a component
    vlab <- connected(L, what="labels")
    vlab <- factor(vlab)
    #' assign each segment to a component
    slab <- vlab[L$from]
    #' total length of each component
    slen <- lengths_psp(as.psp(L))
    lenY <- tapplysum(slen, list(slab))
    #' assign points of X to components
    xlab <- slab[coords(X)$seg]
    wY <- tapplysum(weights, list(xlab))
    #' intensity of X in each component
    valY <- switch(what,
                   estimate = wY/lenY,
                   var      = wY/lenY^2,
                   se       = sqrt(wY)/lenY)
    #' function returning intensity estimate on relevant component
    fff <- function(x, y, seg, tp) { valY[ slab[seg] ] }
  }
  result <- linfun(fff, L)
  return(result)
}

flatdensityatpointslpp <- function(X, ...,
                                   leaveoneout=TRUE,
                                   disconnect=TRUE, weights=NULL,
                                   what=c("estimate", "var", "se")) {
  stopifnot(is.lpp(X))
  trap.extra.arguments(...)
  what <- match.arg(what)
  L <- domain(X)
  nX <- npoints(X)
  if(nX == 0) return(numeric(0))
  if(is.null(weights)) { weights <- rep(1, nX) } else check.nvector(weights, nX)
  if(!disconnect) {
    #' constant intensity across entire network
    totlen <- volume(L)
    numX <- rep(sum(weights), nX)
    if(leaveoneout) numX <- numX - weights
    valX <- switch(what,
                   estimate = numX/totlen,
                   var      = numX/totlen^2,
                   se       = sqrt(numX)/totlen)
  } else {
    #' divide L into connected components and assign each vertex to a component
    vlab <- connected(L, what="labels")
    vlab <- factor(vlab)
    #' assign each segment to a component
    slab <- vlab[L$from]
    #' total length of each component
    slen <- lengths_psp(as.psp(L))
    lenY <- tapplysum(slen, list(slab))
    #' assign points of X to components
    Xlab <- slab[coords(X)$seg]
    #' number of points in each component (or total weight in each component)
    sumY <- tapplysum(weights, list(Xlab))
    #' look up relevant values for each point of X
    numX <- sumY[ Xlab ]
    lenX <- lenY[ Xlab ]
    #' subtract contribution from point itself
    if(leaveoneout)
      numX <- numX - weights
    #' intensity in each component
    valX <- switch(what,
                   estimate = numX/lenX,
                   var      = numX/lenX^2,
                   se       = sqrt(numX)/lenX)
  }
  return(valX)
}
