#'
#'    density.lpp.R
#'
#'    Method for 'density' for lpp objects
#'
#'    Copyright (C) 2017 Greg McSwiggan and Adrian Baddeley
#'

density.lpp <- function(x, sigma=NULL, ...,
                        weights=NULL,
                        distance=c("path", "euclidean"),
                        kernel="gaussian", 
                        continuous=TRUE,
                        epsilon=1e-6,
                        verbose=TRUE, debug=FALSE, savehistory=TRUE,
                        old=FALSE) {
  stopifnot(inherits(x, "lpp"))
  distance <- match.arg(distance)

  if(distance == "euclidean") 
    return(densityQuick.lpp(x, sigma, ..., kernel=kernel, weights=weights))
  
  kernel <- match.kernel(kernel)
  if(continuous && (kernel == "gaussian") && !old)
     return(PDEdensityLPP(x, sigma, ..., weights=weights))

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
  if(savehistory)
    attr(out, "history") <- history
  return(out)
}

density.splitppx <- function(x, sigma=NULL, ...) {
  if(!all(sapply(x, is.lpp)))
    stop("Only implemented for patterns on a linear network")
  solapply(x, density.lpp, sigma=sigma, ...)
}

PDEdensityLPP <- function(x, sigma, ..., weights=NULL, 
                          dx=NULL, dt=NULL, 
                          iterMax=1e6,
                          fun=FALSE, 
                          finespacing=FALSE, finedata=finespacing) {
  stopifnot(is.lpp(x))
  L <- as.linnet(x)
  check.1.real(sigma)
  check.finite(sigma)
  if(!is.null(weights)) 
    check.nvector(weights, npoints(x))
  if(is.null(dx)) {
    #' default rule for spacing of sample points
    lenths <- lengths_psp(as.psp(L))
    lbar <- mean(lenths)
    nseg <- length(lenths)
    ltot <- lbar * nseg
    if(finespacing) {
      #' specify 30 steps per segment, on average
      dx <- lbar/30
    } else {
      #' use pixel size
      argh <- list(...)
      W <- Frame(x)
      eps <- if(!is.null(argh$eps)) {
               min(argh$eps)
             } else if(!is.null(argh$dimyx)) {
               min(sidelengths(W)/argh$dimyx)
             } else if(!is.null(argh$xy)) {
               with(as.mask(W, xy=argh$xy), min(xstep, ystep))
             } else min(sidelengths(W)/spatstat.options("npixel"))
      dx <- max(eps/1.4, lbar/30)
    }
    D <- ceiling(ltot/dx)
    D <- min(D, .Machine$integer.max)
    dx <- ltot/D
  }
  verdeg <- vertexdegree(L)
  amb <- max(verdeg[L$from] + verdeg[L$to])
  dtmax <- min(0.95 * (dx^2)/amb, sigma^2/(2 * 10), sigma * dx/6)
  if(is.null(dt)) {
    dt <- dtmax
  } else if(dt > dtmax) {
    stop(paste("dt is too large: maximum value", dtmax),
         call.=FALSE)
  }
  a <- FDMKERNEL(lppobj=x, sigma=sigma, dtx=dx, dtt=dt,
                 weights=weights,
                 iterMax=iterMax, sparse=TRUE,
                 stepnames=list(time="dt", space="dx"))
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

# Greg's code 
FDMKERNEL <- function(lppobj, sigma, dtt, weights=NULL, iterMax=5000, 
	              sparse=FALSE, dtx,
                      stepnames=list(time="dtt", space="dtx")) {
  net2 <- as.linnet(lppobj)
#  ends1 <- net2$lines$ends
  lenfs <- lengths_psp(as.psp(net2))
  seg_in_lengths <- pmax(1, round(lenfs/dtx))
  new_lpp <- lixellate(lppobj, nsplit=seg_in_lengths)
  net_nodes <- as.linnet(new_lpp)
  vvv <- as.data.frame(vertices(net_nodes)) 
  vertco_new <- vvv[, c("x", "y")]
  vertseg_new <- vvv$segcoarse # marks
  verttp_new <- vvv$tpcoarse   # marks
  if(npoints(lppobj) == 0) {
    U0 <- numeric(npoints(net_nodes$vertices))
  } else {
    tp1 <- as.numeric(new_lpp$data$tp)
    tp2 <- as.vector(rbind(1 - tp1, tp1))
    newseg <- as.integer(new_lpp$data$seg)
    vert_init_events1 <- as.vector(rbind(net_nodes$from[newseg],
                                         net_nodes$to[newseg]))
    highest_vert <- npoints(net_nodes$vertices)
    vert_numbers <- seq_len(highest_vert)
    ff <- factor(vert_init_events1, levels=vert_numbers)
    ww <- if(is.null(weights)) tp2 else (rep(weights, each=2) * tp2)
    ww <- ww/dtx
    U0 <- tapply(ww, ff, sum)
    U0[is.na(U0)] <- 0
  } 
  M <- round((sigma^2)/(2*dtt))
  if(M < 10)
    stop(paste("No of time iterations must be > 10; decrease time step",
               stepnames[["time"]]))
  if(M > iterMax)
    stop(paste0("No of time iterations = ", M,
                " exceeds maximum number iterMax = ", iterMax, 
                "; increase time step ", stepnames[["time"]],
                ", or increase iterMax"))

  alpha <- dtt/(dtx^2)

  A1 <- net_nodes$m *1
#  ml <- nrow(net_nodes$m)

  degree <- colSums(A1)
  dmax <- max(degree)

  A2 <- A1 * alpha
  diag(A2) <- 1 - alpha * degree
  
  if(1 - dmax*alpha < 0)
    stop(paste0("Algorithm is unstable: alpha = ",
                stepnames[["time"]], "/", stepnames[["space"]], "^2 = ", alpha,
                " does not satisfy (dmax * alpha <= 1)",
                " where DMAX = highest vertex degree = ", dmax,
                "; decrease time step ", stepnames[["time"]],
                ", or increase spacing ", stepnames[["space"]]))

  if(npoints(lppobj) > 0) {
    v <- as.numeric(U0)
    for(j in 1:M)
      v <- A2 %*% v
    finalU <- as.numeric(v)
  } else finalU <- U0

  vert_new <- cbind(vertco_new, vertseg_new, verttp_new)
  colnames(vert_new) <- c("x", "y", "seg", "tp")
  Nodes <- lpp(vert_new, net2, check=FALSE)
  nodemap <- nnfun(Nodes)
  interpUxyst <- function(x, y, seg, tp) {
    finalU[nodemap(x,y,seg,tp)]
  }
  interpU <- linfun(interpUxyst, net2)
  df <- cbind(vert_new, data.frame(values=finalU))
  out <- list(kernel_fun = interpU,
              df         = df, 
              deltax     = dtx,
              deltat     = dtt)
  return(out)
}


resolve.heat.steps <-
  function(sigma, ...,
           ## main parameters (all are optional)
           ## A=adjustable by code, F=fixed, F*=adjustable only if allow.adjust=TRUE
           dx=NULL, # spacing of sample points (A)
           dt=NULL, # time step (A)
           niter=NULL,  # number of iterations (F*)
           iterMax=100000, # maximum number of iterations (can be Inf) (F)
           nsave=1, # number of time points for which data should be saved (F)
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
  check.1.real(sigma)  # infinite sigma is allowed
  check.1.integer(nsave)
  stopifnot(nsave >= 1)
  dx.given    <- !is.null(dx) && check.1.real(dx)
  dt.given    <- !is.null(dt) && check.1.real(dt)
  niter.given <- !is.null(niter) && check.1.integer(niter)
  nsave.given <- (nsave > 1)
  
  one <- 1 + .Machine$double.eps # tolerance for comparisons

  if(verbose) {
    if(dx.given) splat("Given: dx =", dx)
    if(dt.given) splat("Given: dt =", dx)
    if(niter.given) splat("Given: niter =", niter)
    if(nsave.given) splat("Given: nsave =", nsave)
  }

  ## ---------- CHARACTERISTICS OF NETWORK ------------------
  if(is.null(L)) {
    check.nvector(seglengths, nsegments(L), things="edges of the network")
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
  lmin <- min(seglengths)
  lbar <- mean(seglengths)
  ltot <- lbar * nseg
  if(verbose) {
    splat(" Network:")
    splat("    total length =", ltot)
    splat("    number of edges =", nseg)
    splat("    average edge length = ", lbar)
    splat("    shortest edge length = ", lmin)
  }

  ## ----------- NUMBER OF ITERATIONS ---------------------------------
  if(niter.given) {
    if(verbose) splat(" Validating niter =", niter)
    stopifnot(niter >= 10)
    stopifnot(niter <= iterMax)
    if(nsave.given && ((niter < nsave) || (niter %% nsave != 0))) {
      if(!allow.adjust)
        stop(paste("niter =", niter, "is not a multiple of nsave =", nsave),
             call.=FALSE)
      niterOLD <- niter
      niter <- nsave * max(1L, floor(as.double(niter)/nsave))
      if(warn.adjust || verbose) {
        comment <- paste("niter was adjusted from", niterOLD, "to", niter,
                         "to ensure it is a multiple of nsave =", nsave)
        if(warn.adjust) warning(comment)
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
                      if(nsave.given) "and nsave" else NULL)
    niter <- nsave * round(sigma^2/(nsave*2*dt))
    if(niter > iterMax) {
      problem <- paste("Time step dt =", dt,
                       "implies number of iterations =", niter,
                       "exceeds maximum iterMax =", iterMax)
      if(!allow.adjust)
        stop(paste0(problem, "; increase dt or increase iterMax"),
             call.=FALSE)
      niter <- iterMax
      if(nsave.given)
        niter <- nsave * floor(as.double(niter)/nsave)
      dt <- sigma^2/(2 * niter)
      if(warn.adjust || verbose) {
        comment <- paste0(problem,
                          "; niter reduced to iterMax and dt increased to ", dt)
        if(warn.adjust) warning(comment, call.=FALSE)
        if(verbose) splat(comment)
      }
    } 
    if(verbose) splat(" niter =", niter)
  }

  ## check dt satisfies basic constraint
  if((dt.known <- dt.given || niter.given)) {
    if(verbose) splat(" Validating dt")
    dxmax <- lmin/3
    dtmax <- min(0.95 * (dxmax^2)/AMbound, sigma^2/(2 * 10), sigma * dxmax/6)
    niterMin <- round(sigma^2/(2 * dtmax))
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
          niter <- round(sigma^2/(2 * dt))
          if(nsave.given)
            niter <- nsave * floor(as.double(niter)/nsave)
          comment <- paste("niter adjusted to", niter)
          if(warn.adjust) warning(comment, call.=FALSE)
          if(verbose) splat(comment)
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
      stop(paste("dx must not exceed (shortest edge length)/3 =", lmin/3),
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
      splat(" Mean Edge Length/", fineNsplit, "=", lbar/fineNsplit)
      splat(" Total Network Length/", fineNlixels, "=", ltot/fineNlixels)
      splat(" Min Edge Length/3 = ", lmin/3)
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
    } else splat("Fine spacing rule") 
    splat(" dx = ", dx)
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
      splat(" dt (determined by constraints) = ", dt)
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
        niter <- round(sigma^2/(2 * dt))
        if(nsave.given)
          niter <- nsave * floor(as.double(niter)/nsave)
        comment <- paste("niter adjusted to", niter)
        if(warn.adjust) warning(comment, call.=FALSE)
        if(verbose) splat(comment)
      }
    }
  }

  #' finally determine the number of iterations, if not already done.

  if(is.null(niter)) {
    niter <- nsave * round(sigma^2/(nsave * 2 * dt))
    dt <- sigma^2/(2 * niter)
    if(verbose) {
      splat(" Number of iterations (determined from dt) =", niter)
      splat(" Updated dt =", dt)
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
