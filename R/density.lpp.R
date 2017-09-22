#'
#'    density.lpp.R
#'
#'    Method for 'density' for lpp objects
#'
#'    Copyright (C) 2017 Greg McSwiggan and Adrian Baddeley
#'

density.lpp <- function(x, sigma, ...,
                        weights=NULL,
                        kernel="gaussian", 
                        continuous=TRUE,
                        epsilon=1e-6,
                        verbose=TRUE, debug=FALSE, savehistory=TRUE,
                        old=FALSE) {
  stopifnot(inherits(x, "lpp"))
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
  Llengths <- lengths.psp(Llines)
  # initialise stack
  stack <- data.frame(seg=integer(0), from=logical(0), 
                  distance=numeric(0), weight=numeric(0))
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
                              weight = rep(weights[i], 2L)),
                   stack)
  }
  Lfrom <- L$from
  Lto   <- L$to
  if(verbose)
    niter <- 0
  if(savehistory)
    history <- data.frame(iter=integer(0), qlen=integer(0),
                          totmass=numeric(0), maxmass=numeric(0))
  # process the stack
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
                                weight = Jweight),
                     stack)
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

density.splitppx <- function(x, sigma, ...) {
  if(!all(sapply(x, is.lpp)))
    stop("Only implemented for patterns on a linear network")
  solapply(x, density.lpp, sigma=sigma, ...)
}

PDEdensityLPP <- function(x, sigma, ..., weights=NULL, 
                          dx=NULL, dt=NULL,
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
    lenths <- lengths.psp(as.psp(L))
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
               with(as.mask(W, xy=xy), min(xstep, ystep))
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
                 iterMax=1e6, sparse=TRUE)
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
	              sparse=FALSE, dtx) {
  net2 <- as.linnet(lppobj)
#  ends1 <- net2$lines$ends
  lenfs <- lengths.psp(as.psp(net2))
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
  if(M < 10) stop("No of time iterations must be > 10, decrease dtt")
  if(M > iterMax)
    stop("No of time iterations exceeds iterMax; increase dtt or increase iterMax")

  alpha <- dtt/(dtx^2)

  A1 <- net_nodes$m *1
#  ml <- nrow(net_nodes$m)

  degree <- colSums(A1)
  dmax <- max(degree)

  A2 <- A1 * alpha
  diag(A2) <- 1 - alpha * degree
  
  if(1 - dmax*alpha < 0)
     stop("alpha must satisfy (1 - HIGHEST VERTEX DEGREE * ALPHA) > 0; decrease dtt or decrease D")

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
