#
# kernel smoothing on linear network
# (Okabe algorithms)
#  Computationally expensive unless sigma is very small
#  You Have Been Warned
#
#   $Revision: 1.6 $  $Date: 2016/02/17 00:30:43 $
#

density.lpp <- function(x, sigma, ...,
                        weights=NULL,
                        kernel="gaussian", 
                        continuous=TRUE,
                        epsilon=1e-6,
                        verbose=TRUE, debug=FALSE, savehistory=TRUE) {
  stopifnot(inherits(x, "lpp"))
  kernel <- match.kernel(kernel)
  L <- as.linnet(x)
  # weights
  np <- npoints(x)
  if(is.null(weights)) {
    weights <- rep(1, np)
  } else {
    stopifnot(is.numeric(weights))
    check.nvector(weights, np, oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, np) 
  }
  # pixellate linear network
  Llines <- as.psp(L)
  linemask <- as.mask.psp(Llines, ...)
  lineimage <- as.im(linemask)
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
    # Gaussian density on segment containing x[i]
    relevant <- (projmap$mapXY == segi)
    values[relevant] <- values[relevant] +
      dkernel(len * (projmap$tp[relevant] - tpi),
              kernel=kernel, sd=sigma)
    # push the two tails onto the stack
    stack <- rbind(data.frame(seg = c(segi, segi),
                              from  = c(TRUE, FALSE), 
                              distance = len * c(tpi, 1-tpi),
                              weight = rep(weights[i], 2)),
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
                       data.frame(iter=nrow(history)+1,
                                  qlen=nrow(stack),
                                  totmass=totmass,
                                  maxmass=maxmass))
    if(verbose) {
      niter <- niter + 1
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
    H  <- stack[1, , drop=FALSE]
    stack <- stack[-1, , drop=FALSE]
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
