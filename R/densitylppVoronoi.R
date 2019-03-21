#'
#'    densitylppVoronoi.R
#'
#'    densityVoronoi.lpp
#'
#'    $Revision: 1.10 $  $Date: 2019/03/21 03:01:49 $
#' 

densityVoronoi.lpp <- function(X, f = 1, ..., nrep = 1, verbose = TRUE){
  # Check input
  stopifnot(is.lpp(X))
  check.1.real(f)
  if(badprobability(f))
    stop("f should be a probability between 0 and 1")
  check.1.integer(nrep)
  stopifnot(nrep >= 1)

  #' secret argument
  what <- resolve.1.default(list(what="image"), list(...))
  what <- match.arg(what, c("image", "function"))

  if(f == 0 || npoints(X) == 0) {
    #' uniform estimate
    lambdabar <- intensity(unmark(X))
    fun <- function(x, y, seg, tp) { rep(lambdabar, length(seg)) }
    g <- linfun(fun, domain(X))
    if(what == "image") g <- as.linim(g, ...)
    return(g)
  }
  if(f == 1) {
    #' Voronoi estimate
    tes <- lineardirichlet(X)
    v <- tile.lengths(tes)
    g <- as.linfun(tes, values=1/v, navalue=0)
    if(what == "image") g <- as.linim(g, ...)
    return(g)
  }

  #' Smoothed Voronoi estimate.
  #' For each repetition calculate Dirichlet tessellation;
  #' save information in a list of dataframes; and save the
  #' corresponding intensity values (i.e. inverse tile lengths)
  #' in a list of vectors.
  dflist <- tilevalueslist <- vector("list", nrep)
  blankentry <- data.frame(seg = integer(0),
                           t0 = numeric(0), t1 = numeric(0),
                           tile = integer(0))
  for (i in 1:nrep) {
    Xthin <- rthin(X, f)
    if(npoints(Xthin) == 0){
      tilevalueslist[[i]] <- 0
      dflist[[i]]         <- blankentry
    } else {
      rslt <- lineardirichlet(Xthin)
      tilevalueslist[[i]] <- 1/tile.lengths(rslt)
      dflist[[i]] <- rslt$df
    }
  }
  #' Make the result into a function on the linear network 
  fun <- function(x, y, seg, tp) {
    result <- numeric(length(seg))
    for(j in 1:nrep){
      dfj <- dflist[[j]]
      if(nrow(dfj) > 0) {
        #' classify query points by tessellation
        k <- lineartileindex(seg, tp, dfj)
        #' add Voronoi estimate
        lamj <- tilevalueslist[[j]]
        if(!anyNA(k)) {
          result <- result + lamj[k]
        } else {
          ok <- !is.na(k)
          result[ok] <- result[ok] + lamj[k[ok]]
        }
      }
    }
    return(result/(nrep*f))
  }
  g <- linfun(fun, domain(X))
  if(what == "image") g <- as.linim(g, ...)
  return(g)
}

bw.voronoi <- function(X, ..., probrange = c(0.2,0.8), nprob = 10,
                       prob = NULL, nrep = 100, verbose = TRUE){
  stopifnot(is.lpp(X))
  trap.extra.arguments(..., .Context="in bw.voronoi")
  if(!is.null(prob)) {
    stopifnot(is.numeric(prob) && is.vector(prob))
    nprob <- length(prob)
  } else {
    check.range(probrange)
    prob <- seq(from=probrange[1L], to=probrange[2L], length.out=nprob)
  }
  check.1.integer(nrep)

  nX <- npoints(X)
  cooX <- coords(X)
  segX <- cooX$seg
  tpX <- cooX$tp

  if(nX == 0) return(max(prob))
  
  if(verbose) {
    cat("Performing", nrep, "replicates... ")
    pstate <- list()
  }

  lamhat <- array(, dim=c(nX, nprob, nrep))

  for(irep in seq_len(nrep)) {
    if(verbose) pstate <- progressreport(irep, nrep, state=pstate)
    U <- runif(nX)
    for(j in seq_len(nprob)) {
      pj <- prob[j]
      retain <- (U <= pj)
      if(any(retain)) {
        Xp <- X[retain]
        #' compute leave-one-out estimates for points in Xp
        lamhat[retain, j, irep] <- looVoronoiLPP(Xp)/pj
        #' compute leave-one-out estimates for other points
        if(any(extra <- !retain)) {
          tess <- lineardirichlet(Xp)
          idx <- as.integer(lineartileindex(segX[extra], tpX[extra], tess))
          lamhat[extra, j, irep] <- 1/(pj * tile.lengths(tess)[idx])
        }
      } else lamhat[,j,irep] <- 0
    }
  }
  lamhat <- apply(lamhat, c(1,2), mean)
  cv <- colSums(log(lamhat))
  result <- bw.optim(cv, prob, iopt=which.max(cv),
                     creator="bw.voronoi",
                     criterion="Likelihood Cross-Validation",
                     unitname=NULL)
  return(result)
}


looVoronoiLPP <- function(X) {
  #' Compute leave-one-out Voronoi intensity estimate
  #' Hacked from 'lineardirichlet'
  nX <- npoints(X)
  if(nX == 0) return(numeric(0))
  #' Unique points, remembering original sequence
  ii <- which(!duplicated(X))
  uX <- X[ii]
  nuX <- npoints(uX)
  #' trivial case
  if(nuX <= 1) 
    return(rep(1/volume(domain(X)), nX))
  #' local coordinates
  coUX <- coords(uX)[, c("seg", "tp")]
  #' add label from original sequence index
  coUX$lab <- ii
  #' reorder
  oo <- with(coUX, order(seg, tp))
  coUXord <- coUX[oo, , drop=FALSE]
  seg <- coUXord$seg
  tp  <- coUXord$tp
  #' nearest neighbour of each data point, in sorted unique pattern
  nnid <- nnwhich(uX[oo])
  #' for each data point Y[i] in the sorted pattern Y,
  #' find the label of the tile that will cover Y[i] when Y[i] is removed
  neighlab <- coUXord$lab[nnid]
  #' network data
  L <- domain(X)
  nv <- nvertices(L)
  ns <- nsegments(L)
  seglen <- lengths.psp(as.psp(L))
  from <- L$from
  to   <- L$to
  #' upper bound on interpoint distance
  huge <- sum(seglen)
  #' numerical tolerance for nnwhich
  tol <- max(sqrt(.Machine$double.eps), diameter(Frame(L))/2^20)
  #' For each vertex of network, find nearest and second-nearest data points
  a <- vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax=2)
  vnndist  <- a$vnndist
  vnnwhich <- a$vnnwhich
  vnnlab <- coUXord$lab[vnnwhich] # index into original data pattern
  vnnlab <- matrix(vnnlab, ncol=2)
  #' compute result for each unique point
  lenf <- numeric(nuX)
  for(i in seq_len(nuX)) {
    #' compute Dirichlet tessellation WITHOUT point i
    coo.i <- coUXord[-i, , drop=FALSE]
    usenearest <- (vnnwhich[,1] != i)
    vnd <- ifelse(usenearest, vnndist[,1], vnndist[,2])
    vnw <- ifelse(usenearest, vnnwhich[,1], vnnwhich[,2])
    vnl <- ifelse(usenearest, vnnlab[,1], vnnlab[,2])
    adjust <- (vnw > i)
    vnw[adjust] <- vnw[adjust] - 1L
    df <- ldtEngine(nv, ns, from, to, seglen, huge,
                    coo.i, 
                    vnd, vnw, vnl)
    #' tile label of nearest neighbour
    neigh <- neighlab[i]
    #' find tile length associated with nearest neighbour of point i
    lenf[i] <- with(df, sum((tile == neigh) * seglen[seg] * (t1-t0)))
  }
  #' put back in correct place
  result <- numeric(npoints(X))
  result[ii[oo]] <- 1/lenf
  return(result)
}
