#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.12 $   $Date: 2020/01/05 02:46:04 $
#

rpoislpp <- function(lambda, L, ..., nsim=1, drop=TRUE) {
  if(missing(L) || is.null(L)) {
    if(inherits(lambda, c("linim", "linfun"))) {
      L <- as.linnet(lambda)
    } else if(all(sapply(lambda, inherits, what=c("linim", "linfun")))) {
      L <- unique(lapply(lambda, as.linnet))
      if(length(L) > 1)
        stop("All entries of lambda must be defined on the same network")
      L <- L[[1L]]
    } else stop("L is missing", call.=FALSE)
  } else verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  S <- as.psp(L)
  bugout <- (nsim == 1) && drop
  for(i in seq_len(nsim)) {
    X <- datagen.rpoisppOnLines(lambda, S, ...)
    Y <- lpp(X, L)
    if(bugout) return(Y)
    result[[i]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

runiflpp <- function(n, L, nsim=1, drop=TRUE) {
  verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  S <- as.psp(L)
  bugout <- (nsim == 1) && drop
  for(i in seq_len(nsim)) {
    X <- datagen.runifpointOnLines(n, S)
    Y <- lpp(X, L)
    if(bugout) return(Y)
    result[[i]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

rlpp <- function(n, f, ..., nsim=1, drop=TRUE) {
  if(inherits(f, "linfun")) 
    f <- as.linim(f, ...)
  ismulti <- FALSE
  if(length(n) > 1 && inherits(f, "linim")) {
    f <- rep(list(f), length(n))
    ismulti <- TRUE
  } else if(!inherits(f, "linim") && is.list(f) &&
            all(sapply(f, inherits, what=c("linim", "linfun")))) {
    #' f is a list of densities for each type of point
    if(length(n) == 1) {
      n <- rep(n, length(f))
    } else stopifnot(length(n) == length(f))
    ismulti <- TRUE
  }
  if(ismulti) {
    Y <- mapply(rlpp, n=as.list(n), f=f,
                MoreArgs=list(nsim=nsim, drop=FALSE, ...),
                SIMPLIFY=FALSE)
    names(Y) <- names(f) %orifnull% as.character(seq(along=f))
    Z <- do.call(mapply, c(list(superimpose), Y, list(SIMPLIFY=FALSE)))
    result <- simulationresult(Z, nsim, drop)
    return(result)
  }
  if(!inherits(f, "linim"))
    stop("f should be a linfun or linim object")
  if(length(n) > 1) {
    flist <- rep(list(f), length(n))
    return(rlpp(n, flist, nsim=nsim, drop=drop, ...))
  }
  check.1.integer(nsim)
  if(nsim <= 0) return(list())
  #' extract data
  L <- as.linnet(f)
  df <- attr(f, "df")
  seglen <- lengths.psp(as.psp(L))
  #' sort into segments, left-to-right within segments
  df <- df[order(df$mapXY, df$tp), , drop=FALSE]
  nr <- nrow(df)
  fvals <- df$values
  if(anyNA(fvals)) stop("f has some NA values")
  if(min(fvals) < 0) stop("f has some negative values")
  #' find interval corresponding to each sample point
  sameseg <- (diff(df$mapXY) == 0)
  sharenext     <- c(sameseg, FALSE)
  shareprevious <- c(FALSE, sameseg)
  tcur   <- df$tp
  tnext  <- c(tcur[-1], NA)
  tprev  <- c(NA, tcur[-nr])
  tleft  <- ifelse(shareprevious, (tcur + tprev)/2, 0)
  tright <- ifelse(sharenext,     (tcur + tnext)/2, 1)
  #' compute probability of each interval
  probs <- fvals * (tright - tleft) * seglen[df$mapXY]
  probs <- probs/sum(probs)
  #' 
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    #' sample intervals and place point uniformly in each interval
    ii <- sample.int(nr, size=n, replace=TRUE, prob=probs)
    seg <- df[ii, "mapXY"]
    tp  <- runif(n, tleft[ii], tright[ii])
    result[[isim]] <- as.lpp(seg=seg, tp=tp, L=L)
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}
