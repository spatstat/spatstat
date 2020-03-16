#
# randomOnLines.R
#
# $Revision: 1.10 $  $Date: 2020/03/16 10:28:51 $
#
# Generate random points on specified lines
#

runifpointOnLines <- function(n, L, nsim=1, drop=TRUE) {
  if(!is.numeric(n) || any(n < 0) || any(n %% 1 != 0))
    stop("n should be a nonnegative integer or integers")
  if(!is.psp(L))
    L <- as.psp(L)
  W <- as.owin(L)
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    X <- datagen.runifpointOnLines(n, L)
    Y <- ppp(X$x, X$y, marks=X$marks, window=W, check=FALSE)
    result[[i]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

datagen.runifpointOnLines <- function(n, L) {
  stopifnot(is.psp(L))
  m <- length(n)
  ismarked <- (m > 1)
  if(m == 0 || (m == 1 && n == 0))
    return(data.frame(x=numeric(0),
                      y=numeric(0),
                      seg=integer(0),
                      tp=numeric(0)))
  # extract segment information
  len <- lengths_psp(L)
  sumlen <- sum(len)
  cumlen <- cumsum(len)
  cum0len <- c(0, cumlen)
  Ldf <- as.data.frame(L)
  x0 <- with(Ldf, x0)
  y0 <- with(Ldf, y0)
  dx <- with(Ldf, x1-x0)
  dy <- with(Ldf, y1-y0)
  # determine mark space
  if(ismarked) {
    markvalues <- names(n)
    if(sum(nzchar(markvalues)) < m)
      markvalues <- paste(1:m)
  }
  # initialise output data.frame
  out <- data.frame(x=numeric(0), y=numeric(0), seg=integer(0), tp=numeric(0))
  if(ismarked) 
    out <- cbind(out, data.frame(marks=character(0)))
  # generate points of each mark in turn
  for(j in 1:m) {
    if(n[[j]] > 0) {
      # generate random positions
      uu <- runif(n[[j]], min=0, max=sumlen)
      # identify segment for each point
      kk <- findInterval(uu, cum0len, rightmost.closed=TRUE, all.inside=TRUE)
      # parametric position along segment
      tt <- (uu - cum0len[kk])/len[kk]
      tt[!is.finite(tt)] <- 0
      # convert to (x,y)
      x <- x0[kk] + tt * dx[kk]
      y <- y0[kk] + tt * dy[kk]
      # assemble result
      if(!ismarked) {
        out <- data.frame(x=x, y=y, seg=kk, tp=tt)
      } else {
        outj <- data.frame(x=x, y=y, seg=kk, tp=tt, marks=markvalues[j])
        out <- rbind(out, outj)
      }
    }
  }
  if(ismarked) out$marks <- factor(out$marks, levels=markvalues)
  return(out)
}

runifpoisppOnLines <- function(lambda, L, nsim=1, drop=TRUE) {
  if(!is.numeric(lambda) || !all(is.finite(lambda) && (lambda >= 0)))
    stop("lambda should be a finite, nonnegative number or numbers")
  if(!is.psp(L))
    L <- as.psp(L)
  W <- as.owin(L)
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    X <- datagen.runifpoisppOnLines(lambda, L)
    Y <- ppp(X$x, X$y, marks=X$marks, window=W, check=FALSE)
    result[[i]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

datagen.runifpoisppOnLines <- function(lambda, L) {
  stopifnot(is.psp(L))
  mu <- lambda * sum(lengths_psp(L))
  n <- rpois(rep.int(1, length(mu)), mu)
  if(length(n) > 1)
    names(n) <- names(lambda)
  df <- datagen.runifpointOnLines(n, L)
  return(df)
}

rpoisppOnLines <- function(lambda, L, lmax=NULL, ..., nsim=1, drop=TRUE) {
  if(!is.psp(L))
    L <- as.psp(L)
  W <- as.owin(L)
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    X <- datagen.rpoisppOnLines(lambda, L, lmax=lmax, ...)
    Y <- ppp(X$x, X$y, marks=X$marks, window=W, check=FALSE)
    result[[i]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

datagen.rpoisppOnLines <- function(lambda, L, lmax=NULL, ..., check=TRUE)  {
  stopifnot(is.psp(L))
  if(is.numeric(lambda)) 
    return(datagen.runifpoisppOnLines(lambda, L))
  # ensure lambda is a list
  if(is.function(lambda) || is.im(lambda))
    lambda <- list(lambda)
  m <- length(lambda)
  # determine type of argument
  argtype <-
    if(all(unlist(lapply(lambda, is.im)))) "im" else
    if(all(unlist(lapply(lambda, is.function)))) "function" else
    stop(paste(sQuote("lambda"),
               "must be a numeric vector, a function, an image,",
               "a list of functions, or a list of images"))
  # check values of lambda
  if(argtype == "im") {
    for(j in seq_len(m)) {
      lamj <- lambda[[j]]
      if(!(lamj$type %in% c("real", "integer")))
        stop("lambda must be numeric-valued or integer-valued")
      lrange <- range(lamj)
      if(any(is.infinite(lrange)))
        stop("Infinite pixel values not permitted")
      if(lrange[1] < 0)
        stop("Negative pixel values not permitted")
    }
  }
  # determine uniform bound
  if(!is.null(lmax)) {
    stopifnot(is.numeric(lmax))
    if(length(lmax) != m) {
      if(length(lmax) == 1) {
        lmax <- rep.int(lmax, m)
      } else stop("Length of lmax does not match length of lambda")
    }
  } else {
    # compute lmax
    lmax <- numeric(m)
    for(j in seq_len(m)) {
      lamj <- lambda[[j]]
      if(is.function(lamj)) {
        X <- pointsOnLines(L, np=10000)
        lambdaX <- lamj(X$x, X$y, ...)
        lmax[j] <- max(lambdaX, na.rm=TRUE)
      } else if(is.im(lamj)) 
        lmax[j] <- max(lamj)
    }
    if(!all(is.finite(lmax)))
      stop("Infinite values of lambda obtained")
    if(any(lmax < 0))
      stop("Negative upper bound for lambda obtained")
    names(lmax) <- names(lambda)
  } 
  # Lewis-Shedler (rejection) method
  Y <- datagen.runifpoisppOnLines(lmax, L)
  n <- nrow(Y)
  if(n == 0)
    return(Y)
  # evaluate lambda at each simulated point
  if(m == 1) {
    lambda <- lambda[[1]]
    markindex <- 1
    if(is.function(lambda)) 
      lambdaY <- lambda(Y$x, Y$y, ...)
    else
      lambdaY <- safelookup(lambda, as.ppp(Y, W=as.owin(L)))
  } else {
    lambdaY <- numeric(n)
    markindex <- as.integer(Y$marks)
    for(j in seq_len(m)) {
      lamj <- lambda[[j]]
      jrows <- (markindex == j)
      Yj <- Y[jrows, , drop=FALSE]
      if(is.function(lamj)) 
        lambdaY[jrows] <- lamj(Yj$x, Yj$y, ...)
      else
        lambdaY[jrows] <- safelookup(lamj, as.ppp(Yj, W=as.owin(L)))
    }
  }
  lambdaY[is.na(lambdaY)] <- 0
  # accept/reject
  pY <- lambdaY/lmax[markindex]
  if(check) {
    if(any(pY < 0))
      warning("Negative values of lambda obtained")
    if(any(pY > 1))
      warning("lmax is not an upper bound for lambda")
  }
  retain <- (runif(n) < pY)
  Y <- Y[retain, , drop=FALSE]
  return(Y)
}

      
  
  
  
