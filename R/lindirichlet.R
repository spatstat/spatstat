#'   lindirichlet.R
#'
#'   Dirichlet tessellation on a linear network
#'
#'   $Revision: 1.9 $  $Date: 2017/11/04 03:49:18 $

lineardirichlet <- function(X) {
  stopifnot(is.lpp(X))
  #' unique points, remembering original sequence
  ii <- which(!duplicated(X))
  uX <- X[ii]
  #' local coordinates
  coUX <- coords(uX)[, c("seg", "tp")]
  #' add label from original sequence index
  coUX$lab <- ii
  #' reorder
  oo <- with(coUX, order(seg, tp))
  coUXord <- coUX[oo, , drop=FALSE]
  seg <- coUXord$seg
  tp  <- coUXord$tp
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
  #' Find data point (in sorted pattern) nearest to each vertex of network
  a <- vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol)
  vnndist  <- a$vnndist
  vnnwhich <- a$vnnwhich
  #' index back into original data pattern
  vnnlab <- coUXord$lab[vnnwhich] 
  #' compute Dirichlet tessellation
  df <- ldtEngine(nv, ns, from, to, seglen, huge,
                  coUXord,
                  vnndist, vnnwhich, vnnlab)
  return(lintess(L, df))
}

ldtEngine <- function(nv, ns, from, to, seglen, huge,  # network
                      coUXord,  # point coordinates, sorted
                      vnndist, vnnwhich, # nearest data point for each vertex
                      vnnlab) {
  #' initialise tessellation data
  seg <- integer(0)
  t0 <- numeric(0)
  t1 <- numeric(0)
  tile <- integer(0)
  #' split point data by segment, discarding segments which contain no points
  fseg <- factor(coUXord$seg, levels=1:ns)
  blist <- split(coUXord, fseg, drop=TRUE)
  #' process each segment containing data points
  for(b in blist) {
    n <- nrow(b)
    #' which segment?
    sygmund <- b$seg[[1L]]
    lenf <- seglen[sygmund]
    #' segment endpoints
    A <- from[sygmund]
    B <- to[sygmund]
    #' data points (from X) closest to endpoints
    jA <- vnnlab[A]
    jB <- vnnlab[B]
    dA <- vnndist[A]
    dB <- vnndist[B]
    #' data points (along segment) closest to endpoints
    iA <- b$lab[1L]
    iB <- b$lab[n]
    #' splits between consecutive data points
    btp <- b$tp
    tcut <- if(n < 2) numeric(0) else (btp[-1] + btp[-n])/2
    labs <- b$lab
    #' consider left endpoint
    if(jA == iA) {
      #' leftmost data point covers left endpoint
      tcut <- c(0, tcut)
    } else {
      #' cut between left endpoint and leftmost data point
      dA1 <- lenf * btp[1L]
      dx <- (dA1 - dA)/2
      if(dx > 0) {
        #' expected!
        tx <- dx/lenf
        tcut <- c(0, tx, tcut)
        labs <- c(jA, labs)
      } else {
        #' unexpected
        tcut <- c(0, tcut)
      }
    }
    #' consider right endpoint
    if(jB == iB) {
      #' rightmost data point covers right endpoint
      tcut <- c(tcut, 1)
    } else {
      #' cut between right endpoint and rightmost data point
      dB1 <- lenf * (1 - btp[n])
      dx <- (dB1 - dB)/2
      if(dx > 0) {
        #' expected!
        tx <- 1 - dx/lenf
        tcut <- c(tcut, tx, 1)
        labs <- c(labs, jB)
      } else {
        #' unexpected
        tcut <- c(tcut, 1)
      }
    }
    m <- length(tcut)
    seg <- c(seg, rep(sygmund, m-1L))
    t0 <- c(t0, tcut[-m])
    t1 <- c(t1, tcut[-1L])
    tile <- c(tile, labs)
  }
  df <- data.frame(seg=seg, t0=t0, t1=t1, tile=tile)
  #' now deal with segments having no data points
  unloved <- (table(fseg) == 0)
  if(any(unloved)) {
    unlovedt0 <- rep(0, 2*sum(unloved))
    unlovedt1 <- rep(1, 2*sum(unloved))
    unlovedseg <- unlovedtile <- rep(-1, 2*sum(unloved))
    counter <- 0
    for(sygmund in which(unloved)) {
      counter <- counter + 1
      lenf <- seglen[sygmund]
      #' segment endpoints
      A <- from[sygmund]
      B <- to[sygmund]
      #' data points (from X) closest to endpoints
      jA <- vnnlab[A]
      jB <- vnnlab[B]
      dA <- vnndist[A]
      dB <- vnndist[B]
      if(is.na(jA) || is.na(jB) || jA == jB) {
        #' entire segment is covered by one tile
        unlovedtile[counter] <- if(is.na(jA)) jB else jA
        unlovedseg[counter] <- sygmund
      } else {
        #' split somewhere
	      tx <- (dB - dA + lenf)/(2 * lenf)
	      if(tx >= 0 && tx <= 1) {
	        unlovedseg[counter] <- sygmund
	        unlovedtile[counter] <- jA
	        unlovedt1[counter] <- tx
	        counter <- counter + 1
	        unlovedseg[counter] <- sygmund
	        unlovedtile[counter] <- jB
	        unlovedt0[counter] <- tx
	      } else if(tx < 0) {
    	    # weird
	        unlovedseg[counter] <- sygmund
	        unlovedtile[counter] <- jB
	      } else {
	        # weird
	        unlovedseg[counter] <- sygmund
	        unlovedtile[counter] <- jA
	      }
      }
    }
    newdf <- data.frame(seg = unlovedseg[1:counter],
                        t0 = unlovedt0[1:counter],
                        t1 = unlovedt1[1:counter],
                        tile = unlovedtile[1:counter])
    df <- rbind(df, newdf)
  }
  return(df)
}
