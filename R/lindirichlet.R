#'   lindirichlet.R
#'
#'   Dirichlet tessellation on a linear network
#'

lineardirichlet <- function(X) {
  stopifnot(is.lpp(X))
  #' unique points, remembering original sequence
  ii <- which(!duplicated(X))
  uX <- X[ii]
  nuX <- npoints(uX)
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
  #' upper point on interpoint distance
  huge <- sum(seglen)
  #' numerical tolerance for nnwhich
  tol <- max(sqrt(.Machine$double.eps), diameter(Frame(L))/2^20)
  #' Find data point nearest to each vertex of network
  from0 <- from - 1L
  to0   <- to - 1L
  seg0  <- seg - 1L
  z <- .C("Clinvwhichdist",
          np = as.integer(nuX),
	  sp = as.integer(seg0),
	  tp = as.double(tp),
	  nv = as.integer(nv),
	  ns = as.integer(ns),
	  from = as.integer(from0),
	  to   = as.integer(to0),
	  seglen = as.double(seglen),
	  huge = as.double(huge),
	  tol = as.double(tol),
	  dist = as.double(numeric(nv)),
	  which = as.integer(integer(nv)),
	  PACKAGE = "spatstat")
   vnndist <- z$dist
   vnnwhich <- z$which + 1L  # index into sorted unique point pattern
   vnnwhich[vnnwhich == 0] <- NA # possible if network is disconnected
   vnnlab <- coUXord$lab[vnnwhich] # index into original data pattern
   #' initialise tessellation data
   df <- data.frame(seg=integer(0),
                    t0=numeric(0),
		    t1=numeric(0),
		    tile=integer(0))
   #' split point data by segment, discarding segments which contain no points
   fseg <- factor(seg, levels=1:ns)
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
     newdf <- data.frame(seg=sygmund, t0=tcut[-m], t1=tcut[-1L], tile=labs)
     df <- rbind(df, newdf)
   }
   #' now deal with segments having no data points
   unloved <- (table(fseg) == 0)
   if(any(unloved)) {
     for(sygmund in which(unloved)) {
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
        thetile <- if(is.na(jA)) jB else jA
	newdf <- data.frame(seg=sygmund, t0=0.0, t1=1.0, tile=thetile)
      } else {
        #' split somewhere
	tx <- (dB - dA + lenf)/(2 * lenf)
	if(tx >= 0 && tx <= 1) {
  	  newdf <- data.frame(seg=sygmund,
	                      t0=c(0,tx), t1=c(tx,1), tile=c(jA, jB))
	} else if(tx < 0) {
	  # weird
	  newdf <- data.frame(seg=sygmund, t0=0.0, t1=1.0, tile=jB)
	} else {
	  # weird
	  newdf <- data.frame(seg=sygmund, t0=0.0, t1=1.0, tile=jA)
	}
      }
      df <- rbind(df, newdf)
     }
   }
   return(lintess(L, df))
}