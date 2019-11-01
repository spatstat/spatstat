#'
#'       lintessmakers.R
#'
#'   Creation of linear tessellations
#'   and intersections between lintess objects
#'
#'   $Revision: 1.2 $  $Date: 2019/11/01 09:53:54 $
#' 

divide.linnet <- local({
  
  #'  Divide a linear network into tiles demarcated by
  #' the points of a point pattern

  divide.linnet <- function(X) {
    stopifnot(is.lpp(X))
    L <- as.linnet(X)
    coo <- coords(X)
    #' add identifiers of endpoints
    coo$from <- L$from[coo$seg]
    coo$to   <- L$to[coo$seg]
    #' group data by segment, sort by increasing 'tp'
    coo <- coo[with(coo, order(seg, tp)), , drop=FALSE]
    bits <- split(coo, coo$seg)
    #' expand as a sequence of intervals
    bits <- lapply(bits, expanddata)
    #' reassemble as data frame
    df <- Reduce(rbind, bits)
    #' find all undivided segments
    other <- setdiff(seq_len(nsegments(L)), unique(coo$seg))
    #' add a single line for each undivided segment
    if(length(other) > 0)
      df <- rbind(df, data.frame(seg=other, t0=0, t1=1,
                                 from=L$from[other], to=L$to[other]))
    #' We now have a tessellation 
    #' Sort again
    df <- df[with(df, order(seg, t0)), , drop=FALSE]
    #' Now identify connected components
    #' Two intervals are connected if they share an endpoint
    #' that is a vertex of the network.
    nvert <- nvertices(L)
    nbits <- nrow(df)
    iedge <- jedge <- integer(0)
    for(iv in seq_len(nvert)) {
      joined <- with(df, which(from == iv | to == iv))
      njoin <- length(joined)
      if(njoin > 1)
        iedge <- c(iedge, joined[-njoin])
      jedge <- c(jedge, joined[-1L])
    }
    nedge <- length(iedge)
    zz <- .C("cocoGraph",
             nv = as.integer(nbits),
             ne = as.integer(nedge), 
             ie = as.integer(iedge - 1L),
             je = as.integer(jedge - 1L),
             label = as.integer(integer(nbits)), 
             status = as.integer(integer(1L)),
             PACKAGE = "spatstat")
    if (zz$status != 0) 
      stop("Internal error: connectedness algorithm did not converge")
    lab <- zz$label + 1L
    lab <- as.integer(factor(lab))
    df <- df[,c("seg", "t0", "t1")]
    df$tile <- lab
    return(lintess(L, df))
  }

  expanddata <- function(z) {
    df <- with(z,
               data.frame(seg=c(seg[1L], seg),
                          t0 = c(0, tp),
                          t1 = c(tp, 1),
                          from=NA_integer_,
                          to=NA_integer_))
    df$from[1L] <- z$from[1L]
    df$to[nrow(df)] <- z$to[1L]
    return(df)
  }

  divide.linnet
})


intersect.lintess <- function(X, Y) {
  # common refinement of two tessellations on linear network
  verifyclass(X, "lintess")
  verifyclass(Y, "lintess")
  if(!identical(as.linnet(X), as.linnet(Y)))
    stop("X and Y must be defined on the same linear network")
  L <- as.linnet(X)
  ns <- nsegments(L)
  marX <- marks(X)
  marY <- marks(Y)
  X <- X$df
  Y <- Y$df
  XY <- data.frame(seg=integer(0), t0=numeric(0), t1=numeric(0),
                   tile=character(0))
  for(seg in seq_len(ns)) {
    xx <- X[X$seg == seg, , drop=FALSE]
    yy <- Y[Y$seg == seg, , drop=FALSE]
    nxx <- nrow(xx)
    nyy <- nrow(yy)
    if(nxx > 0 && nyy > 0) {
      for(i in 1:nxx) {
        xxi <- xx[i,,drop=FALSE]
        xr <- with(xxi, c(t0, t1))
        for(j in 1:nyy) {
          yyj <- yy[j,,drop=FALSE]
          yr <- with(yyj, c(t0, t1))
          zz <- intersect.ranges(xr, yr, fatal=FALSE)
          if(!is.null(zz)) {
            XY <- rbind(XY, 
                        data.frame(seg=seg, t0=zz[1], t1=zz[2],
                                   tile=paste0(xxi$tile, ":", yyj$tile)))
          }
        }
      }
    }
  }
  out <- lintess(L, XY)
  if(!is.null(marX) || !is.null(marY)) {
    ## associate marks with TILES
    XYtiles <- levels(out$df$tile)
    posstiles <- outer(levels(X$tile), levels(Y$tile), paste, sep=":")
    m <- match(XYtiles, as.character(posstiles))
    if(anyNA(m)) stop("Internal error in matching tile labels")
    xid <- as.integer(row(posstiles))[m]
    yid <- as.integer(col(posstiles))[m]
    marXid <- marksubset(marX, xid)
    marYid <- marksubset(marY, yid)
    if(is.null(marX)) {
      marks(out) <- marYid
    } else if(is.null(marY)) {
      marks(out) <- marXid
    } else {
      if(identical(ncol(marX), 1L)) colnames(marXid) <- "marksX"
      if(identical(ncol(marY), 1L)) colnames(marYid) <- "marksY"
      marks(out) <- data.frame(marksX=marXid, marksY=marYid)
    }
  }
  return(out)
}


chop.linnet <- function(X, L) {
  X <- as.linnet(X)
  verifyclass(L, "infline")
  ## convert infinite lines to segments
  LS <- clip.infline(L, Frame(X))
  linemap <- marks(LS) # line which generated each segment
  ## extract segments of network
  XS <- as.psp(X)
  ## find crossing points (remembering provenance)
  Y <- crossing.psp(LS, XS, fatal=FALSE, details=TRUE)
  ## initialise tessellation
  Tess <- lintess(X)
  if(is.null(Y) || npoints(Y) == 0)
    return(Tess)
  ## extract info about network
  V <- vertices(X)
  startvertex <- X$from
  nXS <- nsegments(XS)
  segseq <- seq_len(nXS)
  ## allocate vertices to halfplanes defined by lines
  Vin <- whichhalfplane(L, V)
  ## group crossing-points by the infinite line that made them
  M <- marks(Y) # column names: iA, tA, jB, tB
  MM <- split(M, linemap[M$iA], drop=FALSE)
  #' for each infinite line,
  #' make the tessellation induced by this line
  for(i in seq_along(MM)) {
    Mi <- MM[[i]]
    if(is.data.frame(Mi) && (ni <- nrow(Mi)) > 0) {
      #' for each segment, determine which end is in lower left halfplane
      startsinside <- Vin[i, startvertex ]
      if(anyNA(startsinside))
        browser()
      #' find segments of X that are split, and position of split
      jj <- Mi$jB
      tt <- Mi$tB
      ss <- startsinside[jj]
      #' assemble data for these segments: 2 entries for each split segment
      inside <- paste0(i, ifelse(ss, "-", "+"))
      outside <- paste0(i, ifelse(ss, "+", "-"))
      df <- data.frame(seg=rep(jj, 2),
                       t0=c(rep(0, ni), tt),
                       t1=c(tt, rep(1, ni)),
                       tile=c(inside, outside))
      #' segments not split
      otherseg <- segseq[-jj]
      #' segments entirely inside
      otherin <- startsinside[otherseg]
      #' tack on
      if(any(otherin))
        df <- rbind(df,
                    data.frame(seg=otherseg[otherin],
                               t0=0, t1=1, tile=paste0(i, "-")))
      if(any(!otherin))
        df <- rbind(df, 
                    data.frame(seg=otherseg[!otherin],
                               t0=0, t1=1, tile=paste0(i, "+")))
      #' make it a tessellation
      Tessi <- lintess(X, df)
      #' intersect with existing
      Tess <- intersect.lintess(Tess, Tessi)
    }
  }
  return(Tess)
}

