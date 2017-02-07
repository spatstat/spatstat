#
# linequad.R
#
#  $Revision: 1.13 $ $Date: 2017/02/07 08:12:05 $
#
# create quadscheme for a pattern of points lying *on* line segments

linequad <- function(X, Y, ..., eps=NULL, nd=1000, random=FALSE) {
  epsgiven <- !is.null(eps)
  if(is.lpp(X)) {
    # extract local coordinates from lpp object
    coo <- coords(X)
    mapXY <- coo$seg
    tp    <- coo$tp
    Xproj <- as.ppp(X)
    if(!missing(Y) && !is.null(Y))
      warning("Argument Y ignored when X is an lpp object")
    Y <- as.psp(X)
  } else if(is.ppp(X)) {
    # project data points onto segments
    stopifnot(is.psp(Y))
    v <- project2segment(X, Y)
    Xproj <- v$Xproj
    mapXY <- v$mapXY
    tp    <- v$tp
  } else stop("X should be an object of class lpp or ppp")
  
  # handle multitype
  ismulti <- is.multitype(X)
  if(is.marked(X) && !ismulti)
    stop("Not implemented for marked patterns")
  if(ismulti) {
    marx <- marks(X)
    flev <- factor(levels(marx))
  }
  #
  win <- as.owin(Y)
  len <- lengths.psp(Y)
  nseg <- length(len)
  if(is.null(eps)) {
    stopifnot(is.numeric(nd) && length(nd) == 1L & is.finite(nd) && nd > 0)
    eps <- sum(len)/nd
  } else
  stopifnot(is.numeric(eps) && length(eps) == 1L && is.finite(eps) && eps > 0)
  ##
  if(is.lpp(X) && spatstat.options('Clinequad')) {
    L <- as.linnet(X)
    W <- Frame(L)
    V <- vertices(L)
    nV <- npoints(V)
    coordsV <- coords(V)
    coordsX <- coords(X)
    nX <- npoints(X)
    ooX <- order(coordsX$seg)
    ndumeach <- ceiling(len/eps) + 1L
    ndummax <- sum(ndumeach)
    maxdataperseg <- max(table(factor(coordsX$seg, levels=1:nsegments(L))))
    maxscratch <- max(ndumeach) + maxdataperseg
    if(!ismulti) {
      if(!random) {
        z <- .C("Clinequad",
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ndat  = as.integer(nX),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch))
      } else {
        z <- .C("ClineRquad",
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ndat  = as.integer(nX),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch))
      }
      seqdum <- seq_len(z$ndum)
      dum <- with(z, ppp(xdum[seqdum], ydum[seqdum], window=W, check=FALSE))
      wdum <- z$wdum[seqdum]
      wdat <- numeric(nX)
      wdat[ooX] <- z$wdat
      dat <- as.ppp(X)
    } else {
      ntypes <- length(flev)
      ndummax <- ntypes * (ndummax + nX)
      maxscratch <- ntypes * maxscratch
      if(!random) {
        z <- .C("ClineMquad",
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ntypes = as.integer(ntypes),
                ndat  = as.integer(nX),
                xdat  = as.double(coordsX$x),
                ydat  = as.double(coordsX$y),
                mdat  = as.integer(as.integer(marx)-1L),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                mdum  = as.integer(integer(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch))
      } else {
        z <- .C("ClineRMquad",
                ns    = as.integer(nseg),
                from  = as.integer(L$from-1L),
                to    = as.integer(L$to-1L), 
                nv    = as.integer(nV),
                xv    = as.double(coordsV$x),
                yv    = as.double(coordsV$y), 
                eps   = as.double(eps),
                ntypes = as.integer(ntypes),
                ndat  = as.integer(nX),
                xdat  = as.double(coordsX$x),
                ydat  = as.double(coordsX$y),
                mdat  = as.integer(as.integer(marx)-1L),
                sdat  = as.integer(coordsX$seg[ooX]-1L),
                tdat  = as.double(coordsX$tp[ooX]),
                wdat  = as.double(numeric(nX)),
                ndum  = as.integer(integer(1L)),
                xdum  = as.double(numeric(ndummax)),
                ydum  = as.double(numeric(ndummax)),
                mdum  = as.integer(integer(ndummax)),
                sdum  = as.integer(integer(ndummax)),
                tdum  = as.double(numeric(ndummax)),
                wdum  = as.double(numeric(ndummax)),
                maxscratch = as.integer(maxscratch))
      }
      seqdum <- seq_len(z$ndum)
      marques <- factor(z$mdum[seqdum] + 1L, labels=flev)
      dum <- with(z, ppp(xdum[seqdum], ydum[seqdum], marks=marques,
                         window=W, check=FALSE))
      wdum <- z$wdum[seqdum]
      wdat <- numeric(nX)
      wdat[ooX] <- z$wdat
      dat <- as.ppp(X)
    }      
  } else {
    ## older, interpreted code
    ## initialise quad scheme 
    dat <- dum <- ppp(numeric(0), numeric(0), window=win)
    wdat <- wdum <- numeric(0)
    if(ismulti)
      marks(dat) <- marks(dum) <- marx[integer(0)]
    ## consider each segment in turn
    YY    <- as.data.frame(Y)
    for(i in 1:nseg) {
      ## divide segment into pieces of length eps
      ## with shorter bits at each end
      leni <- len[i]
      nwhole <- floor(leni/eps)
      if(leni/eps - nwhole < 0.5 && nwhole > 2)
        nwhole <- nwhole - 1
      rump <- (leni - nwhole * eps)/2
      brks <- c(0, rump + (0:nwhole) * eps, leni)
      nbrks <- length(brks)
      ## dummy points at middle of each piece
      sdum <- (brks[-1L] + brks[-nbrks])/2
      x <- with(YY, x0[i] + (sdum/leni) * (x1[i]-x0[i]))
      y <- with(YY, y0[i] + (sdum/leni) * (y1[i]-y0[i]))
      newdum <- list(x=x, y=y)
      ndum <- length(sdum)
      IDdum <- 1:ndum
      ## relevant data points
      relevant <- (mapXY == i)
      newdat <- Xproj[relevant]
      sdat   <- leni * tp[relevant]
      IDdat  <- findInterval(sdat, brks,
                             rightmost.closed=TRUE, all.inside=TRUE)
      ## determine weights
      w <- countingweights(id=c(IDdum, IDdat), areas=diff(brks))
      wnewdum <- w[1:ndum]
      wnewdat <- w[-(1:ndum)]
      ##
      if(!ismulti) {
        ## unmarked pattern
        dat <- superimpose(dat, newdat, W=win, check=FALSE)
        dum <- superimpose(dum, newdum, W=win, check=FALSE)
        wdat <- c(wdat, wnewdat)
        wdum <- c(wdum, wnewdum)
      } else {
        ## marked point pattern
        ## attach correct marks to data points
        marks(newdat) <- marx[relevant]
        dat <- superimpose(dat, newdat, W=win, check=FALSE)
        wdat <- c(wdat, wnewdat)
        newdum <- as.ppp(newdum, W=win, check=FALSE)
        ## replicate dummy points with each mark
        ## also add points at data locations with other marks
        for(k in seq_len(length(flev))) {
          le <- flev[k]
          avoid <- (marks(newdat) != le)
          dum <- superimpose(dum,
                             newdum %mark% le,
                             newdat[avoid] %mark% le,
                             W=win, check=FALSE)
          wdum <- c(wdum, wnewdum, wnewdat[avoid])
        }
      }
    }
  }
  ## save parameters
  dmethod <- paste("Equally spaced along each segment at spacing eps =",
                    signif(eps, 4),
                    summary(unitname(X))$plural)
  if(!epsgiven)
    dmethod <- paste0(dmethod, "\nOriginal parameter nd = ", nd)
  wmethod <- "Counting weights based on segment length"
  param <- list(dummy = list(method=dmethod),
                weight = list(method=wmethod))
  ## make quad scheme
  Qout <- quad(dat, dum, c(wdat, wdum), param=param)
  ## silently attach lines
  attr(Qout, "lines") <- Y
  return(Qout)
}

