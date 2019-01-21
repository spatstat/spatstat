#
#  lohboot.R
#
#  $Revision: 1.19 $   $Date: 2019/01/21 10:34:36 $
#
#  Loh's bootstrap CI's for local pcf, local K etc
#

lohboot <-
  function(X,
           fun=c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom"),
           ...,
           block=FALSE,
           global=FALSE,
           basicboot=FALSE,
           Vcorrection=FALSE,
           confidence=0.95,
           nx = 4, ny = nx,
           nsim=200,
           type=7) {
  stopifnot(is.ppp(X))
  fun.name <- short.deparse(substitute(fun))
  if(is.character(fun)) {
    fun <- match.arg(fun)
  } else if(is.function(fun)) {
    flist <- list(pcf=pcf, Kest=Kest, Lest=Lest,
                  pcfinhom=pcfinhom, Kinhom=Kinhom, Linhom=Linhom)
    id <- match(list(fun), flist)
    if(is.na(id))
      stop(paste("Loh's bootstrap is not supported for the function",
                 sQuote(fun.name)))
    fun <- names(flist)[id]
  } else stop("Unrecognised format for argument fun")
  # validate confidence level
  stopifnot(confidence > 0.5 && confidence < 1)
  alpha <- 1 - confidence
  if(!global) {
    probs <- c(alpha/2, 1-alpha/2)
    rank <- nsim * probs[2L]
  } else {
    probs <- 1-alpha
    rank <- nsim * probs
  }
  if(abs(rank - round(rank)) > 0.001)
    warning(paste("confidence level", confidence,
                  "corresponds to a non-integer rank", paren(rank),
                  "so quantiles will be interpolated"))
  n <- npoints(X)
  # compute local functions
  localfun <- switch(fun,
                     pcf=localpcf,
                     Kest=localK,
                     Lest=localL,
                     pcfinhom=localpcfinhom,
                     Kinhom=localKinhom,
                     Linhom=localLinhom)
  f <- localfun(X, ...)
  theo <- f$theo
  # parse edge correction info
  correction <- attr(f, "correction")
  switch(correction,
         none      = { ctag <- "un";    cadj <- "uncorrected" },
         border    = { ctag <- "bord";  cadj <- "border-corrected" },
         translate = { ctag <- "trans"; cadj <- "translation-corrected" },
         isotropic = { ctag <- "iso";   cadj <- "Ripley isotropic corrected" })
  # first n columns are the local pcfs (etc) for the n points of X
  y <- as.matrix(as.data.frame(f))[, 1:n]
  nr <- nrow(y)
    
  ## ---------- Modification by Christophe Biscio -----------------
  ##                (some re-coding by Adrian)
  
  if(!block) {
    ## Adrian's wrong code
    ## average local statistics
    ymean <- .rowMeans(y, na.rm=TRUE, nr, n)
    ## resample
    ystar <- matrix(, nrow=nr, ncol=nsim)
    for(i in 1:nsim) {
      ## resample n points with replacement
      ind <- sample(n, replace=TRUE)
      ## average their local statistics
      ystar[,i] <- .rowMeans(y[,ind], nr, n, na.rm=TRUE)
    }
  } else {
    ## Correct block bootstrap as described by Loh.
    W <- Window(X)
    GridTess <- quadrats(boundingbox(W), nx = nx, ny =ny)
    ## Classify points of X into grid tiles
    BlockIndex <- tileindex(X$x, X$y, GridTess)
    ## Use only 'full' blocks
    if(!is.rectangle(W)) {
      blocks <- tiles(GridTess)
      fullblocks <- sapply(blocks, is.subset.owin, B = W)
      if(sum(fullblocks)<2)
        stop("Not enough blocks are fully contained in the window", call.=FALSE)
      warning(paste("For non-rectangular windows,", 
                    "only blocks fully contained in the window are used:",
                    paste(sum(fullblocks), "were used and",
                          sum(!fullblocks),
                          "were ignored.")
                    ),
              call.=FALSE)
      ## blocks <- blocks[fullblocks]
      ## adjust classification of points of X
      indexmap <- cumsum(fullblocks)
      indexmap[!fullblocks] <- NA 
      BlockIndex <- indexmap[BlockIndex]
      ## adjust total number of points 
      n <- sum(!is.na(BlockIndex))
      BlockFactor <- factor(BlockIndex, levels=unique(indexmap[!is.na(indexmap)]))
    } else BlockFactor <- factor(BlockIndex)
    nmarks <- length(levels(BlockFactor))
    ## Average the local function values in each block
    ymarks <- by(t(y), BlockFactor, colSums, na.rm=TRUE, simplify=FALSE)
    ## Ensure empty data yield zero
    if(any(isempty <- sapply(ymarks, is.null))) 
      ymarks[isempty] <- rep(list(numeric(nr)), sum(isempty))
    ymarks <- as.matrix(do.call(cbind, ymarks)) * nmarks/n
    ## average all the marks
    ymean <- .rowMeans(ymarks, na.rm=TRUE, nr, nmarks)
    ## Average the marks in each block
    ystar <- matrix(, nrow=nr, ncol=nsim)
    for(i in 1:nsim) {
      ## resample nblocks blocks with replacement
      ind <- sample( nmarks , replace=TRUE)
      ## average their local function values
      ystar[,i] <-  .rowMeans(ymarks[,ind], nr, nmarks, na.rm=TRUE)
    }
  }
    
  ## compute quantiles
  if(!global) {
    ## pointwise quantiles
    hilo <- apply(ystar, 1, quantile,
                  probs=probs, na.rm=TRUE, type=type)
    
    ## Ripley's K function correction proposed by Loh
    if(Vcorrection && (fun=="Kest" || fun=="Kinhom")) { 
      Vcov=sqrt(1+2*pi*n*(f$r)^2/area.owin(W)) 
      hilo[1L,] <- ymean+(ymean-hilo[1L,]) / Vcov
      hilo[2L,] <- ymean+(ymean-hilo[2L,]) / Vcov
      hilo <- hilo[2:1,] # switch index so hilo[1,] is lower bound
      basicboot <- FALSE # The basic bootstrap interval is already used. Ensure that I do not modify hilo
    }
    ## So-called "basic bootstrap interval" proposed in Loh's paper;
    ## the intervals are asymptotically the same
    if(basicboot) { 
      hilo[1L,] <-  2*ymean-hilo[1L,]   
      hilo[2L,] <-  2*ymean-hilo[2L,]
      hilo <- hilo[c(2,1),] # switch index so hilo[1,] is lower bound
    } 
  } else {
    ## quantiles of deviation
    ydif <- sweep(ystar, 1, ymean)
    ydev <- apply(abs(ydif), 2, max, na.rm=TRUE)
    crit <- quantile(ydev, probs=probs, na.rm=TRUE, type=type)
    hilo <- rbind(ymean - crit, ymean + crit)
  }
    
  ## =============  End Modification by Christophe Biscio ===================
    
  # create fv object
  df <- data.frame(r=f$r,
                   theo=theo,
                   ymean,
                   lo=hilo[1L,],
                   hi=hilo[2L,])
  colnames(df)[3L] <- ctag
  CIlevel <- paste(100 * confidence, "%% confidence", sep="")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            paste(cadj, "estimate of %s"),
            paste("lower", CIlevel, "limit for %s"),
            paste("upper", CIlevel, "limit for %s"))
  clabl <- paste("hat(%s)[", ctag, "](r)", sep="")
  labl <- c("r", "%s[pois](r)", clabl, "%s[loCI](r)", "%s[hiCI](r)")
  switch(fun,
         pcf={ fname <- "g" ; ylab <- quote(g(r)) },
         Kest={ fname <- "K" ; ylab <- quote(K(r)) },
         Lest={ fname <- "L" ; ylab <- quote(L(r)) },
         pcfinhom={ fname <- "g[inhom]" ; ylab <- quote(g[inhom](r)) },
         Kinhom={ fname <- "K[inhom]" ; ylab <- quote(K[inhom](r)) },
         Linhom={ fname <- "L[inhom]" ; ylab <- quote(L[inhom](r)) })
  g <- fv(df, "r", ylab, ctag, , c(0, max(f$r)), labl, desc, fname=fname)
  formula(g) <- . ~ r
  fvnames(g, ".") <- c(ctag, "theo", "hi", "lo")
  fvnames(g, ".s") <- c("hi", "lo")
  unitname(g) <- unitname(X)
  g
}


    
  
  
  
  
