#
#  lohboot.R
#
#  $Revision: 1.16 $   $Date: 2017/10/02 07:54:12 $
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
    
    ###### Modification by Christophe Biscio
    
  if(!block) {  # This is the loop in the former version of code
  # average them
  ymean <- .rowMeans(y, na.rm=TRUE, nr, n)
  # resample
  ystar <- matrix(, nrow=nr, ncol=nsim)
  for(i in 1:nsim) {
    # resample n points with replacement
    ind <- sample(n, replace=TRUE)
    # average their local pcfs
    ystar[,i] <- .rowMeans(y[,ind], nr, n, na.rm=TRUE)
  }
  } else{ # Block bootstrap as described by Loh.
      # Block creation for the bootstrap
      W <- Window(X)
      blocks <- tiles(quadrats(boundingbox(W), nx = nx, ny =ny))
      if(!is.rectangle(W)){
        fullblocks <- sapply(blocks, is.subset.owin, B = W)
        if(sum(fullblocks)<2){
          stop("Not enough blocks are fully contained in the window.")
        }
        warning("For non-rectangular windows only blocks fully contained in the window are used:",
                paste(sum(fullblocks), "are used, and ", sum(!fullblocks), "are omitted."))
        blocks <- blocks[fullblocks]
        n <- sum(sapply(blocks, function(w) npoints(X[w])))
      }
      
      # Average the marks in each blocks
      nmarks <- length(blocks)  # same as the number of columns in ymarks
      Xinblocks <-  lapply( 1:nmarks, FUN = function(i) {which(inside.owin(X, w=blocks[[i]]))}) # which point is in which block
      ymarks <- lapply(1:nmarks, FUN = function(i) { if(length(Xinblocks[[i]])==0) { rep(0,nr) }
        else {.rowSums(y[,Xinblocks[[i]]], nr , length(Xinblocks[[i]]), na.rm=TRUE)*nmarks/n  } } )
      ymarks <- do.call(cbind,ymarks)
      
      # average all the marks
      ymean <- .rowMeans(ymarks, na.rm=TRUE, nr, nmarks)
      
      # Average the marks in each blocks
      ystar <- matrix(, nrow=nr, ncol=nsim)
      for(i in 1:nsim) {
        # resample nblocks blocks with replacement
        ind <- sample( nmarks , replace=TRUE)
        # average their local pcfs
        ystar[,i] <-  .rowMeans(ymarks[,ind], nr, nmarks, na.rm=TRUE)
      }
    }
    
  # compute quantiles
  if(!global) {
    # pointwise quantiles
    hilo <- apply(ystar, 1, quantile,
                  probs=probs, na.rm=TRUE, type=type)
    
    # Ripley's K function correction proposed by Loh
    if(Vcorrection & (fun=="Kest" | fun=="Kinhom")) { 
      Vcov=sqrt(1+2*pi*n*(f$r)^2/area.owin(W)) 
      hilo[1L,] <-   ymean+(ymean-(hilo[1L,]) ) / Vcov
      hilo[2L,] <-   ymean+(ymean-(hilo[2L,]) ) / Vcov
      hilo <- hilo[2:1,] # switch of the index to have hilo[1,] as the lower bound and hilo[2,] as the upper bound
      basicboot <- FALSE # The basic bootstrap interval is already used. Ensure that I do not modified hilo
    }
    # So called "basic bootstrap interval" proposed in Loh's paper, the intervals are asymptotically the same
    if(basicboot)  { 
      hilo[1L,] <-  2*ymean-(hilo[1L,])   
      hilo[2L,] <-  2*ymean-(hilo[2L,])   
      hilo <- hilo[c(2,1),] # switch of the index to have hilo[1,] as the lower bound and hilo[2,] as the upper bound
    } 
    
  } else {
    # quantiles of deviation
    ydif <- sweep(ystar, 1, ymean)
    ydev <- apply(abs(ydif), 2, max, na.rm=TRUE)
    crit <- quantile(ydev, probs=probs, na.rm=TRUE, type=type)
    hilo <- rbind(ymean - crit, ymean + crit)
  }
    
    ####### End Modification by Christophe Biscio        
    
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


    
  
  
  
  
