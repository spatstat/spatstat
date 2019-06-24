#
#  lohboot.R
#
#  $Revision: 1.24 $   $Date: 2019/06/24 03:15:26 $
#
#  Loh's bootstrap CI's for local pcf, local K etc
#


spatstatLocalFunctionInfo <- function(key) {
  ## This table has to be built on the fly.
  TheTable <- list(
    pcf          = list(Global=pcf,
                        Local=localpcf,
                        L=FALSE,  inhom=FALSE,  indices=0),
    Kest         = list(Global=Kest,
                        Local=localK,
                        L=FALSE,  inhom=FALSE,  indices=0),
    Lest         = list(Global=Lest,
                        Local=localK,  # stet!
                        L=TRUE,   inhom=FALSE,  indices=0),
    pcfinhom     = list(Global=pcfinhom,
                        Local=localpcfinhom,
                        L=FALSE,  inhom=TRUE,   indices=0),
    Kinhom       = list(Global=Kinhom,
                        Local=localKinhom,
                        L=FALSE,  inhom=TRUE,   indices=0),
    Linhom       = list(Global=Linhom,
                        Local=localKinhom,  # stet!
                        L=TRUE,   inhom=TRUE,   indices=0),
    Kcross       = list(Global=Kcross,
                        Local=localKcross,
                        L=FALSE,  inhom=FALSE,  indices=2),
    Lcross       = list(Global=Lcross,
                        Local=localKcross,  # stet!
                        L=TRUE,   inhom=FALSE,  indices=2), 
    Kdot         = list(Global=Kdot,
                        Local=localKdot,
                        L=FALSE,  inhom=FALSE,  indices=1),
    Ldot         = list(Global=Ldot,
                        Local=localKdot,  # stet!
                        L=TRUE,   inhom=FALSE,  indices=1), 
    Kcross.inhom = list(Global=Kcross.inhom,
                        Local=localKcross.inhom,
                        L=FALSE,  inhom=TRUE,   indices=2),
    Lcross.inhom = list(Global=Lcross.inhom,
                        Local=localKcross.inhom,  # stet!
                        L=TRUE,   inhom=TRUE,   indices=2)
  )
  if(length(key) != 1)
    stop("Argument must be a single character string or function", call.=FALSE)
  nama <- names(TheTable)
  pos <- if(is.character(key)) {
           match(key, nama)
         } else if(is.function(key)) {
           match(list(key), lapply(TheTable, getElement, name="Global"))
         } else NULL
  if(is.na(pos)) return(NULL)
  out <- TheTable[[pos]]
  out$GlobalName <- nama[pos]
  return(out)
}

lohboot <-
  function(X,
           fun=c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom",
                 "Kcross", "Lcross", "Kdot", "Ldot",
                 "Kcross.inhom", "Lcross.inhom"),
           ...,
           block=FALSE,
           global=FALSE,
           basicboot=FALSE,
           Vcorrection=FALSE,
           confidence=0.95,
           nx = 4, ny = nx,
           nsim=200,
           type=7)
{
  stopifnot(is.ppp(X))
  ## validate 'fun'
  fun.name <- short.deparse(substitute(fun))
  if(is.character(fun)) 
    fun <- match.arg(fun)
  info <- spatstatLocalFunctionInfo(fun)
  if(is.null(info))
    stop(paste("Loh's bootstrap is not supported for the function",
               sQuote(fun.name)),
         call.=FALSE)
  fun      <- info$GlobalName
  localfun <- info$Local
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
  ## compute local functions
  f <- localfun(X, ...)
  theo <- f$theo
  ## parse edge correction info
  correction <- attr(f, "correction")
  switch(correction,
         none = {
           ckey <- clab <- "un"
           cadj <- "uncorrected"
         },
         border = {
           ckey <- "border"
           clab <- "bord"
           cadj <- "border-corrected"
         },
         translate = {
           ckey <- clab <- "trans"
           cadj <- "translation-corrected"
         },
         isotropic = {
           ckey <- clab <- "iso"
           cadj <- "Ripley isotropic corrected"
         })
  ## determine indices for Kcross etc
  types <- levels(marks(X))
  from <- resolve.1.default(list(from=types[1]), list(...))
  to <- resolve.1.default(list(to=types[2]), list(...))
  fromName <- make.parseable(paste(from))
  toName <- make.parseable(paste(to))
  ## TEMPORARY HACK for cross/dot functions.
  ## Uses a possibly temporary attribute to overwrite X with only "from" points.
  if(info$indices > 0) {
    X <- attr(f, "Xfrom")
  }
  # first n columns are the local pcfs (etc) for the n points of X
  n <- npoints(X)
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

  ## Transform to L function if required
  if(info$L) {
    theo <- sqrt(theo/pi)
    ymean <- sqrt(ymean/pi)
    hilo <- sqrt(hilo/pi)
    warn.once("lohbootLfun",
              "The calculation of confidence intervals for L functions",
              "in lohboot() has changed in spatstat 1.60-0 and later;",
              "they are now computed by transforming the confidence intervals",
              "for the corresponding K functions.")
  }

  ## create fv object
  df <- data.frame(r=f$r,
                   theo=theo,
                   ymean,
                   lo=hilo[1L,],
                   hi=hilo[2L,])
  
  colnames(df)[3L] <- ckey
  CIlevel <- paste(100 * confidence, "%% confidence", sep="")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            paste(cadj, "estimate of %s"),
            paste("lower", CIlevel, "limit for %s"),
            paste("upper", CIlevel, "limit for %s"))
  switch(fun,
         pcf={
           fname <- "g"
           yexp <- ylab <- quote(g(r))
         },
         Kest={
           fname <- "K"
           yexp <- ylab <- quote(K(r))
         },
         Lest={
           fname <- "L"
           yexp <- ylab <- quote(L(r))
         },
         pcfinhom={
           fname <- c("g", "inhom") 
           yexp <- ylab <- quote(g[inhom](r))
         },
         Kinhom={
           fname <- c("K", "inhom")
           yexp <- ylab <- quote(K[inhom](r))
         },
         Linhom={
           fname <- c("L", "inhom")
           yexp <- ylab <- quote(L[inhom](r))
         },
         Kcross={
           fname <- c("K", paste0("list(", fromName, ",", toName, ")"))
           ylab <- substitute(K[fra,til](r),
                              list(fra=fromName,til=toName))
           yexp <- substitute(K[list(fra,til)](r),
                              list(fra=fromName,til=toName))
         },
         Lcross={
           fname <- c("L", paste0("list(", fromName, ",", toName, ")"))
           ylab <- substitute(L[fra,til](r),
                              list(fra=fromName,til=toName))
           yexp <- substitute(L[list(fra,til)](r),
                              list(fra=fromName,til=toName))
         },
         Kdot={
           fname <- c("K", paste0(fromName, "~ symbol(\"\\267\")"))
           ylab <- substitute(K[fra ~ dot](r),
                              list(fra=fromName))
           yexp <- substitute(K[fra ~ symbol("\267")](r),
                              list(fra=fromName))
         },
         Ldot={
           fname <- c("L", paste0(fromName, "~ symbol(\"\\267\")"))
           ylab <- substitute(L[fra ~ dot](r),
                              list(fra=fromName))
           yexp <- substitute(L[fra ~ symbol("\267")](r),
                              list(fra=fromName))
         },
         Kcross.inhom={
           fname <- c("K", paste0("list(inhom,", fromName, ",", toName, ")"))
           ylab <- substitute(K[inhom,fra,til](r),
                              list(fra=fromName,til=toName))
           yexp <- substitute(K[list(inhom,fra,til)](r),
                              list(fra=fromName,til=toName))
         },
         Lcross.inhom={
           fname <- c("L", paste0("list(inhom,", fromName, ",", toName, ")"))
           ylab <- substitute(L[inhom,fra,til](r),
                              list(fra=fromName,til=toName))
           yexp <- substitute(L[list(inhom,fra,til)](r),
                              list(fra=fromName,til=toName))
         })
  labl <- c("r",
            makefvlabel(NULL, NULL, fname, "pois"),
            makefvlabel(NULL, "hat", fname, clab),
            makefvlabel(NULL, "hat", fname, "loCI"),
            makefvlabel(NULL, "hat", fname, "hiCI"))
  g <- fv(df, "r", ylab=ylab,
          ckey, , c(0, max(f$r)), labl, desc,
          fname=fname, yexp=yexp)
  formula(g) <- . ~ r
  fvnames(g, ".") <- c(ckey, "theo", "hi", "lo")
  fvnames(g, ".s") <- c("hi", "lo")
  unitname(g) <- unitname(X)
  g
}

