#
#  slrm.R
#
#  Spatial Logistic Regression
#
#  $Revision: 1.27 $   $Date: 2016/02/11 10:17:12 $
#

slrm <- function(formula, ..., data=NULL, offset=TRUE, link="logit",
                 dataAtPoints=NULL, splitby=NULL) {
  
  # remember call
  CallInfo <- list(callstring = short.deparse(sys.call()),
                   cl = match.call(),
                   formula = formula,
                   offset=offset,
                   link=link,
                   splitby=splitby,
                   dotargs=list(...))
  if(!(link %in% c("logit", "cloglog")))
    stop(paste("Unrecognised link", dQuote(link)))

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(formula, "formula"))
    stop(paste("Argument", dQuote("formula"), "should be a formula"))

  # check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Argument", sQuote("formula"),
               "must have a left hand side"))
  Yname <- formula[[2]]
  trend <- rhs <- formula[c(1,3)]
  if(!is.name(Yname))
    stop("Left hand side of formula should be a single name")
  Yname <- paste(Yname)
  if(!inherits(trend, "formula"))
    stop("Internal error: failed to extract RHS of formula")

  varnames <- unique(variablesinformula(trend))
  specials <- c("x", "y", "logpixelarea")
  covnames <- varnames[!(varnames %in% specials)]

  # add 'splitby' to covariate names
  if(!is.null(splitby)) {
    if(!is.character(splitby) || length(splitby) != 1)
      stop("splitby should be a single character string")
    covnames <- unique(c(covnames, splitby))
  }

  CallInfo$responsename <- Yname
  CallInfo$varnames     <- varnames
  CallInfo$covnames     <- covnames
  
  # Parent environment
  parenv <- environment(formula)

  ########  FIND DATA AND RESHAPE #######################

  Data <- slr.prepare(CallInfo, parenv, data, dataAtPoints, splitby)

#  W  <- Data$W
  df <- Data$df
  
  ########  FIT MODEL ###############################

  dformula <- formula
  if(offset) {
    # insert offset term in formula
    rhs <- paste(as.character(rhs), collapse=" ")
    rhs <- paste(c(rhs, "offset(logpixelarea)"), collapse="+")
    dformula <- as.formula(paste(Yname, rhs))
  }

  linkname <- link
  FIT  <- glm(dformula, family=binomial(link=linkname),
              data=df, na.action=na.exclude)

  result <- list(call     = CallInfo$cl,
                 CallInfo = CallInfo,
                 Data     = Data,
                 Fit      = list(FIT=FIT, dformula=dformula),
                 terms    = terms(formula))

  class(result) <- c("slrm", class(result))
  return(result)
}

################ UTILITY TO FIND AND RESHAPE DATA #################

slr.prepare <- function(CallInfo, envir, data,
                        dataAtPoints=NULL, splitby=NULL,
                        clip=TRUE) {
  # CallInfo is produced by slrm()
  # envir is parent environment of model formula
  # data  is 'data' argument that takes precedence over 'envir'
  # 'clip' is TRUE if the data should be clipped to the domain of Y
  Yname    <- CallInfo$responsename
#  varnames <- CallInfo$varnames
  covnames <- CallInfo$covnames
  dotargs  <- CallInfo$dotargs
  #
  getobj <- function(nama, env, dat) {
    if(!is.null(dat) && !is.null(x <- dat[[nama]]))
      return(x)
    else return(get(nama, envir=env))
  }
  # Get the response point pattern Y 
  Y <- getobj(Yname, envir, data)
  if(!is.ppp(Y))
    stop(paste("The response", sQuote(Yname), "must be a point pattern"))
  #
  if(!is.null(dataAtPoints)) {
    dataAtPoints <- as.data.frame(dataAtPoints)
    if(nrow(dataAtPoints) != npoints(Y))
      stop(paste("dataAtPoints should have one row for each point in",
                 dQuote(Yname)))
  }
  # Find the covariates
  ncov <- length(covnames)
  covlist <- lapply(as.list(covnames), getobj, env = envir, dat=data)
  names(covlist) <- covnames
  # Each covariate should be an image, a window, a function, or a single number
  if(ncov == 0) {
    isim <- isowin <- ismask <- isfun <- isnum <- isspatial <- israster <- logical(0)
  } else {
    isim  <- sapply(covlist, is.im)
    isowin  <- sapply(covlist, is.owin)
    ismask  <- sapply(covlist, is.mask)
    isfun  <- sapply(covlist, is.function)
    isspatial <- isim | isowin | isfun
    israster <- isim | ismask
    isnum <- sapply(covlist, is.numeric) & (lengths(covlist) == 1)
  }
  if(!all(ok <- (isspatial | isnum))) {
    n <- sum(!ok)
    stop(paste(ngettext(n, "The argument", "Each of the arguments"),
               commasep(sQuote(covnames[!ok])),
               "should be either an image, a window, or a single number"))
  }
  # 'splitby' 
  if(!is.null(splitby)) {
    splitwin <- covlist[[splitby]]
    if(!is.owin(splitwin))
      stop("The splitting covariate must be a window")
    # ensure it is a polygonal window
    covlist[[splitby]] <- splitwin <- as.polygonal(splitwin)
    # delete splitting covariate from lists to be processed
    issplit <- (covnames == splitby)
    isspatial[issplit] <- FALSE
    israster[issplit] <- FALSE
  }
  # 
#  nnum <- sum(isnum)
#  nspatial <- sum(isspatial)
  nraster <- sum(israster)
  #
  numlist <- covlist[isnum]
  spatiallist <- covlist[isspatial]
  rasterlist <- covlist[israster]
  #
  numnames <- names(numlist)
  spatialnames <- names(spatiallist)
#  rasternames <- names(rasterlist)
  #
  
  ########  CONVERT TO RASTER DATA  ###############################

  convert <- function(x,W) {
    if(is.im(x) || is.function(x)) return(as.im(x,W))
    if(is.owin(x)) return(as.im(x, W, value=TRUE, na.replace=FALSE))
    return(NULL)
  }

  # determine spatial domain & common resolution: convert all data to it
  if(length(dotargs) > 0 || nraster == 0) {
    # Pixel resolution is determined by explicit arguments
    if(clip) {
      # Window extent is determined by response point pattern
      D <- as.owin(Y)
    } else {
      # Window extent is union of domains of data
      domains <- lapply(append(spatiallist, list(Y)), as.owin)
      D <- do.call(union.owin, domains)
    }
    # Create template mask
    W <- do.call(as.mask, append(list(D), dotargs))
    # Convert all spatial objects to this resolution
    spatiallist <- lapply(spatiallist, convert, W=W)
  } else {
    # Pixel resolution is determined implicitly by covariate data
    W <- do.call(commonGrid, rasterlist)
    if(clip) {
      # Restrict data to spatial extent of response point pattern
      W <- intersect.owin(W, as.owin(Y))
    }
    # Adjust spatial objects to this resolution
    spatiallist <- lapply(spatiallist, convert, W=W)
  }
  # images containing coordinate values
  xcoordim <- as.im(function(x,y){x}, W=W)
  ycoordim <- as.im(function(x,y){y}, W=W)
  #
  # create a list of covariate images, with names as in formula
  covimages <- append(list(x=xcoordim, y=ycoordim), spatiallist)

  basepixelarea <- W$xstep * W$ystep

  ########  ASSEMBLE DATA FRAME  ###############################

  if(is.null(splitby)) {
    df <- slrAssemblePixelData(Y, Yname, W,
                               covimages, dataAtPoints, basepixelarea)
    sumYloga <- Y$n * log(basepixelarea)
    serial <- attr(df, "serial")
  } else {
    # fractional pixel areas
    pixsplit <- pixellate(splitwin, W)
    splitpixelarea <- as.vector(as.matrix(pixsplit))
    # determine which points of Y are inside/outside window
    ins <- inside.owin(Y$x, Y$y, splitwin)
    # split processing
    dfIN <- slrAssemblePixelData(Y[ins], Yname, W, covimages,
                                 dataAtPoints[ins, ], splitpixelarea)
    serialIN <- attr(dfIN, "serial")
    dfIN[[splitby]] <- TRUE
    dfOUT <- slrAssemblePixelData(Y[!ins], Yname, W, covimages,
                                  dataAtPoints[!ins, ],
                                  basepixelarea - splitpixelarea)
    serialOUT <- attr(dfOUT, "serial")
    dfOUT[[splitby]] <- FALSE
    df <- rbind(dfIN, dfOUT)
    serial <- c(serialIN, serialOUT)
    # sum of log pixel areas associated with points
    Ysplit <- pixsplit[Y]
    sumYloga <- sum(log(ifelseXY(ins, Ysplit, basepixelarea - Ysplit)))
  }
  
  # tack on any numeric values
  df <- do.call(cbind, append(list(df), numlist))
  
  ### RETURN ALL 
  Data <- list(response=Y,
               covariates=covlist,
               spatialnames=spatialnames,
               numnames=numnames,
               W=W,
               df=df,
               serial=serial,
               sumYloga=sumYloga,
               dataAtPoints=dataAtPoints)
  return(Data)
}

#  
slrAssemblePixelData <- function(Y, Yname, W,
                                 covimages, dataAtPoints, pixelarea) {
  # pixellate point pattern
  Z <- pixellate(Y, W=W)
  Z <- eval.im(as.integer(Z>0))
  # overwrite pixel entries for data points using exact values
  # coordinates
  xcoordim <- covimages[["x"]]
  ycoordim <- covimages[["y"]]
  xcoordim[Y] <- Y$x
  ycoordim[Y] <- Y$y
  covimages[["x"]] <- xcoordim
  covimages[["y"]] <- ycoordim
  # overwrite pixel entries
  if(!is.null(dataAtPoints)) {
    enames <- colnames(dataAtPoints)
    relevant <- enames %in% names(covimages)
    for(v in enames[relevant]) {
      cova <- covimages[[v]]
      cova[Y] <- dataAtPoints[, v, drop=TRUE]
      covimages[[v]] <- cova
    }
  }
  # assemble list of all images
  Ylist <- list(Z)
  names(Ylist) <- Yname
  allimages <- append(Ylist, covimages)
  # extract pixel values of each image
  pixelvalues <-
    function(z) {
      v <- as.vector(as.matrix(z))
      if(z$type != "factor") return(v)
      lev <- levels(z)
      return(factor(v, levels=seq_along(lev), labels=lev))
    }
  pixdata <- lapply(allimages, pixelvalues)
  df <- as.data.frame(pixdata)
  serial <- seq_len(nrow(df))
  # add log(pixel area) column
  if(length(pixelarea) == 1) {
    df <- cbind(df, logpixelarea=log(pixelarea))
  } else {
    ok <- (pixelarea > 0)
    df <- cbind(df[ok, ], logpixelarea=log(pixelarea[ok]))
    serial <- serial[ok]
  }
  attr(df, "serial") <- serial
  return(df)
}

is.slrm <- function(x) {
  inherits(x, "slrm")
}

coef.slrm <- function(object, ...) {
  coef(object$Fit$FIT)
}

print.slrm <- function(x, ...) {
  lk <- x$CallInfo$link
  switch(lk,
         logit= {
           splat("Fitted spatial logistic regression model")
         },
         cloglog= {
           splat("Fitted spatial regression model (complementary log-log)")
         },
         {
           splat("Fitted spatial regression model")
           splat("Link =", dQuote(lk))
         })
  cat("Formula:\t")
  print(x$CallInfo$formula)
  splat("Fitted coefficients:")
  print(coef(x))
  return(invisible(NULL))
}

logLik.slrm <- function(object, ..., adjust=TRUE) {
  FIT  <- object$Fit$FIT
  ll <- -deviance(FIT)/2
  if(adjust) {
    sumYloga <- object$Data$sumYloga
    ll <- ll - sumYloga
  }
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  return(ll)
}

fitted.slrm <- function(object, ...) {
  if(length(list(...)) > 0)
    warning("second argument (and any subsequent arguments) ignored")
  predict(object, type="probabilities")
}

predict.slrm <- function(object, ..., type="intensity",
                         newdata=NULL, window=NULL) {
  type <- pickoption("type", type,
                     c(probabilities="probabilities",
                       link="link",
                       intensity="intensity",
                       lambda="intensity"))
  
  FIT  <- object$Fit$FIT
  link <- object$CallInfo$link
  W    <- object$Data$W
  df   <- object$Data$df
  loga <- df$logpixelarea

  if(is.null(newdata) && is.null(window)) {
    # fitted values from existing fit
    switch(type,
           probabilities={
             values <- fitted(FIT)
           },
           link={
             values <- predict(FIT, type="link")
           },
           intensity={
             # this calculation applies whether an offset was included or not
             if(link == "cloglog") {
               linkvalues <- predict(FIT, type="link")
               values <- exp(linkvalues - loga)
             } else {
               probs <- fitted(FIT)
               values <- -log(1-probs)/exp(loga)
             }
           }
           )
    out <- im(values, xcol=W$xcol, yrow=W$yrow, unitname=unitname(W))
    return(out)
  } else {
    # prediction using new values
    # update arguments that may affect pixel resolution
    CallInfo <- object$CallInfo
    CallInfo$dotargs <- resolve.defaults(list(...), CallInfo$dotargs)
    #
    if(!is.null(window)) {
      # insert fake response in new window
      if(is.null(newdata)) newdata <- list()
      window <- as.owin(window)
      newdata[[CallInfo$responsename]] <- ppp(numeric(0), numeric(0),
                                            window=window)
    }
    # process new data
    newData <- slr.prepare(CallInfo, environment(CallInfo$formula), newdata,
                           clip=!is.null(window))
    newdf   <- newData$df
    newW    <- newData$W
    newloga <- newdf$logpixelarea
    # avoid NA etc
    npixel <- nrow(newdf)
    ok <- complete.cases(newdf)
    if(!all(ok)) {
      newdf   <- newdf[ok, , drop=FALSE]
      newloga <- newloga[ok]
    }
    # compute link values
    linkvalues <- predict(FIT, newdata=newdf, type="link")
    # transform to desired scale
    linkinv <- family(FIT)$linkinv
    switch(type,
           probabilities={
             values <- linkinv(linkvalues)
           },
           link={
             values <- linkvalues
           },
           intensity={
             # this calculation applies whether an offset was included or not
             if(link == "cloglog") {
               values <- exp(linkvalues - newloga)
             } else {
               probs <- linkinv(linkvalues)
               values <- -log(1-probs)/exp(newloga)
             }
           }
           )
    # form image
    v <- rep.int(NA_real_, npixel)
    v[ok] <- values
    out <- im(v, xcol=newW$xcol, yrow=newW$yrow, unitname=unitname(W))
    return(out)
  }
}

plot.slrm <- function(x, ..., type="intensity") {
  xname <- short.deparse(substitute(x))
  y <- predict(x, type=type)
  do.call(plot.im, resolve.defaults(list(x=y), list(...), list(main=xname)))
}

formula.slrm <- function(x, ...) {
  f <- x$CallInfo$formula
  return(f)
}

terms.slrm <- function(x, ...) {
  terms(formula(x), ...)
}

labels.slrm <- function(object, ...) {
  # extract fitted trend coefficients
  co <- coef(object)
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  if(length(lab) == 0)
    return(character(0))
  # model matrix
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # 'ass' associates coefficients with model terms
  # except ass == 0 for the Intercept
  coef.ok <- is.finite(co)
  relevant <- (ass > 0) 
  okterms <- unique(ass[coef.ok & relevant])
  return(lab[okterms])
}

extractAIC.slrm <- function (fit, scale = 0, k = 2, ...)
{
    edf <- length(coef(fit))
    aic <- AIC(fit)
    c(edf, aic + (k - 2) * edf)
}

model.matrix.slrm <- function(object,..., keepNA=TRUE) {
  FIT <- object$Fit$FIT
  mm <- model.matrix(FIT, ...)
  if(!keepNA)
    return(mm)
  df <- object$Data$df
  comp <- complete.cases(df)
  if(all(comp))
    return(mm)
  if(sum(comp) != nrow(mm))
      stop("Internal error in patching NA's")
  mmplus <- matrix(NA, nrow(df), ncol(mm))
  mmplus[comp, ] <- mm
  return(mmplus)
}

model.images.slrm <- function(object, ...) {
  mm <- model.matrix(object, ...)
  mm <- as.data.frame(mm)
  Data <- object$Data
  W      <- Data$W
  serial <- Data$serial
  splitby <- object$CallInfo$splitby
  blank   <- as.im(NA_real_, W)
  assignbyserial <- function(values, serial, template) {
    Z <- template
    Z$v[serial] <- values
    return(Z)
  }
  if(is.null(splitby)) {
    result <- lapply(as.list(mm), assignbyserial, serial=serial, template=blank)
  } else {
    df <- Data$df
    IN <- as.logical(df[[splitby]])
    OUT <- !IN
    mmIN <- mm[IN, , drop=FALSE]
    mmOUT <- mm[OUT, , drop=FALSE]
    resultIN <- lapply(as.list(mmIN), assignbyserial,
                       serial=serial[IN], template=blank)
    resultOUT <- lapply(as.list(mmOUT), assignbyserial,
                       serial=serial[OUT], template=blank)
    names(resultIN) <- paste(names(resultIN), splitby, "TRUE", sep="")
    names(resultOUT) <- paste(names(resultOUT), splitby, "FALSE", sep="")
    result <- c(resultIN, resultOUT)
  }
  return(as.solist(result))
}

update.slrm <- function(object, ..., evaluate=TRUE, env=parent.frame()) {
  e <- update.default(object, ..., evaluate=FALSE)
  if(evaluate)
    e <- eval(e, envir=env)
  return(e)
}

anova.slrm <- local({

  anova.slrm <- function(object, ..., test=NULL) {
    objex <- append(list(object), list(...))
    if(!all(unlist(lapply(objex, is.slrm))))
      stop("Some arguments are not of class slrm")
    fitz <- lapply(objex, getFIT)
    do.call(anova, append(fitz, list(test=test)))
  }

  getFIT <- function(z) {z$Fit$FIT}

  anova.slrm
})


vcov.slrm <- function(object, ..., what=c("vcov", "corr", "fisher", "Fisher")) {
  stopifnot(is.slrm(object))
  what <- match.arg(what)
  vc <- vcov(object$Fit$FIT)
  result <- switch(what,
                   vcov = vc,
                   corr = {
                     sd <- sqrt(diag(vc))
                     vc / outer(sd, sd, "*")
                   },
                   fisher=,
                   Fisher={
                     solve(vc)
                   })
  return(result)
}

unitname.slrm <- function(x) {
  return(unitname(x$Data$response))
}

"unitname<-.slrm" <- function(x, value) {
  unitname(x$Data$response) <- value
  return(x)
}

is.stationary.slrm <- function(x) {
  fo <- formula(x)
  trend <- fo[c(1,3)]
  return(identical.formulae(trend, ~1))
}

is.poisson.slrm <- function(x) { TRUE }


simulate.slrm <- function(object, nsim=1, seed=NULL, ...,
                          window=NULL, covariates=NULL, 
                          verbose=TRUE, drop=FALSE) {
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  starttime <- proc.time()
  
  # determine simulation window and compute intensity
  if(!is.null(window))
    stopifnot(is.owin(window))
  lambda <- predict(object, type="intensity", newdata=covariates, window=window)

  # max lambda (for efficiency)
  summ <- summary(lambda)
  lmax <- summ$max + 0.05 * diff(summ$range)

  # run
  out <- list()
  verbose <- verbose && (nsim > 1)
  if(verbose) {
    cat(paste("Generating", nsim, "simulations... "))
    pstate <- list()
  }
  for(i in 1:nsim) {
    out[[i]] <- rpoispp(lambda, lmax=lmax)
    if(verbose) pstate <- progressreport(i, nsim, state=pstate)
  }
  # pack up
  if(nsim == 1 && drop) {
    out <- out[[1]]
  } else {
    out <- as.solist(out)
    if(nsim > 0)
      names(out) <- paste("Simulation", 1:nsim)
  }
  out <- timed(out, starttime=starttime)
  attr(out, "seed") <- RNGstate
  return(out)
}
