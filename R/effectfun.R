#
#  effectfun.R
#
#   $Revision: 1.23 $ $Date: 2019/07/10 08:58:46 $
#

effectfun <- local({

  okclasses <- c("ppm", "kppm", "lppm", "dppm", "rppm", "profilepl")

effectfun <-  function(model, covname, ..., se.fit=FALSE, nvalues=256) {
  if(!inherits(model, okclasses))
    stop(paste("First argument 'model' should be a fitted model of class",
               commasep(sQuote(okclasses), " or ")),
	 call.=FALSE)
  orig.model <- model	 
  model <- as.ppm(model)
  dotargs <- list(...)
  #' determine names of covariates involved
  intern.names <-
    if(is.marked.ppm(model)) c("x", "y", "marks") else c("x", "y")
  needed.names <- variablesinformula(rhs.of.formula(formula(model)))
  #' check for clashes/quirks
  if("lambda" %in% needed.names) {
    if(is.dppm(orig.model) && (
       identical.formulae(formula(model), ~offset(log(lambda))-1) ||
       identical.formulae(formula(model), ~log(lambda)-1)
       ))
      stop("effectfun is not defined for a DPP model with fixed intensity",
           call.=FALSE)
    intensityname <- setdiff(c("Lambda", "intensity"), needed.names)[1]
  } else intensityname <- "lambda"
  ## validate the relevant covariate 
  if(missing(covname) || is.null(covname)) {
    mc <- model.covariates(model)
    if(length(mc) == 1) covname <- mc else stop("covname must be provided")
  }
  if(!(covname %in% c(intern.names, needed.names)))
    stop(paste("model does not have a covariate called", sQuote(covname)),
         call.=FALSE)
  # check that fixed values for all other covariates are provided 
  given.covs <- names(dotargs)
  if(any(uhoh <- !(needed.names %in% c(given.covs, covname)))) {
    nuh <- sum(uhoh)
    stop(paste(ngettext(nuh,
                        "A value for the covariate",
                        "Values for the covariates"),
               commasep(dQuote(needed.names[uhoh])),
               "must be provided (as",
               ngettext(nuh, "an argument", "arguments"),
               "to effectfun)"))
  }
  #' establish type and range of covariate values
  check.1.integer(nvalues)
  stopifnot(nvalues >= 128)
  N0 <- nvalues
  if(covname == "x") {
    covtype <- "real"
    W <- as.owin(data.ppm(model))
    Zr <- W$xrange
    Zvals <- seq(from=Zr[1L], to=Zr[2L], length.out=N0)
  } else if(covname == "y") {
    covtype <- "real"
    W <- as.owin(data.ppm(model))
    Zr <- W$yrange
    Zvals <- seq(from=Zr[1L], to=Zr[2L], length.out=N0)
  } else if(covname == "marks") {
    covtype <- "factor"
    Zvals <- levels(marks(data.ppm(model)))
  } else {
    # covariate is external
    if(is.data.frame(covdf <- model$covariates) && (covname %in% names(covdf))) {
      Z <- covdf[,covname]
      covtype <- typeof(Z)
      if(covtype == "double")
        covtype <- "real"
      switch(covtype,
             real={
               Zr <- range(Z)
               Zvals <- seq(from=Zr[1L], to=Zr[2L], length.out=N0)
             },
             integer={
               Zr <- range(Z)
               Zvals <- seq(from=Zr[1L], to=Zr[2L], by=ceiling((diff(Zr)+1)/N0))
             },
             factor={
               Zvals <- levels(Z)
             },
             logical={
               Zvals <- c(FALSE, TRUE)
             },
             stop(paste("Cannot handle covariate of type", dQuote(covtype)))
             )
    } else {
      Z <- getdataobjects(covname,
                          environment(formula(model)),
                          model$covariates)[[1L]]
      if(is.null(Z))
        stop(paste("Cannot find covariate", sQuote(covname)),
             call.=FALSE)
      # convert to image
      if(!is.im(Z))
        Z <- as.im(Z, W=as.owin(model))
      covtype <- Z$type
      switch(covtype,
             real={
               Zr <- summary(Z)$range
               Zvals <- seq(from=Zr[1L], to=Zr[2L], length.out=N0)
             },
             factor={
               Zvals <- levels(Z)
             },
             logical={
               Zvals <- c(FALSE, TRUE)
             },
             stop(paste("Cannot handle covariate of type", dQuote(covtype)))
             )
    }
  }
  # set up data frames of fake data for predict method
  # First set up default, constant value for each covariate
  N <- length(Zvals)
  fakeloc <- resolve.defaults(dotargs,
                              list(x=0, y=0))[c("x","y")]
  if(is.marked.ppm(model)) {
    if("marks" %in% names(dotargs)) {
      fakeloc$marks <- dotargs$marks
      dotargs <- dotargs[names(dotargs) != "marks"]
    } else {
      lev <- levels(marks(data.ppm(model)))
      fakeloc$marks <- lev[1L]
    }
  }
  fakeloc <- lapply(fakeloc, padout, N=N)
  fakecov <- lapply(dotargs, padout, N=N)
  # Overwrite value for covariate of interest
  if(covname %in% intern.names)
    fakeloc[[covname]] <- Zvals
  else fakecov[[covname]] <- Zvals
  # convert to data frame
  fakeloc <- do.call(data.frame, fakeloc)
  fakecov <- if(length(fakecov) > 0) do.call(data.frame, fakecov) else NULL
  #
  # Now predict
  pred <- predict(orig.model, locations=fakeloc, covariates=fakecov, se=se.fit)
  if(!se.fit) lambda <- pred else {
    lambda <- pred$estimate
    se     <- pred$se
    sedf <- data.frame(se =se,
                       hi = lambda + 2 * se,
                       lo = lambda - 2 * se)
  }
  #
  dfin <- if(!is.null(fakecov)) cbind(fakeloc, fakecov) else fakeloc 
  dfin <- dfin[covname]
  dflam <- data.frame(lambda=lambda)
  names(dflam) <- intensityname
  df <- cbind(dfin, dflam)
  #
  if(covtype != "real") {
    result <- df
    if(se.fit) result <- cbind(result, sedf)
  } else {
    bc <- paren(covname)
    result <- fv(df, argu=covname, 
                 ylab=substitute(lambda(X),
		                 list(X=as.name(covname),
				      lambda=as.name(intensityname))),
                 labl=c(covname,
                   paste("hat(%s)", bc)),
                 valu=intensityname, alim=Zr,
                 desc=c(paste("value of covariate", covname),
                   "fitted intensity"),
                 fname=intensityname)
    if(se.fit) {
      result <- bind.fv(result, sedf,
                        labl=c(paste("se[%s]", bc),
                          paste("%s[hi]", bc),
                          paste("%s[lo]", bc)),
                        desc=c("standard error of fitted trend",
                          "upper limit of pointwise 95%% CI for trend",
                          "lower limit of pointwise 95%% CI for trend"))
      fvnames(result, ".") <- c(intensityname, "hi", "lo")
      fvnames(result, ".s") <- c("hi", "lo")
      formula(result) <- paste(". ~ ", covname)
    }
  }
  return(result)
}

 padout <- function(x,N) { rep.int(x[1L],N) }

 effectfun
})
