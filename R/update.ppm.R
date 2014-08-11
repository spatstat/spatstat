#
#  update.ppm.R
#
#
#  $Revision: 1.40 $    $Date: 2014/04/11 09:47:09 $
#
#
#

update.ppm <- local({

  newformula <- function(old, change, eold, enew) {
    old <- if(is.null(old)) ~1 else eval(old, eold)
    change <- if(is.null(change)) ~1 else eval(change, enew)
    old <- as.formula(old, env=eold)
    change <- as.formula(change, env=enew)
    update.formula(old, change)
  }
  
  update.ppm <- function(object, ...,
                         fixdummy=TRUE, use.internal=NULL,
                         envir=environment(terms(object))) {
    verifyclass(object, "ppm")
    aargh <- list(...)

    call <- getCall(object)
    if(!is.call(call))
      stop(paste("Internal error - getCall(object) is not of class",
                 sQuote("call")))
    callframe <- object$callframe
    callstring <- short.deparse(sys.call())
  
    ## Special cases 
    ## (1) no new information
    if(length(aargh) == 0) {
      result <- eval(call, as.list(envir), enclos=callframe)
      result$callframe <- callframe
      return(result)
    }

    ## (2) model can be updated using existing covariate data frame
    ismpl <- with(object, method == "mpl" && !is.null(fitter) && fitter %in% c("gam", "glm"))
    only.fmla <- length(aargh) == 1 && inherits(fmla <- aargh[[1]], "formula")
    if(ismpl && only.fmla) {
      ## This is a dangerous hack! 
      glmdata <- object$internal$glmdata
      ## check whether data for new variables are available
      ## (this doesn't work with things like 'pi')
      vars.available <- c(colnames(glmdata), names(object$covfunargs))
      if(all(variablesinformula(fmla) %in% c(".", vars.available))) {
        ## we can update using internal data
        FIT <- object$internal$glmfit
        orig.env <- environment(FIT$terms)
        ## ensure fmla has only a right hand side
        if(!is.null(lh <- lhs.of.formula(fmla)))
          warning("Ignored left side of formula")
        fmla <- rhs.of.formula(fmla)
        ## update formulae using "." rules
        trend <- newformula(object$trend, fmla, callframe, envir)
        fmla  <- newformula(formula(FIT), fmla, callframe, envir)
        ## update GLM/GAM fit 
        upd.glm.call <- update(FIT, fmla, evaluate=FALSE)
        FIT <- eval(upd.glm.call, envir=orig.env)
        environment(FIT$terms) <- orig.env
        object$internal$glmfit <- FIT
        ## update entries of object
        object$trend <- trend
        object$terms <- terms(fmla)
        object$coef <- co <- FIT$coef
        object$callstring <- callstring
        object$internal$fmla <- fmla
        if(is.finite(object$maxlogpl)) {
          ## Update maxlogpl provided it is finite
          ## (If the likelihood is infinite, this is due to the interaction;
          ## if we update the trend, the likelihood will remain infinite.)
          W <- glmdata$.mpl.W
          SUBSET <- glmdata$.mpl.SUBSET        
          Z <- is.data(object$Q)
          object$maxlogpl <- -(deviance(FIT)/2 +
                               sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
        }
        ## update the model call
        upd.call <- call
        upd.call$trend <- trend
        object$call <- upd.call
        ## update fitted interaction (depends on coefficients, if not Poisson)
        if(!is.null(inter <- object$interaction) && !is.poisson(inter)) 
          object$fitin <-
            fii(inter, co, object$internal$Vnames, object$internal$IsOffset)
        return(object)
      }
    }
  
    ## general case.
    undecided <- is.null(use.internal) || !is.logical(use.internal)
    force.int   <- !undecided && use.internal
    force.ext   <- !undecided && !use.internal
    if(!force.int) {
      ## check for validity of format
      badformat <- damaged.ppm(object)
    }
    if(undecided) {
      use.internal <- badformat
      if(badformat)
        message("object format corrupted; repairing it")
    } else if(force.ext && badformat)
      warning("object format corrupted; try update(object, use.internal=TRUE)")
  
    if(use.internal) {
      ## reset the main arguments in the call using the internal data
      call$Q <- data.ppm(object)
      namobj <- names(call)
      if("trend" %in% namobj)
        call$trend <- newformula(call$trend, object$trend, callframe, envir)
      if("interaction" %in% namobj) call$interaction <- object$interaction
      if("covariates" %in% namobj) call$covariates <- object$covariates
    }
  
    Q.is.new <- FALSE
  
    ## split named and unnamed arguments
    nama <- names(aargh)
    named <- if(is.null(nama)) rep.int(FALSE, length(aargh)) else (nama != "")
    namedargs <- aargh[named]
    unnamedargs <- aargh[!named]
    nama <- names(namedargs)
  
    if(any(named)) {
      ## any named arguments that were also present in the original call
      ## override their original values
      existing <- !is.na(match(nama, names(call)))
      for (a in nama[existing]) call[[a]] <- aargh[[a]]

      ## add any named arguments not present in the original call
      if (any(!existing)) {
        call <- c(as.list(call), namedargs[!existing])
        call <- as.call(call)
      }
      ## is the point pattern or quadscheme new ?
      if("Q" %in% nama)
        Q.is.new <- TRUE
    }
    if(any(!named)) {
      ## some objects identified by their class
      if(n <- sp.foundclasses(c("ppp", "quad"), unnamedargs, "Q", nama)) {
        call$Q <- unnamedargs[[n]]
        Q.is.new <- TRUE
      }
      if(n<- sp.foundclass("interact", unnamedargs, "interaction", nama))
        call$interaction <- unnamedargs[[n]]
      if(n <- sp.foundclasses(c("data.frame", "im"), unnamedargs, "covariates", nama))
        call$covariates <- unnamedargs[[n]]
      if(n<- sp.foundclass("formula", unnamedargs, "trend", nama))
        call$trend <- newformula(call$trend, unnamedargs[[n]], callframe, envir)
      else if(n <- sp.foundclass("character", unnamedargs, "trend", nama)) {
        ## string that might be interpreted as a formula
        strg <- unnamedargs[[n]]
        if(!is.na(charmatch("~", strg))) {
          fo <- as.formula(strg)
          call$trend <- newformula(call$trend, fo, callframe, envir)
        } 
      }
    }
  
    ## *************************************************************
    ## ****** Special action when Q is a point pattern *************
    ## *************************************************************
    if(Q.is.new && fixdummy && inherits((X <- eval(call$Q)), "ppp")) {
      ## Instead of allowing default.dummy(X) to occur,
      ## explicitly create a quadrature scheme from X,
      ## using the same dummy points and weight parameters
      ## as were used in the fitted model 
      Qold <- quad.ppm(object)
      if(is.marked(Qold)) {
        dpar <- Qold$param$dummy
        wpar <- Qold$param$weight
        Qnew <- do.call("quadscheme", append(list(X), append(dpar, wpar)))
      } else {
        Dum <- Qold$dummy
        wpar <- Qold$param$weight
        Qnew <- do.call("quadscheme", append(list(X, Dum), wpar))
      }
      ## replace X by new Q
      call$Q <- Qnew
    }

    ## finally call ppm
    return(eval(call, as.list(envir), enclos=callframe))
  }

  update.ppm
})

sp.foundclass <- function(cname, inlist, formalname, argsgiven) {
  ok <- unlist(lapply(inlist, inherits, what=cname))
  nok <- sum(ok)
  if(nok > 1)
    stop(paste("I am confused: there are two unnamed arguments",
               "of class", sQuote(cname)))
  if(nok == 0) return(0)
  absent <- !(formalname %in% argsgiven)
  if(!absent)
    stop(paste("I am confused: there is an unnamed argument",
               "of class", sQuote(cname), "which conflicts with the",
               "named argument", sQuote(formalname)))
  theposition <- seq_along(ok)[ok]
  return(theposition)
}

sp.foundclasses <- function(cnames, inlist, formalname, argsgiven) {
  ncn <- length(cnames)
  pozzie <- logical(ncn)
  for(i in seq_len(ncn))
    pozzie[i] <- sp.foundclass(cnames[i],  inlist, formalname, argsgiven)
  found <- (pozzie > 0)
  nfound <- sum(found)
  if(nfound == 0)
    return(0)
  else if(nfound == 1)
    return(pozzie[found])
  else
    stop(paste("I am confused: there are", nfound,
               "unnamed arguments of different classes (",
               paste(sQuote(cnames(pozzie[found])), collapse=", "),
               ") which could be interpreted as",
               sQuote(formalname)))
}
    

damaged.ppm <- function(object) {
  ## guess whether the object format has been damaged
  ## e.g. by dump/restore
  gf <- getglmfit(object)
  badfit <- !is.null(gf) && !inherits(gf$terms, "terms")
  if(badfit)
    return(TRUE)
  ## escape clause for fake models
  if(identical(object$fake, TRUE))
    return(FALSE)
  ## otherwise it was made by ppm
  Qcall <- object$call$Q
  cf <- object$callframe
  if(is.null(cf)) {
    ## Old format of ppm objects
    if(is.name(Qcall) && !exists(paste(Qcall)))
      return(TRUE)
    Q <- eval(Qcall)
  } else {
    ## New format of ppm objects
    if(is.name(Qcall) && !exists(paste(Qcall), cf))
      return(TRUE)
    Q <- eval(Qcall, cf)
  }
  badQ <- is.null(Q) || !(inherits(Q, "ppp") || inherits(Q,"quad"))
  return(badQ)
}
