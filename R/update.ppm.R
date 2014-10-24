#
#  update.ppm.R
#
#
#  $Revision: 1.55 $    $Date: 2014/10/24 00:22:30 $
#
#
#

update.ppm <- local({

  ## update formula and expand polynomial
  newformula <- function(old, change, eold, enew) {
    old <- if(is.null(old)) ~1 else eval(old, eold)
    change <- if(is.null(change)) ~1 else eval(change, enew)
    old <- as.formula(old, env=eold)
    change <- as.formula(change, env=enew)
    answer <- update.formula(old, change)
    if(spatstat.options("expand.polynom")) 
      answer <- expand.polynom(answer)
    return(answer)
  }

  ## update point pattern dataset using either data or formula
  newpattern <- function(oldpattern, lhs, callframe, envir) {
    eval(eval(substitute(substitute(l, list("."=Q)),
                         list(l=lhs,
                              Q=oldpattern)),
              envir=as.list(envir), enclos=callframe),
         envir=as.list(envir), enclos=callframe)
  }
  
  update.ppm <- function(object, ...,
                         fixdummy=TRUE, use.internal=NULL,
                         envir=environment(terms(object))) {
    verifyclass(object, "ppm")
    new.callstring <- short.deparse(sys.call())
    aargh <- list(...)
    
    if(inherits(object, "ippm")) {
      call <- object$dispatched$call
      callframe <- object$dispatched$callframe
    } else {
      call <- getCall(object)
      if(!is.call(call))
        stop(paste("Internal error - getCall(object) is not of class",
                   sQuote("call")))
      callframe <- object$callframe
    }
    
    callfun <- as.character(call[[1]])
    newstyle <- (callfun == "ppm.formula")
    oldstyle <- !newstyle
    
    ## Special cases 
    ## (1) no new information given
    if(length(aargh) == 0 && !identical(use.internal, TRUE)) {
      result <- eval(call, as.list(envir), enclos=callframe)
      result$callframe <- callframe
      return(result)
    }

    ## (2) model can be updated using existing covariate data frame
    if(!identical(use.internal, FALSE) &&
       ## single argument which is a formula
       (length(aargh) == 1) &&
       inherits(fmla <- aargh[[1]], "formula") &&
       is.null(lhs.of.formula(fmla)) &&
       ## not a ppm.formula call
       oldstyle &&
       ## fitted by mpl using glm/gam
       with(object,
            method == "mpl" &&
            !is.null(fitter) &&
            fitter %in% c("gam", "glm"))) {
      ## This is a dangerous hack! 
      glmdata <- object$internal$glmdata
      ## check whether data for new variables are available
      ## (this doesn't work with things like 'pi')
      vars.available <- c(colnames(glmdata), names(object$covfunargs))
      if(all(variablesinformula(fmla) %in% c(".", vars.available))) {
        ## we can update using internal data
        FIT <- object$internal$glmfit
        orig.env <- environment(FIT$terms)
        ## update formulae using "." rules
        trend <- newformula(object$trend, fmla, callframe, envir)
        fmla  <- newformula(formula(FIT), fmla, callframe, envir)
        ## expand polynom() in formula
        if(spatstat.options("expand.polynom")) {
          fmla <- expand.polynom(fmla)
          trend <- expand.polynom(trend)
        }
        ## update GLM/GAM fit 
        upd.glm.call <- update(FIT, fmla, evaluate=FALSE)
        FIT <- eval(upd.glm.call, envir=orig.env)
        environment(FIT$terms) <- orig.env
        object$internal$glmfit <- FIT
        ## update entries of object
        object$trend <- trend
        object$terms <- terms(fmla)
        object$coef <- co <- FIT$coef
        object$callstring <- new.callstring
        object$internal$fmla <- fmla
        ##
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
        ##
        if(is.stationary(object) && !is.marked(object)) {
          ## uniform Poisson
          if(eval(call$rename.intercept) %orifnull% TRUE) {
             names(object$coef) <- "log(lambda)"
          }
        }
        return(object)
      }
    }

    ## (3) Need to use internal data   
    if(oldstyle) {
      ## decide whether to use internal data
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
        call$Q <- quad.ppm(object)
        namobj <- names(call)
        if("trend" %in% namobj)
          call$trend <- newformula(call$trend, object$trend, callframe, envir)
        if("interaction" %in% namobj) call$interaction <- object$interaction
        if("covariates" %in% namobj) call$covariates <- object$covariates
      }
    }

    ## General case.
    X.is.new <- FALSE
    
    ## First split named and unnamed arguments
    nama <- names(aargh)
    named <- if(is.null(nama)) rep.int(FALSE, length(aargh)) else nzchar(nama)
    namedargs <- aargh[named]
    unnamedargs <- aargh[!named]
    nama <- names(namedargs)

    ## Find the argument 'Q' by name or implicitly by class
    ##   (including detection of conflicts)
    argQ <- NULL
    if(n <- sp.foundclasses(c("ppp", "quad"), unnamedargs, "Q", nama)) {
      argQ <- unnamedargs[[n]]
      unnamedargs <- unnamedargs[-n]
    }
    if("Q" %in% nama) {
      argQ <- namedargs$Q
      nama <- setdiff(nama, "Q")
      namedargs <- namedargs[nama]
    }
    ## Deal with argument 'Q' which has several possible forms
    if(!is.null(argQ)) {
      X.is.new <- TRUE
      if(inherits(argQ, "formula")) {
        ## Q = X ~ trend
        if(newstyle) {
          ## update the formula
          call$Q <- newformula(call$Q, argQ, callframe, envir)
        } else {
          ## split into Q = X and trend = ~trend
          if(!is.null(lhs <- lhs.of.formula(argQ)))
            call$Q <- newpattern(call$Q, lhs, callframe, envir)
          call$trend <- newformula(call$trend,
                                   rhs.of.formula(eval(argQ)),
                                   callframe, envir)
        }
      } else {
        ## Q = X
        if(newstyle) {
          ## convert old call to old style
          fo <- as.formula(call$Q)
          Yexpr <- lhs.of.formula(fo)
          trend <- rhs.of.formula(fo)
          newcall <- call("ppm", Q=Yexpr, trend=trend)
          if(length(call) > 2) {
            whichQ <- which(names(call) == "Q")
            morecall <- call[-c(1, whichQ)]
            if((mc <- length(morecall)) > 0) {
              newcall[3 + 1:mc] <- morecall
              names(newcall)[3 + 1:mc] <- names(call)[-c(1, whichQ)]
            }
          }
          call <- newcall
          newstyle <- FALSE
          oldstyle <- TRUE
        }
        ## Now update the dataset
        call$Q <- argQ
      }
    }

    ## Find any formula arguments 
    ##   (including detection of conflicts)
    argfmla <- NULL
    if(n <- sp.foundclass("formula", unnamedargs, "trend", nama)) {
      argfmla <- unnamedargs[[n]]
      unnamedargs <- unnamedargs[-n]
    } else if(n <- sp.foundclass("character", unnamedargs, "trend", nama)) {
      ## string that might be interpreted as a formula
      strg <- unnamedargs[[n]]
      if(!is.na(charmatch("~", strg))) {
        argfmla <- as.formula(strg)
        unnamedargs <- unnamedargs[-n]
      }
    }
    if("trend" %in% nama) {
      argfmla <- namedargs$trend
      nama <- setdiff(nama, "trend")
      namedargs <- namedargs[nama]
    }
    ## Handle new formula
    if(!is.null(argfmla)) {
      lhs <- lhs.of.formula(argfmla)
      if(newstyle) {
        ## ppm.formula: update the formula
        if(is.null(lhs)) {
          argfmla <- as.formula(paste(".", deparse(argfmla)))
        } else X.is.new <- TRUE
        call$Q <- newformula(call$Q, argfmla, callframe, envir)
      } else {
        ## ppm.ppp: update the trend and possibly the data
        if(is.null(lhs)) {
          ## assign new trend
          call$trend <- newformula(call$trend, argfmla, callframe, envir)
        } else {
          ## split into Q = X and trend = ~trend
          X.is.new <- TRUE
          call$Q <- newpattern(call$Q, lhs, callframe, envir)
          call$trend <- newformula(call$trend,
                                   rhs.of.formula(argfmla),
                                   callframe, envir)
        }
      } 
    }
    
    if(length(namedargs) > 0) {
      ## any other named arguments that were also present in the original call
      ## override their original values.
      existing <- !is.na(match(nama, names(call)))
      for (a in nama[existing]) call[[a]] <- aargh[[a]]

      ## add any named arguments not present in the original call
      if (any(!existing)) {
        call <- c(as.list(call), namedargs[!existing])
        call <- as.call(call)
      }
    }
    if(length(unnamedargs) > 0) {
      ## some further objects identified by their class
      if(n<- sp.foundclass("interact", unnamedargs, "interaction", nama)) {
        call$interaction <- unnamedargs[[n]]
        unnamedargs <- unnamedargs[-n]
      }
      if(n <- sp.foundclasses(c("data.frame", "im"),
                              unnamedargs, "covariates", nama)) {
        call$covariates <- unnamedargs[[n]]
        unnamedargs <- unnamedargs[-n]
      }
    }
  
    ## *************************************************************
    ## ****** Special action when Q is a point pattern *************
    ## *************************************************************
    if(X.is.new && fixdummy && oldstyle &&
       inherits((X <- eval(call$Q, as.list(envir), enclos=callframe)), "ppp")) {
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
    call[[1]] <- as.name('ppm')
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
  badQ <- is.null(Q) || !(inherits(Q, c("ppp", "quad", "formula")))
  return(badQ)
}
