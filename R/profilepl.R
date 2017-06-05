#
# profilepl.R
#
#  $Revision: 1.43 $  $Date: 2017/06/05 10:31:58 $
#
#  computes profile log pseudolikelihood
#

profilepl <- local({

  ## Determine edge correction
  ## with partial matching, avoiding collisions with
  ## other arguments to ppm that have similar names.
  getppmcorrection <- function(..., correction = "border",
           covariates = NULL, covfunargs = NULL, control = NULL) {
    return(correction)
  }
  isSingleNA <- function(x) { length(x) == 1 && is.na(x) }
  
  profilepl <- function(s, f, ..., aic=FALSE, rbord=NULL, verbose=TRUE) {
    callenv <- parent.frame()
    s <- as.data.frame(s)
    n <- nrow(s)
    stopifnot(is.function(f))
    ## validate 's'
    parms <- names(s)
    fargs <- names(formals(f))
    if(!all(fargs %in% parms)) {
      bad <- !(fargs %in% parms)
      forgiven <- sapply(formals(f)[bad], isSingleNA)
      if(!all(forgiven)) {
        slecht <- fargs[bad[!forgiven]]
        nsl <- length(slecht)
        stop(paste(ngettext(nsl, "Argument", "Arguments"),
                   commasep(sQuote(slecht)),
                   ngettext(nsl, "is", "are"),
                   "not provided in the data frame s"))
      }
    }
    ## extra columns in 's' are assumed to be parameters of covariate functions
    is.farg <- parms %in% fargs
    pass.cfa <- any(!is.farg)
    got.cfa <- "covfunargs" %in% names(list(...))
    if(pass.cfa && got.cfa)
      stop("Some columns in s are superfluous")
    ##
    criterion <- numeric(n)
    ## make a fake call
    pseudocall <- match.call()
    pseudocall[[1]] <- as.symbol("ppm")
    namcal <- names(pseudocall)
    ## remove arguments 's' and 'verbose'
    retain <- !(namcal %in% c("s", "verbose"))
    pseudocall <- pseudocall[retain]
    namcal <- namcal[retain]
    ## place 'f' argument third 
    np <- length(pseudocall)
    fpos <- (1:np)[namcal == "f"]
    indices <- (1:np)[-fpos]
    if(length(indices) < 3) {
      indices <- c(indices, fpos)
    } else {
      indices <- c(indices[1:3], fpos, indices[-(1:3)])
    }
    pseudocall <- pseudocall[indices]
    namcal <- names(pseudocall)
    namcal[namcal=="f"] <- "interaction"
    names(pseudocall) <- namcal
    ## get correction
    correction <- getppmcorrection(...)
    if(correction == "border") {
      ## determine border correction distance
      if(is.null(rbord)) {
        ## compute rbord = max reach of interactions
        if(verbose) message("(computing rbord)")
        for(i in 1:n) {
          fi <- do.call(f, as.list(s[i, is.farg, drop=FALSE]))
          if(!inherits(fi, "interact"))
            stop(paste("f did not yield an object of class",
                       sQuote("interact")))
          re <- reach(fi)
          if(is.null(rbord))
            rbord <- re
          else if(rbord < re)
            rbord <- re
        }
      }
    } 
    ## determine whether computations can be saved
    if(pass.cfa || got.cfa) {
      savecomp <- FALSE
    } else {
      Q <- do.call(ppm,
                   append(list(...), list(rbord=rbord, justQ=TRUE)),
                   envir=callenv)
      savecomp <- !oversize.quad(Q)
    }
    ## go
    gc()
    if(verbose) {
      message(paste("comparing", n, "models..."))
      pstate <- list()
    }
    for(i in 1:n) {
      if(verbose)
        pstate <- progressreport(i, n, state=pstate)
      fi <- do.call(f, as.list(s[i, is.farg, drop=FALSE]))
      if(!inherits(fi, "interact"))
        stop(paste("f did not yield an object of class", sQuote("interact")))
      if(pass.cfa)
        cfai <- list(covfunargs=as.list(s[i, !is.farg, drop=FALSE])) 
      ## fit model
      if(i == 1) {
        ## fit from scratch
        arg1 <- list(...,
                     interaction=fi, 
                     rbord=rbord, savecomputed=savecomp,
                     warn.illegal=FALSE,
                     callstring="",
                     skip.border=TRUE)
        if(pass.cfa) arg1 <- append(arg1, cfai)
        fiti <- do.call(ppm, arg1, envir=callenv)
        ## save intermediate computations (pairwise distances, etc)
        precomp <- fiti$internal$computed
        savedargs <- list(...,
                          rbord=rbord, precomputed=precomp,
                          warn.illegal=FALSE,
                          callstring="",
                          skip.border=TRUE)
      } else {
        ## use precomputed data
        argi <- append(savedargs, list(interaction=fi))
        if(pass.cfa) argi <- append(argi, cfai)
        fiti <- do.call(ppm, argi, envir=callenv)
      }
      ## save log pl for each fit
      criterion[i] <-
          if(aic) -AIC(fiti) else as.numeric(logLik(fiti, warn=FALSE))
      ## save fitted coefficients for each fit
      co <- coef(fiti)
      if(i == 1) {
        allcoef <- data.frame(matrix(co, nrow=1))
        names(allcoef) <- names(co)
      } else
        allcoef <- rbind(allcoef, co)
    }
    if(verbose) message("fitting optimal model...")
    opti <- which.max(criterion)
    gc()
    optint <- do.call(f, as.list(s[opti, is.farg, drop=FALSE]))
    optarg <- list(..., interaction=optint, rbord=rbord)
    if(pass.cfa) {
      optcfa <- as.list(s[opti, !is.farg, drop=FALSE])
      attr(optcfa, "fitter") <- "profilepl"
      optarg <- append(optarg, list(covfunargs=optcfa))
    }
    optfit <- do.call(ppm, optarg, envir=callenv)
    if(verbose) message("done.")
    critname <- if(aic) "-AIC" else
                if(is.poisson(optfit)) "log l" else
                if(optfit$method == "logi") "log CL" else "log PL"
    result <- list(param=s,
                   prof=criterion,
                   critname=critname,
                   iopt=opti,
                   fit=optfit,
                   rbord=rbord,
                   fname=as.interact(optfit)$name,
                   allcoef=allcoef,
                   otherstuff=list(...),
                   pseudocall=pseudocall)
    class(result) <- c("profilepl", class(result))
    return(result)
  }

  profilepl
})

##
##   print method
##

print.profilepl <- function(x, ...) {
  head1 <- "profile log pseudolikelihood"
  head2 <- "for model: "
  psc <- paste(unlist(strsplitretain(format(x$pseudocall))),
               collapse=" ")
  if(nchar(psc) + nchar(head2) + 1 <= getOption('width')) {
    splat(head1)
    splat(head2, psc)
  } else {
    splat(head1, head2)
    splat(psc)
  }
  nparm <- ncol(x$param)
  if(waxlyrical('extras')) {
    corx <- x$fit$correction
    if(identical(corx, "border") && !is.null(x$rbord))
      splat("fitted with rbord =", x$rbord)
    splat("interaction:", x$fname)
    splat("irregular",
          ngettext(nparm, "parameter:", "parameters:\n"),
          paste(names(x$param),
                "in",
                unlist(lapply(lapply(as.list(x$param), range), prange)),
                collapse="\n"))
  }
  popt <- x$param[x$iopt,, drop=FALSE]
  splat("optimum",
        ngettext(nparm, "value", "values"),
        "of irregular",
        ngettext(nparm, "parameter: ", "parameters:\n"),
        commasep(paste(names(popt), "=", popt)))
  invisible(NULL)
}

##
##   summary method
##

summary.profilepl <- function(object, ...) {
  print(object)
  cat("\n\noptimal model:\n")
  print(object$fit)
}

as.ppm.profilepl <- function(object) {
  object$fit
}

predict.profilepl <- function(object, ...) {
  predict(object$fit, ...)
}

##
##  plot method 
##

plot.profilepl <- local({

  plot.profilepl <- function(x, ..., add=FALSE, main=NULL,
                             tag=TRUE, coeff=NULL, xvariable=NULL,
                             col=1, lty=1, lwd=1,
                             col.opt="green", lty.opt=3, lwd.opt=1) {
    para <- x$param
    ## graphics arguments may be expressions involving parameters
    if(ncol(para) > 1) {
      col <- eval(substitute(col), para)
      lwd <- eval(substitute(lwd), para)
      lty <- eval(substitute(lty), para)
      px <- cbind(para, col, lwd, lty, stringsAsFactors=FALSE)
      col <- px$col
      lwd <- px$lwd
      lty <- px$lty
    }
    ## strip any column that is entirely na
    nacol <- sapply(para, none.finite)
    para <- para[, !nacol, drop=FALSE]
    ## 
    npara <- ncol(para)
    ## main header
    if(is.null(main))
      main <- short.deparse(x$pseudocall)
    ## x variable for plot
    if(is.null(xvariable)) {
      xvalues <- para[,1]
      xname <- names(para)[1]
    } else {
      stopifnot(is.character(xvariable))
      if(!(xvariable %in% names(para)))
        stop("there is no irregular parameter named", sQuote(xvariable))
      xvalues <- para[[xvariable]]
      xname <- xvariable
    }
    ## y variable for plot                  
    if(is.null(coeff)) {
      yvalues <- x$prof
      ylab <- x$critname %orifnull% "log pl"
    } else {
      stopifnot(is.character(coeff))
      allcoef <- x$allcoef
      if(!(coeff %in% names(allcoef)))
        stop(paste("there is no coefficient named", sQuote(coeff),
                   "in the fitted model"))
      yvalues <- allcoef[[coeff]]
      ylab <- paste("coefficient:", coeff)
    }
    ## start plot
    if(!add)
      do.call.matched(plot.default,
                      resolve.defaults(list(x=range(xvalues), y=range(yvalues)),
                                       list(type="n", main=main),
                                       list(...),
                                       list(ylab=ylab, xlab=xname)),
                      extrargs=graphicsPars("plot"))

    linepars <- graphicsPars("lines")
  
    if(npara == 1) {
      ## single curve
      do.call.matched(lines.default,
                      resolve.defaults(list(x=xvalues, y=yvalues, ...),
                                       spatstat.options("par.fv")),
                      extrargs=linepars)
    } else {
      ## multiple curves
      other <- para[, -1, drop=FALSE]
      tapply(1:nrow(para),
             as.list(other),
             plotslice, 
             xvalues=xvalues, yvalues=yvalues, other=other,
             tag=tag, ...,
             col=col, lwd=lwd, lty=lty,
             lineargs=linepars)
    }

    
    ## show optimal value
    do.call.matched(abline,
                    resolve.defaults(list(v = xvalues[x$iopt]),
                                     list(...),
                                     list(lty=lty.opt, lwd=lwd.opt,
                                          col=col.opt)),
                    extrargs=linepars)
    return(invisible(NULL))
  }

  plotslice <- function(z, xvalues, yvalues, other, tag=TRUE, ...,
                        lty=1, col=1, lwd=1, lineargs) {
    fz <- xvalues[z]
    pz <- yvalues[z]
    n <- length(xvalues)
    if(length(lty) == n) lty <- unique(lty[z])[1]
    if(length(col) == n) col <- unique(col[z])[1]
    if(length(lwd) == n) lwd <- unique(lwd[z])[1]
    do.call.matched(lines.default,
                    resolve.defaults(list(x=fz, y=pz,
                                          col=col, lwd=lwd, lty=lty),
                                     list(...)),
                    extrargs=lineargs)
    if(tag) {
      oz <- other[z, , drop=FALSE]
      uniques <- apply(oz, 2, unique)
      labels <- paste(names(uniques), "=", uniques, sep="")
      label <- paste(labels, sep=",")
      ii <- which.max(pz)
      do.call.matched(text.default,
                      list(x=fz[ii], y=pz[ii], labels=label,
                           col=col, ...),
                      funargs=graphicsPars("text"))
    }
    return(NULL)
  }

  none.finite <- function(x) all(!is.finite(x))
  
  plot.profilepl
})


simulate.profilepl <- function(object, ...) {
  simulate(as.ppm(object), ...)
}

parameters.profilepl <- function(model, ...) {
  parameters(as.ppm(model))
}
