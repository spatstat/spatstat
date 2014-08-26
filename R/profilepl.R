#
# profilepl.R
#
#  $Revision: 1.22 $  $Date: 2014/08/07 02:25:27 $
#
#  computes profile log pseudolikelihood
#

profilepl <- function(s, f, ..., rbord=NULL, verbose=TRUE) {
  s <- as.data.frame(s)
  n <- nrow(s)
  fname <- paste(short.deparse(substitute(f)), collapse="")
  stopifnot(is.function(f))
  # validate 's'
  parms <- names(s)
  fargs <- names(formals(f))
  if(!all(fargs %in% parms))
    stop("Some arguments of f are not provided in s")
  # extra columns in 's' are assumed to be parameters of covariate functions
  is.farg <- parms %in% fargs
  pass.cfa <- any(!is.farg)
  got.cfa <- "covfunargs" %in% names(list(...))
  if(pass.cfa && got.cfa)
    stop("Some columns in s are superfluous")
  #
  interlist <- list()
  logmpl <- numeric(n)
  # make a fake call
  pseudocall <- match.call()
  pseudocall[[1]] <- as.symbol("ppm")
  namcal <- names(pseudocall)
  # remove arguments 's' and 'verbose'
  retain <- !(namcal %in% c("s", "verbose"))
  pseudocall <- pseudocall[retain]
  namcal <- namcal[retain]
  # place 'f' argument third 
  np <- length(pseudocall)
  fpos <- (1:np)[namcal == "f"]
  indices <- (1:np)[-fpos]
  if(length(indices) < 3)
    indices <- c(indices, fpos)
  else 
    indices <- c(indices[1:3], fpos, indices[-(1:3)])
  pseudocall <- pseudocall[indices]
  namcal <- names(pseudocall)
  namcal[namcal=="f"] <- "interaction"
  names(pseudocall) <- namcal
  # Determine edge correction
  # with partial matching, avoiding collisions with
  # other arguments to ppm that have similar names.
  getppmcorrection <- function(..., correction = "border",
           covariates = NULL, covfunargs = NULL, control = NULL) {
    return(correction)
  }
  correction <- getppmcorrection(...)
  if(correction == "border") {
    # determine border correction distance
    if(is.null(rbord)) {
      # compute rbord = max reach of interactions
      if(verbose) message("(computing rbord)")
      for(i in 1:n) {
        fi <- do.call("f", as.list(s[i, is.farg, drop=FALSE]))
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
  # determine whether computations can be saved
  if(pass.cfa || got.cfa) {
    savecomp <- FALSE
  } else {
    Q <- ppm(..., rbord=rbord, justQ=TRUE)
    savecomp <- !oversize.quad(Q)
  }
  
  # fit one model and extract quadscheme
  if(verbose) message(paste("Comparing", n, "models..."))
  for(i in 1:n) {
    if(verbose)
      progressreport(i, n)
    fi <- do.call("f", as.list(s[i, is.farg, drop=FALSE]))
    if(!inherits(fi, "interact"))
      stop(paste("f did not yield an object of class", sQuote("interact")))
    if(pass.cfa)
      cfai <- list(covfunargs=as.list(s[i, !is.farg, drop=FALSE])) 
    # fit model
    if(i == 1) {
      # fit from scratch
      arg1 <- list(interaction=fi, ...,
                   rbord=rbord, savecomputed=savecomp,
                   warn.illegal=FALSE,
                   callstring="",
                   skip.border=TRUE)
      if(pass.cfa) arg1 <- append(arg1, cfai)
      fiti <- do.call("ppm", arg1)
      # save intermediate computations (pairwise distances, etc)
      precomp <- fiti$internal$computed
      savedargs <- list(...,
                        rbord=rbord, precomputed=precomp,
                        warn.illegal=FALSE,
                        callstring="",
                        skip.border=TRUE)
    } else {
      # use precomputed data
      argi <- append(savedargs, list(interaction=fi))
      if(pass.cfa) argi <- append(argi, cfai)
      fiti <- do.call("ppm", argi)
    }
    # save log PL for each fit
    logmpl[i] <- as.numeric(logLik(fiti, warn=FALSE))
    # save fitted coefficients for each fit
    co <- coef(fiti)
    if(i == 1) {
      allcoef <- data.frame(matrix(co, nrow=1))
      names(allcoef) <- names(co)
    } else
    allcoef <- rbind(allcoef, co)
  }
  if(verbose) message("Fitting optimal model...")
  opti <- which.max(logmpl)
  optint <- do.call("f", as.list(s[opti, is.farg, drop=FALSE]))
  optarg <- list(interaction=optint, ..., rbord=rbord)
  if(pass.cfa) {
    optcfa <- as.list(s[opti, !is.farg, drop=FALSE])
    attr(optcfa, "fitter") <- "profilepl"
    optarg <- append(optarg, list(covfunargs=optcfa))
  }
  optfit <- do.call("ppm", optarg)
  if(verbose) message("done.")
  result <- list(param=s,
                 prof=logmpl,
                 iopt=opti,
                 fit=optfit,
                 rbord=rbord,
                 fname=fname,
                 allcoef=allcoef,
                 otherstuff=list(...),
                 pseudocall=pseudocall)
  class(result) <- c("profilepl", class(result))
  return(result)
}

#
#   print method
#

print.profilepl <- function(x, ...) {
  head1 <- "Profile log pseudolikelihood"
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
    splat("fitted with rbord =", x$rbord)
    splat("Interaction:", x$fname)
    splat("Irregular",
          ngettext(nparm, "parameter:", "parameters:\n"),
          paste(names(x$param),
                "in",
                unlist(lapply(lapply(as.list(x$param), range), prange)),
                collapse="\n"))
  }
  popt <- x$param[x$iopt,, drop=FALSE]
  splat("Optimum",
        ngettext(nparm, "value", "values"),
        "of irregular",
        ngettext(nparm, "parameter: ", "parameters:\n"),
        commasep(paste(names(popt), "=", popt)))
  invisible(NULL)
}

#
#   summary method
#

summary.profilepl <- function(object, ...) {
  print(object)
  cat("\n\nOptimal model:\n")
  print(object$fit)
}

as.ppm.profilepl <- function(object) {
  object$fit
}


#
#  plot method 
#

plot.profilepl <- function(x, ..., add=FALSE, main=NULL,
                           tag=TRUE, coeff=NULL, xvariable=NULL) {
  para <- x$param
  npara <- ncol(para)
  # main header
  if(is.null(main))
    main <- short.deparse(x$pseudocall)
  # x variable for plot
  if(is.null(xvariable)) {
    xvalues <- para[,1]
    xname <- names(para)[1]
  } else {
    stopifnot(is.character(xvariable))
    if(!(xvariable %in% names(para)))
      stop("There is no irregular parameter named", sQuote(xvariable))
    xvalues <- para[[xvariable]]
    xname <- xvariable
  }
  # y variable for plot                  
  if(is.null(coeff)) {
    yvalues <- x$prof
    ylab <- "log PL"
  } else {
    stopifnot(is.character(coeff))
    allcoef <- x$allcoef
    if(!(coeff %in% names(allcoef)))
      stop(paste("There is no coefficient named", sQuote(coeff), "in the fitted model"))
    yvalues <- allcoef[[coeff]]
    ylab <- paste("coefficient:", coeff)
  }
  # start plot
  if(!add)
    do.call.matched("plot.default",
                  resolve.defaults(list(x=range(xvalues), y=range(yvalues)),
                                   list(type="n", main=main),
                                   list(...),
                                   list(ylab=ylab, xlab=xname)))
  # single curve
  if(npara == 1) {
    do.call.matched("lines.default", list(x=xvalues, y=yvalues, ...))
    abline(v = xvalues[x$iopt], lty=2, col="green")
    return(invisible(NULL))
  }

  # multiple curves
  other <- para[, -1, drop=FALSE]
  tapply(1:nrow(para),
         as.list(other),
         function(z, xvalues, yvalues, other, tag) {
           fz <- xvalues[z]
           pz<- yvalues[z]
           lines(fz, pz)
           if(tag) {
             oz <- other[z, , drop=FALSE]
             uniques <- apply(oz, 2, unique)
             labels <- paste(names(uniques), "=", uniques, sep="")
             label <- paste(labels, sep=",")
             ii <- which.max(pz)
             text(fz[ii], pz[ii], label)
           }
         },
         xvalues=xvalues, yvalues=yvalues, other=other, tag=tag)

  abline(v = xvalues[x$iopt], lty=2, col="green")
  return(invisible(NULL))
}

