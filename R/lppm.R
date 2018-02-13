#
#  lppm.R
#
#  Point process models on a linear network
#
#  $Revision: 1.42 $   $Date: 2018/02/13 02:41:14 $
#

lppm <- function(X, ...) {
  UseMethod("lppm")
}


lppm.formula <- function(X, interaction=NULL, ..., data=NULL) {
  ## remember call
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(X, "formula"))
    stop(paste("Argument 'X' should be a formula"))
  formula <- X
  
  if(spatstat.options("expand.polynom"))
    formula <- expand.polynom(formula)

  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- formula[[2L]]
  trend <- formula[c(1L,3L)]
  
  ## FIT #######################################
  thecall <- call("lppm", X=Yexpr, trend=trend,
                  data=data, interaction=interaction)
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
  result <- eval(thecall, parent.frame())

  result$call <- cl
  result$callstring <- callstring

  return(result)
}

lppm.lpp <- function(X, ..., eps=NULL, nd=1000, random=FALSE) {
  Xname <- short.deparse(substitute(X))
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  cl <- match.call()
  nama <- names(list(...))
  resv <- c("method", "forcefit")
  if(any(clash <- resv %in% nama))
    warning(paste(ngettext(sum(clash), "Argument", "Arguments"),
                  commasep(sQuote(resv[clash])),
                  "must not be used"))
  stopifnot(inherits(X, "lpp"))
  Q <- linequad(X, eps=eps, nd=nd, random=random)
  fit <- ppm(Q, ..., method="mpl", forcefit=TRUE)
  if(!is.poisson.ppm(fit))
    warning("Non-Poisson models currently use Euclidean distance")
  out <- list(X=X, fit=fit, Xname=Xname, call=cl, callstring=callstring)
  class(out) <- "lppm"
  return(out)
}

is.lppm <- function(x) { inherits(x, "lppm") }

# undocumented
as.ppm.lppm <- function(object) { object$fit }

fitted.lppm <- function(object, ..., dataonly=FALSE, new.coef=NULL,
                        leaveoneout=FALSE) {
  pfit <- object$fit
  v <- fitted(pfit, dataonly=dataonly, new.coef=new.coef,
              leaveoneout=leaveoneout)
  return(v)
}
  
predict.lppm <- function(object, ..., 
                         type="trend", locations=NULL,
                         new.coef=NULL) {
  type <- pickoption("type", type,
                     c(trend="trend", cif="cif", lambda="cif"))
  X <- object$X
  fit <- object$fit
  L <- as.linnet(X)

  if(!is.null(locations)) {
    #' locations given; return a vector of predicted values
    if(is.lpp(locations)) locations <- as.ppp(locations)
    values <- predict(fit, locations=locations, type=type, new.coef=new.coef)
    return(values)
  }
  
  # locations not given; want a pixel image
  # pixellate the lines
  Llines <- as.psp(L)
  linemask <- as.mask.psp(Llines, ...)
  lineimage <- as.im(linemask)
  # extract pixel centres
  xx <- rasterx.mask(linemask)
  yy <- rastery.mask(linemask)
  mm <- linemask$m
  xx <- as.vector(xx[mm])
  yy <- as.vector(yy[mm])
  pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=FALSE)
  pixdf <- data.frame(xc=xx, yc=yy)
  # project pixel centres onto lines
  p2s <- project2segment(pixelcentres, Llines)
  projloc <- as.data.frame(p2s$Xproj)
  projmap <- as.data.frame(p2s[c("mapXY", "tp")])
  projdata <- cbind(pixdf, projloc, projmap)
  # predict at the projected points
  if(!is.multitype(fit)) {
    values <- predict(fit, locations=projloc, type=type, new.coef=new.coef)
    # map to nearest pixels
    Z <- lineimage
    Z[pixelcentres] <- values
    # attach exact line position data
    df <- cbind(projdata, values)
    out <- linim(L, Z, df=df, restrict=FALSE)
  } else {
    # predict for each type
    lev <- levels(marks(data.ppm(fit)))
    out <- list()
    for(k in seq(length(lev))) {
      markk <- factor(lev[k], levels=lev)
      locnk <- cbind(projloc, data.frame(marks=markk))
      values <- predict(fit, locations=locnk, type=type, new.coef=new.coef)
      Z <- lineimage
      Z[pixelcentres] <- values
      df <- cbind(projdata, values)
      out[[k]] <- linim(L, Z, df=df, restrict=FALSE)
    }
    out <- as.solist(out)
    names(out) <- as.character(lev)
  }
  return(out)
}

coef.lppm <- function(object, ...) {
  coef(object$fit)
}

print.lppm <- function(x, ...) {
  splat("Point process model on linear network")
  print(x$fit)
  terselevel <- spatstat.options('terse')
  if(waxlyrical('extras', terselevel))
    splat("Original data:", x$Xname)
  if(waxlyrical('gory', terselevel))
    print(as.linnet(x))
  return(invisible(NULL))
}

summary.lppm <- function(object, ...) {
  splat("Point process model on linear network")
  print(summary(object$fit))
  terselevel <- spatstat.options('terse')
  if(waxlyrical('extras', terselevel))
    splat("Original data:", object$Xname)
  if(waxlyrical('gory', terselevel))
    print(summary(as.linnet(object)))
  return(invisible(NULL))
}

plot.lppm <- function(x, ..., type="trend") {
  xname <- short.deparse(substitute(x))
  y <- predict(x, type=type)
  do.call(plot, resolve.defaults(list(y),
                                   list(...),
                                   list(main=xname)))
}
  
anova.lppm <- function(object, ..., test=NULL) {
  stuff <- list(object=object, ...)
  if(!is.na(hit <- match("override", names(stuff)))) {
    warning("Argument 'override' is outdated and was ignored")
    stuff <- stuff[-hit]
  }
  #' extract ppm objects where appropriate
  mod <- sapply(stuff, is.lppm)
  stuff[mod] <- lapply(stuff[mod], getElement, name="fit")
  #' analysis of deviance or adjusted composite deviance
  do.call(anova.ppm, append(stuff, list(test=test)))
}

update.lppm <- function(object, ...) {
  stopifnot(inherits(object, "lppm"))
  X <- object$X
  fit <- object$fit
  Xname <- object$Xname
  callframe <- environment(formula(fit))
  aargh <- list(...)
  islpp <- sapply(aargh, is.lpp)
  if(any(islpp)) {
    # trap point pattern argument & convert to quadscheme
    ii <- which(islpp)
    if((npp <- length(ii)) > 1)
      stop(paste("Arguments not understood:", npp, "lpp objects given"))
    X <- aargh[[ii]]
    aargh[[ii]] <- linequad(X)
  }
  isfmla <- sapply(aargh, inherits, what="formula")
  if(any(isfmla)) {
    # trap formula pattern argument, update it, evaluate LHS if required
    jj <- which(isfmla)
    if((nf <- length(jj)) > 1)
      stop(paste("Arguments not understood:", nf, "formulae given"))
    fmla <- aargh[[jj]]
    fmla <- update(formula(object), fmla)
    if(!is.null(lhs <- lhs.of.formula(fmla))) {
      X <- eval(lhs, envir=list2env(list("."=X), parent=callframe))
      Qpos <- if(any(islpp)) ii else (length(aargh) + 1L)
      aargh[[Qpos]] <- linequad(X)
    }
    aargh[[jj]] <- rhs.of.formula(fmla)
  }
  newfit <- do.call(update.ppm,
                    append(list(fit), aargh),
                    envir=callframe)
  if(!is.poisson.ppm(newfit))
    warning("Non-Poisson models currently use Euclidean distance")
  out <- list(X=X, fit=newfit, Xname=Xname)
  class(out) <- "lppm"
  return(out)
}

terms.lppm <- function(x, ...) {
  terms(x$fit, ...)
}

logLik.lppm <- function(object, ...) {
  logLik(object$fit, ...)
}

deviance.lppm <- function(object, ...) {
  as.numeric(-2 * logLik(object, ...))
}

pseudoR2.lppm <- function(object, ..., keepoffset=TRUE) {
  dres <- deviance(object, ..., warn=FALSE)
  nullfmla <- . ~ 1
  if(keepoffset && has.offset.term(object)) {
    off <- attr(model.depends(object), "offset")
    offterms <- row.names(off)[apply(off, 1, any)]
    if(length(offterms)) {
      nullrhs <- paste(offterms, collapse=" + ") 
      nullfmla <- as.formula(paste(". ~ ", nullrhs))
    }
  } 
  nullmod <- update(object, nullfmla, forcefit=TRUE)
  dnul <- deviance(nullmod, warn=FALSE)
  return(1 - dres/dnul)
}

formula.lppm <- function(x, ...) {
  formula(x$fit, ...)
}

extractAIC.lppm <- function(fit, ...) {
  extractAIC(fit$fit, ...)
}

as.owin.lppm <- function(W, ..., fatal=TRUE) {
  stopifnot(inherits(W, "lppm"))
  as.owin(as.linnet(W), ..., fatal=fatal)
}

Window.lppm <- function(X, ...) { as.owin(X) }


model.images.lppm <- local({

  model.images.lppm <- function(object, L=as.linnet(object), ...) {
    stopifnot(inherits(object, "lppm"))
    stopifnot(inherits(L, "linnet"))
    m <- model.images(object$fit, W=as.rectangle(L), ...)
    if(length(m) > 0) {
      ## restrict images to L
      rasta <- as.mask(m[[1L]])
      DL <- as.mask.psp(as.psp(L), xy=rasta)
      ZL <- as.im(DL)
      if(!is.hyperframe) {
        ## list of images: convert to list of linims
        m <- tolinims(m, L=L, imL=ZL)
      } else {
        ## hyperframe, each column being a list of images
        mm <- lapply(as.list(m), tolinims, L=L, imL=ZL)
        m <- do.call(hyperframe, mm)
      }
    }
    return(m)
  }

  tolinim <- function(x, L, imL) linim(L, eval.im(x * imL), restrict=FALSE)
  tolinims <- function(x, L, imL) solapply(x, tolinim, L=L, imL=imL)
  
  model.images.lppm
})

  
model.matrix.lppm <- function(object,
                              data=model.frame(object, na.action=NULL),
                             ..., keepNA=TRUE) {
  stopifnot(is.lppm(object))
  if(missing(data)) data <- NULL
  model.matrix(object$fit, data=data, ..., keepNA=keepNA)
}

model.frame.lppm <- function(formula, ...) {
  stopifnot(inherits(formula, "lppm"))
  model.frame(formula$fit, ...)
}

domain.lppm <- as.linnet.lppm <- function(X, ...) {
  as.linnet(X$X, ...)
}

nobs.lppm <- function(object, ...) {
  npoints(object$X)
}

is.poisson.lppm <- function(x) { is.poisson(x$fit) }

is.stationary.lppm <- function(x) { is.stationary(x$fit) }

is.multitype.lppm <- function(X, ...) { is.multitype(X$fit) }

is.marked.lppm <- function(X, ...) { is.marked(X$fit) }

vcov.lppm <- function(object, ...) {
  if(!is.poisson(object))
    stop("vcov.lppm is only implemented for Poisson models")
  vcov(object$fit, ...)
}

valid.lppm <- function(object, ...) {
  valid(object$fit, ...)
}

emend.lppm <- function(object, ...) {
  object$fit <- emend(object$fit, ...)
  return(object)
}

  
