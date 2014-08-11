#
#  lppm.R
#
#  Point process models on a linear network
#
#  $Revision: 1.21 $   $Date: 2013/09/10 09:57:02 $
#

lppm <- function(X, ...) {
  UseMethod("lppm")
}


lppm.formula <- function(X, interaction=NULL, ..., data=NULL) {
  ## remember call
  callstring <- short.deparse(sys.call())
  ## cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(X, "formula"))
    stop(paste("Argument 'X' should be a formula"))
  formula <- X
  
  if(spatstat.options("expand.polynom"))
    formula <- expand.polynom(formula)

  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- lhs <- formula[[2]]
  trend <- rhs <- formula[c(1,3)]
  
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

  if(!("callstring" %in% names(list(...))))
    result$callstring <- callstring
  
  return(result)
}

lppm.lpp <- function(X, ..., eps=NULL, nd=1000) {
  Xname <- short.deparse(substitute(X))
  nama <- names(list(...))
  resv <- c("method", "forcefit")
  if(any(clash <- resv %in% nama))
    warning(paste(ngettext(sum(clash), "Argument", "Arguments"),
                  commasep(sQuote(resv[clash])),
                  "must not be used"))
  stopifnot(inherits(X, "lpp"))
  Q <- linequad(X, eps=eps, nd=nd)
  fit <- ppm(Q, ..., method="mpl", forcefit=TRUE)
  if(!is.poisson.ppm(fit))
    warning("Non-Poisson models currently use Euclidean distance")
  out <- list(X=X, fit=fit, Xname=Xname)
  class(out) <- "lppm"
  return(out)
}

is.lppm <- function(x) { inherits(x, "lppm") }

predict.lppm <- function(object, ..., 
                         type="trend", locations=NULL) {
  type <- pickoption("type", type,
                     c(trend="trend", cif="cif", lambda="cif"))
  X <- object$X
  fit <- object$fit
  L <- as.linnet(X)

  if(!is.null(locations)) {
    # locations given; return a vector of predicted values
    values <- predict(fit, locations=locations, type=type)
    return(values)
  }
  
  # locations not given; want a pixel image
  # pixellate the lines
  Llines <- as.psp(L)
  linemask <- as.mask.psp(Llines, ...)
  lineimage <- as.im(linemask)
  # extract pixel centres
  xx <- raster.x(linemask)
  yy <- raster.y(linemask)
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
    values <- predict(fit, locations=projloc, type=type)
    # map to nearest pixels
    Z <- lineimage
    Z[pixelcentres] <- values
    # attach exact line position data
    df <- cbind(projdata, values)
    out <- linim(L, Z, df=df)
  } else {
    # predict for each type
    lev <- levels(marks(data.ppm(fit)))
    out <- list()
    for(k in seq(length(lev))) {
      markk <- factor(lev[k], levels=lev)
      locnk <- cbind(projloc, data.frame(marks=markk))
      values <- predict(fit, locations=locnk, type=type)
      Z <- lineimage
      Z[pixelcentres] <- values
      df <- cbind(projdata, values)
      out[[k]] <- linim(L, Z, df=df)
    }
    names(out) <- as.character(lev)
    class(out) <- c("listof", class(out))
  }
  return(out)
}

coef.lppm <- function(object, ...) {
  coef(object$fit)
}

print.lppm <- function(x, ...) {
  cat("Point process model on linear network\n")
  print(x$fit)
  cat("Linear network:\n")
  print(as.linnet(x))
  cat(paste("Original data:", x$Xname, "\n"))
  return(invisible(NULL))
}

plot.lppm <- function(x, ..., type="trend") {
  xname <- short.deparse(substitute(x))
  y <- predict(x, type=type)
  do.call("plot", resolve.defaults(list(y),
                                   list(...),
                                   list(main=xname)))
}
  
anova.lppm <- function(object, ..., test=NULL, override=FALSE) {
  stuff <- list(object=object, ...)
  # extract ppm objects where appropriate
  stuff <- lapply(stuff, function(z) { if(inherits(z, "lppm")) z$fit else z })
  # analysis of deviance for 
  do.call("anova.ppm", append(stuff, list(test=test, override=override)))
}

update.lppm <- function(object, ...) {
  stopifnot(inherits(object, "lppm"))
  X <- object$X
  fit <- object$fit
  Xname <- object$Xname
  aargh <- list(...)
  islpp <- unlist(lapply(aargh, inherits, what="lpp"))
  if(!any(islpp)) {
    # pass arguments through to update.ppm
    newfit <- do.call("update.ppm", append(list(fit), aargh))
    newX <- X
  } else {
    # trap point pattern argument & convert to quadscheme
    if((npp <- sum(islpp)) > 1)
      stop(paste("Arguments not understood:", npp, "lpp objects given"))
    newX <- aargh[[which(islpp)]]
    newQ <- linequad(newX)
    newfit <- do.call("update.ppm",
                      append(list(fit, newQ), aargh[!islpp]))
  } 
  if(!is.poisson.ppm(newfit))
    warning("Non-Poisson models currently use Euclidean distance")
  out <- list(X=newX, fit=newfit, Xname=Xname)
  class(out) <- "lppm"
  return(out)
}

terms.lppm <- function(x, ...) {
  terms(x$fit, ...)
}

logLik.lppm <- function(object, ...) {
  logLik(object$fit, ...)
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

model.images.lppm <- function(object, L=as.linnet(object), ...) {
  stopifnot(inherits(object, "lppm"))
  stopifnot(inherits(L, "linnet"))
  m <- model.images(object$fit, W=as.rectangle(L), ...)
  if(length(m) > 0) {
    # restrict images to L
    rasta <- as.mask(m[[1]])
    DL <- as.mask.psp(as.psp(L), xy=rasta)
    ZL <- as.im(DL)
    m <- lapply(m, function(x, Z) eval.im(x * Z), Z=ZL)
    # convert to linim
    m <- lapply(m, function(x, L) linim(L,x), L=L)
  }
  return(as.listof(m))
}
  
model.matrix.lppm <- function(object, data=model.frame(object),
                             ..., keepNA=TRUE) {
  stopifnot(inherits(object, "lppm"))
  model.matrix(object$fit, data=data, ..., keepNA=keepNA)
}

model.frame.lppm <- function(formula, ...) {
  stopifnot(inherits(formula, "lppm"))
  model.frame(formula$fit, ...)
}

as.linnet.lppm <- function(X, ...) {
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
    stop("vov.lppm is only implemented for Poisson models")
  vcov(object$fit, ...)
}
