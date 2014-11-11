#
#   rho2hat.R
#
#   Relative risk for pairs of covariate values
#
#   $Revision: 1.19 $   $Date: 2014/11/11 02:43:31 $
#

rho2hat <- function(object, cov1, cov2, ..., method=c("ratio", "reweight")) {
  cov1name <- short.deparse(substitute(cov1))
  cov2name <- short.deparse(substitute(cov2))
  callstring <- short.deparse(sys.call())
  method <- match.arg(method)
  # validate model
  if(is.ppp(object) || inherits(object, "quad")) {
    model <- ppm(object, ~1, forcefit=TRUE)
    reference <- "area"
    modelcall <- NULL
  } else if(is.ppm(object)) {
    model <- object
    reference <- "model"
    modelcall <- model$call
    if(is.null(getglmfit(model)))
      model <- update(model, forcefit=TRUE)
  } else stop("object should be a point pattern or a point process model")

  # interpret string "x" or "y" as a coordinate function
  getxyfun <- function(s) {
    switch(s,
           x = { function(x,y) { x } },
           y = { function(x,y) { y } },
           stop(paste("Unrecognised covariate name", sQuote(s))))
  }
  if(is.character(cov1) && length(cov1) == 1) {
    cov1name <- cov1
    cov1 <- getxyfun(cov1name)
  }
  if(is.character(cov2) && length(cov2) == 1) {
    cov2name <- cov2
    cov2 <- getxyfun(cov2name)
  }
  if(   (cov1name == "x" && cov2name == "y")
     || (cov1name == "y" && cov2name == "x")) {
    # spatial relative risk
    isxy <- TRUE
    needflip <- (cov1name == "y" && cov2name == "x")
    X <- data.ppm(model)
    if(needflip) X <- flipxy(X)
    switch(method,
           ratio = {
             # ratio of smoothed intensity estimates
             den <- density(X, ...)
             sigma <- attr(den, "sigma")
             varcov <- attr(den, "varcov")
             W <- as.owin(den)
             rslt <- switch(reference,
                            area = { den },
                            model = {
                              lam <- predict(model, locations=W)
                              if(needflip) lam <- flipxy(lam)
                              lam <- blur(lam, sigma=sigma, varcov=varcov,
                                          normalise=TRUE)
                              eval.im(den/lam)
                            })
           },
           reweight = {
             # smoothed point pattern with weights = 1/reference
             W <- do.call.matched("as.mask",
                                  append(list(w=as.owin(X)), list(...)))
             gstarX <- switch(reference,
                              area = {
                                rep.int(area(W), npoints(X))
                              },
                              model = {
                                lam <- predict(model, locations=W)
                                if(needflip) lam <- flipxy(lam)
                                lam[X]
                              })
             rslt <- density(X, ..., weights=1/gstarX)
             sigma <- attr(rslt, "sigma")
             varcov <- attr(rslt, "varcov")
           })
    Z12points <- X
    r1 <- W$xrange
    r2 <- W$yrange
  } else {
    # general case
    isxy <- FALSE
    # harmonise covariates 
    if(is.function(cov1) && is.im(cov2)) {
      cov1 <- as.im(cov1, W=cov2)
    } else if(is.im(cov1) && is.function(cov2)) {
      cov2 <- as.im(cov2, W=cov1)
    }
    # evaluate each covariate at data points and at pixels
    stuff1 <- evalCovar(model, cov1)
    stuff2 <- evalCovar(model, cov2)
    # unpack
    values1 <- stuff1$values
    values2 <- stuff2$values
    # covariate values at each data point
    Z1X      <- values1$ZX
    Z2X      <- values2$ZX
    # covariate values at each pixel
    Z1values <- values1$Zvalues
    Z2values <- values2$Zvalues
    # model intensity
    lambda  <- values1$lambda
    # ranges of each covariate
    r1 <- range(Z1X, Z1values, finite=TRUE)
    r2 <- range(Z2X, Z2values, finite=TRUE)
    scal <- function(x, r) { (x - r[1])/diff(r) }
    # scatterplot coordinates
    Z12points <- ppp(scal(Z1X, r1), scal(Z2X, r2), c(0,1), c(0,1))
    Z12pixels <- ppp(scal(Z1values, r1), scal(Z2values, r2), c(0,1), c(0,1))
    # normalising constants
#    nX   <- length(Z1X)
    npixel <- length(lambda)
    areaW <- area(Window(model))
    pixelarea <- areaW/npixel
    baseline <- if(reference == "area") rep.int(1, npixel) else lambda
    wts <- baseline * pixelarea
    switch(method,
           ratio = {
             # estimate intensities
             fhat <- density(Z12points, ...)
             sigma <- attr(fhat, "sigma")
             varcov <- attr(fhat, "varcov")
             ghat <- do.call("density.ppp",
                             resolve.defaults(list(Z12pixels, weights=wts),
                                              list(...),
                                              list(sigma=sigma,
                                                   varcov=varcov)))
             # compute ratio of smoothed densities
             rslt <- eval.im(fhat/ghat)
           },
           reweight = {
             # compute smoothed intensity with weight = 1/reference
             ghat <- density(Z12pixels, weights=wts, ...)
             rslt <- density(Z12points, weights=1/ghat[Z12points], ...)
             sigma <- attr(rslt, "sigma")
             varcov <- attr(rslt, "varcov")
           })
  }
  # add scale and label info
  attr(rslt, "stuff") <- list(isxy=isxy,
                              cov1name=cov1name,
                              cov2name=cov2name,
                              r1=r1,
                              r2=r2,
                              reference=reference,
                              modelcall=modelcall,
                              callstring=callstring,
                              Z12points=Z12points,
                              sigma=sigma,
                              varcov=varcov)
  class(rslt) <- c("rho2hat", class(rslt))
  rslt
}

plot.rho2hat <- function(x, ..., do.points=FALSE) {
  xname <- short.deparse(substitute(x))
  s <- attr(x, "stuff")
  # resolve "..." arguments
  rd <- resolve.defaults(list(...),
                         list(add=FALSE, axes=!s$isxy,
                              xlab=s$cov1name, ylab=s$cov2name))
  # plot image
  plotparams <- graphicsPars("plot")
  do.call.matched("plot.im",
                  resolve.defaults(list(x=x, axes=FALSE),
                                   list(...), list(main=xname)),
                  extrargs=c(plotparams, "add", "zlim", "breaks"))
  # add axes 
  if(rd$axes) {
    axisparams <- graphicsPars("axis")
    Axis <- function(..., extrargs=axisparams) {
      do.call.matched("axis", resolve.defaults(list(...)), extrargs=extrargs)
    }
    if(s$isxy) {
      # for (x,y) plots the image is at the correct physical scale
      xr <- x$xrange
      yr <- x$yrange
      spak <- 0.05 * max(diff(xr), diff(yr))
      Axis(side=1, ..., at=pretty(xr), pos=yr[1] - spak)
      Axis(side=2, ..., at=pretty(yr), pos=xr[1] - spak)
    } else {
      # for other plots the image was scaled to the unit square
      rx <- s$r1
      ry <- s$r2
      px <- pretty(rx)
      py <- pretty(ry)
      Axis(side=1, labels=px, at=(px - rx[1])/diff(rx), ...)
      Axis(side=2, labels=py, at=(py - ry[1])/diff(ry), ...)
    }
    title(xlab=rd$xlab)
    title(ylab=rd$ylab)
  }
  if(do.points) {
    do.call.matched("plot.ppp",
                    resolve.defaults(list(x=s$Z12points, add=TRUE),
                                     list(...)),
                    extrargs=c("pch", "col", "cols", "bg", "cex", "lwd", "lty"))
  }
  invisible(NULL)
}

print.rho2hat <- function(x, ...) {
  s <- attr(x, "stuff")
  cat("Scatterplot intensity estimate (class rho2hat)\n")
  cat(paste("for the covariates", s$cov1name, "and", s$cov2name, "\n"))
  switch(s$reference,
         area=cat("Function values are absolute intensities\n"),
         model={
           cat("Function values are relative to fitted model\n")
           print(s$modelcall)
         })
  cat(paste("Call:", s$callstring, "\n"))
  if(s$isxy) {
    cat("Obtained by spatial smoothing of original data\n")
    cat("Smoothing parameters used by density.ppp:\n")
  } else {
    cat("Obtained by transforming to the unit square and smoothing\n")
    cat("Smoothing parameters (on unit square) used by density.ppp:\n")
  }
  if(!is.null(s$sigma)) cat(paste("\tsigma = ", signif(s$sigma, 5), "\n"))
  if(!is.null(s$varcov)) { cat("\tvarcov =\n") ; print(s$varcov) }
  cat("Intensity values:\n")
  NextMethod("print")
}
