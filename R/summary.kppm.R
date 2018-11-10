#'
#'       summary.kppm.R
#'
#'   $Revision: 1.6 $  $Date: 2018/11/10 09:43:05 $
#' 

summary.kppm <- function(object, ..., quick=FALSE) {
  nama <- names(object)
  result <- unclass(object)[!(nama %in% c("X", "po", "call", "callframe"))]
  ## handle old format
  if(is.null(result$isPCP)) result$isPCP <- TRUE
  ## summarise trend component
  result$trend <- summary(as.ppm(object), ..., quick=quick)
  if(identical(quick, FALSE)) {
    theta <- coef(object)
    if(length(theta) > 0) {
      vc <- vcov(object, matrix.action="warn")
      if(!is.null(vc)) {
        se <- if(is.matrix(vc)) sqrt(diag(vc)) else
        if(length(vc) == 1) sqrt(vc) else NULL
      }
      if(!is.null(se)) {
        two <- qnorm(0.975)
        lo <- theta - two * se
        hi <- theta + two * se
        zval <- theta/se
        pval <- 2 * pnorm(abs(zval), lower.tail=FALSE)
        psig <- cut(pval, c(0,0.001, 0.01, 0.05, 1),
                    labels=c("***", "**", "*", "  "),
                    include.lowest=TRUE)
        ## table of coefficient estimates with SE and 95% CI
        result$coefs.SE.CI <- data.frame(Estimate=theta, S.E.=se,
                                         CI95.lo=lo, CI95.hi=hi,
                                         Ztest=psig,
                                         Zval=zval)
      }
    }
  }
  class(result) <- "summary.kppm"
  return(result)
}

coef.summary.kppm <- function(object, ...) {
  return(object$coefs.SE.CI)
}

print.summary.kppm <- function(x, ...) {
  terselevel <- spatstat.options('terse')
  digits <- getOption('digits')
  isPCP <- x$isPCP
  splat(if(x$stationary) "Stationary" else "Inhomogeneous",
        if(isPCP) "cluster" else "Cox",
        "point process model")

  if(waxlyrical('extras', terselevel) && nchar(x$Xname) < 20)
    splat("Fitted to point pattern dataset", sQuote(x$Xname))

  Fit <- x$Fit
  
  if(waxlyrical('gory', terselevel)) {
    switch(Fit$method,
           mincon = {
             splat("Fitted by minimum contrast")
             splat("\tSummary statistic:", Fit$StatName)
             print(Fit$mcfit)
           },
           clik  =,
           clik2 = {
             splat("Fitted by maximum second order composite likelihood")
             splat("\trmax =", Fit$rmax)
             if(!is.null(wtf <- Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
             printStatus(optimStatus(Fit$clfit))
           },
           palm = {
             splat("Fitted by maximum Palm likelihood")
             splat("\trmax =", Fit$rmax)
             if(!is.null(wtf <- Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
             printStatus(optimStatus(Fit$clfit))
           },
           warning(paste("Unrecognised fitting method", sQuote(Fit$method)))
           )
  }

  # ............... trend .........................

  parbreak()
  splat("----------- TREND MODEL -----")
  print(x$trend, ...)

  # ..................... clusters ................

  tableentry <- spatstatClusterModelInfo(x$clusters)
  
  parbreak()
  splat("-----------", 
        if(isPCP) "CLUSTER" else "COX",
        "MODEL",
        "-----------")
  splat("Model:", tableentry$printmodelname(x))
  parbreak()
  
  cm <- x$covmodel
  if(!isPCP) {
    # Covariance model - LGCP only
    splat("\tCovariance model:", cm$model)
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      splat("\tCovariance parameters:",
            paste(tagvalue, collapse=", "))
    }
  }
  pa <- x$clustpar
  if (!is.null(pa)) {
    splat("Fitted",
          if(isPCP) "cluster" else "covariance",
          "parameters:")
    print(pa, digits=digits)
  }

  if(!is.null(mu <- x$mu)) {
    if(isPCP) {
      splat("Mean cluster size: ",
            if(!is.im(mu)) paste(signif(mu, digits), "points") else "[pixel image]")
    } else {
      splat("Fitted mean of log of random intensity:",
            if(!is.im(mu)) signif(mu, digits) else "[pixel image]")
    }
  }
  # table of coefficient estimates with SE and 95% CI
  if(!is.null(cose <- x$coefs.SE.CI)) {
    parbreak()
    splat("Final standard error and CI")
    splat("(allowing for correlation of",
          if(isPCP) "cluster" else "Cox",
          "process):")
    print(cose)
  }
  invisible(NULL)
}
