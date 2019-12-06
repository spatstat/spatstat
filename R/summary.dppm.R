#'
#'  summary.dppm.R
#'
#'  $Revision: 1.4 $ $Date: 2019/12/06 01:35:46 $

summary.dppm <- function(object, ..., quick=FALSE) {
  nama <- names(object)
  result <- unclass(object)[!(nama %in% c("X", "po", "call", "callframe"))]
  ## Fitting information
  result$has.subset <- "subset" %in% names(object$call)
  ## Summarise trend component
  result$trend <- summary(as.ppm(object), ..., quick=quick)
  ## repulsion strength
  result$repul <- mean(repul(object))
  #' pack up
  class(result) <- "summary.dppm"
  return(result)
}

print.summary.dppm <- function(x, ...) {
  terselevel <- spatstat.options('terse')
  digits <- getOption('digits')
  splat(if(x$stationary) "Stationary" else "Inhomogeneous",
        "determinantal point process model")

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

  # ..................... determinantal part  ................

  parbreak()
  splat("---------- DETERMINANTAL STRUCTURE -----------------")
  print(x$fitted)

  parbreak()
  splat(if(x$stationary) "Strength" else "(Average) strength",
        "of repulsion:", signif(x$repul, digits))
  
  return(invisible(NULL))
}
