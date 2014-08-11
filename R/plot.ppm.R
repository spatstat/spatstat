#
#    plot.ppm.S
#
#    $Revision: 2.8 $    $Date: 2008/08/12 08:26:45 $
#
#    plot.ppm()
#         Plot a point process model fitted by ppm().
#        
#
#
plot.ppm <- function(x, ngrid = c(40,40),
		     superimpose = TRUE,
                     trend=TRUE, cif=TRUE, se=TRUE, 
                     pause = interactive(),
                     how=c("persp","image", "contour"),
                     plot.it=TRUE,
                     locations=NULL, covariates=NULL, ...)
{
  model <- x
#       Plot a point process model fitted by ppm().
#
  verifyclass(model, "ppm")
#
#       find out what kind of model it is
#
  mod <- summary(model)
  stationary <- mod$stationary
  poisson    <- mod$poisson
  marked     <- mod$marked
  multitype  <- mod$multitype
  data       <- mod$entries$data
        
  if(marked) {
    if(!multitype)
      stop("Not implemented for general marked point processes")
    else
      mrkvals <- levels(marks(data))
  } else mrkvals <- 1
  ntypes <- length(mrkvals)
        
#
#        Interpret options
#        -----------------
#        
#        Whether to plot trend, cif, se
        
  if(!trend && !cif && !se) {
    cat(paste("Nothing plotted;", sQuote("trend"), ",", sQuote("cif"),
              "and", sQuote("se"), "are all FALSE\n"))
    return(invisible(NULL))
  }
#        Suppress uninteresting plots
#        unless explicitly instructed otherwise
  if(missing(trend))
    trend <- !stationary
  if(missing(cif))
    cif <- !poisson
  if(missing(se))
    se <- poisson && !stationary 
  else if(se && !poisson) {
      warning(paste("standard error calculation",
                  "is only implemented for Poisson models"))
      se <- FALSE
  }
  if(!trend && !cif && !se) {
    cat("Nothing plotted -- all plots selected are flat surfaces.\n")
    return(invisible(NULL))
  }
#
#  style of plot: suppress pseudo-default
#  
    if(missing(how))
      how <- "image"
#
#
#        Do the prediction
#        ------------------

  out <- list()
  surftypes <- c("trend","cif","se")[c(trend,cif,se)]
  ng <- if(missing(ngrid) && !missing(locations)) NULL else ngrid

  for (ttt in surftypes) {
    p <- predict(model,
                   ngrid=ng, locations=locations, covariates=covariates,
                   type = ttt)
    if(is.im(p))
      p <- list(p)
    out[[ttt]] <- p
  }

#        Make it a plotppm object
#        ------------------------  
  
  class(out) <- "plotppm"
  attr(out, "mrkvals") <- mrkvals

#        Actually plot it if required
#        ----------------------------  
  if(plot.it) {
    if(!superimpose)
      data <- NULL
    plot(out,data=data,trend=trend,cif=cif,se=se,how=how,pause=pause, ...)
  }

  
  return(invisible(out)) 
}

