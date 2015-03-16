#
#    QQ plot of smoothed residual field against model
#
#  qqplot.ppm()       QQ plot (including simulation)
#
#  $Revision: 1.25 $   $Date: 2015/03/16 11:02:36 $
#

qqplot.ppm <-
  function(fit, nsim=100, expr=NULL, ..., type="raw", style="mean",
           fast=TRUE, verbose=TRUE, plot.it=TRUE,
           dimyx=NULL, nrep=if(fast) 5e4 else 1e5,
           control=update(default.rmhcontrol(fit), nrep=nrep),
           saveall=FALSE,
           monochrome=FALSE,
           limcol=if(monochrome) "black" else "red",
           maxerr=max(100, ceiling(nsim/10)),
           check=TRUE, repair=TRUE)
{
  verifyclass(fit, "ppm")

  if(check && damaged.ppm(fit)) {
    if(!repair)
      stop("object format corrupted; try update(fit, use.internal=TRUE)")
    message("object format corrupted; repairing it.")
    fit <- update(fit, use.internal=TRUE)
  }
  
  if(fast) {
    oldnpixel <- spatstat.options("npixel")
    if(is.null(dimyx)) 
      dimyx <- pmin(40, rev(oldnpixel))
    spatstat.options(npixel=rev(dimyx))
  } 
    
  ################   How to evaluate residuals ##########################
  
  # Quantiles of the residual field will be computed.
  residualfield <- function(fit, ...) {
    d <- diagnose.ppm(fit, which="smooth",
                      plot.it=FALSE, compute.cts=FALSE, compute.sd=FALSE,
                      check=FALSE, ...)
    return(d$smooth$Z$v)
  }

  # Data values
  dat <- residualfield(fit, type=type, ..., dimyx=dimyx)

  # How to refit the model properly!
  refit <- function(fit, pattern) {
    update.ppm(fit, Q=pattern, use.internal=TRUE)
  }
  
  ##################  How to perform simulations?  #######################

  simulate.from.fit <- is.null(expr)

  how.simulating <- if(simulate.from.fit)
    "simulating from fitted model" else paste("evaluating", sQuote("expr"))  

  if(!simulate.from.fit) {
     # 'expr' will be evaluated 'nsim' times
    if(!is.expression(expr))
      stop(paste("Argument", sQuote("expr"), "should be an expression"))
  } else{
    # We will simulate from the fitted model 'nsim' times
    # and refit the model to these simulations

    # prepare rmh arguments
    rcontrol <- rmhcontrol(control)
    rmodel   <- rmhmodel(fit, control=rcontrol, project=FALSE, verbose=verbose)
    rstart   <- rmhstart(n.start=data.ppm(fit)$n)
    # pre-digest arguments
    rmhinfolist <- rmh(rmodel, rstart, rcontrol, preponly=TRUE, verbose=FALSE)
    
    # expression to be evaluated each time
    expr <- expression(
        refit(fit, 
              rmhEngine(rmhinfolist, verbose=FALSE)))
  }

  ######  Perform simulations
  if(verbose) cat(paste("Simulating", nsim, "realisations... "))
  simul.sizes <- numeric(nsim)
  i <- 0
  ierr <- 0
  repeat {
    # protect from randomly-generated crashes in gam
    ei <- try(eval(expr),silent=!verbose)
    if(inherits(ei, "try-error")) {
      # error encountered in evaluating 'expr'
      ierr <- ierr + 1
      if(ierr > maxerr) 
        stop(paste("Exceeded maximum of", maxerr,
                   "failures in", how.simulating,
                   "after generating only", i, "realisations"))
      else break
    } else {
      # simulation successful
      i <- i + 1
      fiti <-
      if(simulate.from.fit)
        ei
      else if(is.ppm(ei))
        ei
      else if(is.ppp(ei))
        refit(fit, ei)
      else
        stop("result of eval(expr) is not a ppm or ppp object")
      # diagnostic info
      simul.sizes[i] <- data.ppm(fiti)$n
      # compute residual field
      resi <- residualfield(fiti, type=type, ..., dimyx=dimyx)
      if(i == 1)
        sim <- array(, dim=c(dim(resi), nsim))
      sim[,,i] <- resi
      if(verbose) 
        progressreport(i, nsim)
      if(i >= nsim)
        break
    }
  }

  ###### Report diagnostics
  if(ierr > 0)
    cat(paste("\n\n**Alert:",
              ierr, "failures occurred in", how.simulating, "\n\n"))
  nempty <- sum(simul.sizes == 0)
  if(nempty > 0)
    cat(paste("\n\n**Alert:",
              nempty, "out of", nsim,
              "simulated patterns were empty.\n\n"))
  else
    cat(paste("\nDiagnostic info:\n",
              "simulated patterns contained an average of",
              mean(simul.sizes), "points.\n"))
  if(nempty == nsim)
    warning("All simulated patterns were empty")
  ############ Plot them
  switch(style,
         classical = {
           rr <- range(c(dat,sim))
           result <- qqplot(sim, dat, xlim=rr, ylim=rr, asp=1.0,
                            xlab="Quantiles of simulation",
                            ylab="Quantiles of data",plot.it=plot.it)
           title(sub=paste("Residuals:", type))
           abline(0,1, lty=2)
           result <- append(result,
                            list(data=dat,
                                 sim=sim,
                                 xlim=rr,
                                 ylim=rr,
                                 xlab="Quantiles of simulation",
                                 ylab="Quantiles of data",
                                 rtype=type,
                                 nsim=nsim,
                                 fit=fit,
                                 expr=
                                 if(simulate.from.fit) NULL else deparse(expr),
                                 simulate.from.fit=simulate.from.fit
                                 )
                            )
         },
         mean = {
           # compute quantiles corresponding to probabilities p[i]
           # separately in each realisation.
           if(verbose) cat("Calculating quantiles...")
           if(fast) {
             p <- ppoints(min(100,length(dat)), 3/8)
             qsim <- apply(sim, 3, quantile, probs=p, na.rm=TRUE)
           } else {
             qsim <- apply(sim, 3, sort, na.last=TRUE)
           }
           if(verbose) cat("averaging...")
           # sample mean of each quantile
           meanq <- apply(qsim, 1, mean, na.rm=TRUE)
           # et cetera
           varq <- apply(qsim, 1, var, na.rm=TRUE)
           sdq <- sqrt(varq)
           q.025 <- apply(qsim, 1, quantile, probs=0.025, na.rm=TRUE)
           q.975 <- apply(qsim, 1, quantile, probs=0.975, na.rm=TRUE)
  
           rr <- range(c(meanq,dat), na.rm=TRUE)

           dats <- if(fast) quantile(dat, probs=p, na.rm=TRUE) else sort(dat, na.last=TRUE)

           if(verbose) cat("..Done.\n")
           if(plot.it) {
           	plot(meanq, dats,
                     xlab="Mean quantile of simulations", ylab="data quantile",
                     xlim=rr, ylim=rr, asp=1.0)
                abline(0,1)
                lines(meanq, q.025, lty=2, col=limcol)
                lines(meanq, q.975, lty=2, col=limcol)
                title(sub=paste("Residuals:", type))
           }
           result <- list(x=meanq, y=dats, sdq=sdq,
                          q.025=q.025, q.975=q.975,
                          data=dat, sim=sim,
                          xlim=rr, ylim=rr,
                          xlab="Mean quantile of simulations",
                          ylab="data quantile",
                          rtype=type,
                          nsim=nsim,
                          fit=fit,
                          expr=if(simulate.from.fit) NULL else deparse(expr),
                          simulate.from.fit=simulate.from.fit)
         },
         stop(paste("Unrecognised option for", sQuote("style")))
       )

  # Throw out baggage if not wanted         
  if(!saveall) {
    result$fit <- summary(fit, quick=TRUE)
    result$sim <- NULL
  }
         
  # reset npixel
  if(fast)
    spatstat.options(npixel=oldnpixel)
  #
  class(result) <- c("qqppm", class(result))
  return(invisible(result))
}

plot.qqppm <- function(x, ..., limits=TRUE,
                       monochrome=spatstat.options('monochrome'),
                       limcol=if(monochrome) "black" else "red") {
  stopifnot(inherits(x, "qqppm"))
  default.type <- if(length(x$x) > 150) "l" else "p"
  myplot <- function(object,
                     xlab = object$xlab, ylab = object$ylab,
                     xlim = object$xlim, ylim = object$ylim,
                     asp = 1,
                     type = default.type,
                     ..., limits=TRUE) {
    plot(object$x, object$y, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, asp = asp, type = type, ...)
    abline(0, 1)
    
    if(limits) {
      if(!is.null(object$q.025))
        lines(object$x, object$q.025, lty = 2, col=limcol)
      if(!is.null(object$q.975))
        lines(object$x, object$q.975, lty = 2, col=limcol)
    }
    title(sub=paste("Residuals:", object$rtype))
  }
  myplot(x, ..., limits=limits)
  return(invisible(x))
}

  
print.qqppm <- function(x, ...) {
  stopifnot(inherits(x, "qqppm"))
  cat(paste("Q-Q plot of point process residuals ",
            "of type", sQuote(x$rtype), "\n",
            "based on ", x$nsim, " simulations\n",
            sep=""))
  if(x$simulate.from.fit) {
    fit  <- x$fit
    sumfit <- if(is.ppm(fit)) summary(fit, quick=TRUE)
              else if(inherits(fit, "summary.ppm")) fit
              else list(name="(unrecognised format)")
    cat(paste("\nSimulations from fitted model:",
              sumfit$name, "\n"))
  } else {
    cat("Simulations obtained by evaluating the following expression:\n")
    print(x$expr)
  } 
  invisible(NULL)
}
