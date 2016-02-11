#
# simulation of FITTED model
#
#  $Revision: 1.33 $ $Date: 2016/02/11 10:17:12 $
#
#
rmh.ppm <- function(model, start = NULL,
                    control = default.rmhcontrol(model, w=w),
                    ...,
                    w = NULL, project=TRUE, verbose=TRUE,
                    new.coef=NULL) {
  verifyclass(model, "ppm")
  argh <- list(...)

  if(is.null(control)) {
    control <- default.rmhcontrol(model, w=w)
  } else {
    control <- rmhcontrol(control)
  }

  # override 
  if(length(list(...)) > 0)
    control <- update(control, ...)
  
  # convert fitted model object to list of parameters for rmh.default
  X <- rmhmodel(model, w=w, verbose=verbose, project=project, control=control,
                new.coef=new.coef)

  # set initial state

  if(is.null(start)) {
    datapattern <- data.ppm(model)
    start <- rmhstart(n.start=datapattern$n)
  }
  
  # call rmh.default 
  # passing only arguments unrecognised by rmhcontrol
  known <- names(argh) %in% names(formals(rmhcontrol.default))
  fargs <- argh[!known]

  Y <- do.call(rmh.default,
               append(list(model=X, start=start, control=control,
                           verbose=verbose),
                      fargs))
  return(Y)
}

simulate.ppm <- function(object, nsim=1, ...,
                         singlerun=FALSE,
                         start = NULL,
                         control = default.rmhcontrol(object, w=w),
                         w = NULL,
                         project=TRUE,
                         new.coef=NULL,
                         verbose=FALSE,
                         progress=(nsim > 1),
                         drop=FALSE) {
  verifyclass(object, "ppm")
  argh <- list(...)
  if(nsim == 0) return(list())

  starttime = proc.time()
  
  # set up control parameters
  if(missing(control) || is.null(control)) {
    rcontr <- default.rmhcontrol(object, w=w)
  } else {
    rcontr <- rmhcontrol(control)
  }
  if(singlerun) {
    # allow nsave, nburn to determine nrep
    nsave <- resolve.1.default("nsave", list(...), as.list(rcontr),
                               .MatchNull=FALSE)
    nburn <- resolve.1.default("nburn", list(...), as.list(rcontr),
                               list(nburn=nsave),
                               .MatchNull=FALSE)
    if(!is.null(nsave)) {
      nrep <- nburn + (nsim-1) * nsave
      rcontr <- update(rcontr, nrep=nrep, nsave=nsave, nburn=nburn)
    } 
  }
  # other overrides
  if(length(list(...)) > 0)
    rcontr <- update(rcontr, ...)

  # Set up model parameters for rmh
  rmodel <- rmhmodel(object, w=w, verbose=FALSE, project=TRUE, control=rcontr,
                     new.coef=new.coef)
  if(is.null(start)) {
    datapattern <- data.ppm(object)
    start <- rmhstart(n.start=datapattern$n)
  }
  rstart <- rmhstart(start)

  #########
  
  if(singlerun && nsim > 1) {
    # //////////////////////////////////////////////////
    # execute one long run and save every k-th iteration
    if(is.null(rcontr$nsave)) {
      # determine spacing between subsamples
      if(!is.null(rcontr$nburn)) {
        nsave <- max(1, with(rcontr, floor((nrep - nburn)/(nsim-1))))
      } else {
        # assume nburn = 2 * nsave
        nsave <- max(1, with(rcontr, floor(nrep/(nsim+1))))
        nburn <- 2 * nsave
      }
      rcontr <- update(rcontr, nsave=nsave, nburn=nburn)
    }
    # check nrep is enough
    nrepmin <- with(rcontr, nburn + (nsim-1) * nsave)
    if(rcontr$nrep < nrepmin)
      rcontr <- update(rcontr, nrep=nrepmin)
    # OK, run it
    if(progress) {
      cat(paste("Generating", nsim, "simulated patterns in a single run ... ")) 
      flush.console()
    }
    Y <- rmh(rmodel, rstart, rcontr, verbose=verbose)
    if(progress)
      cat("Done.\n")
    # extract sampled states
    out <- attr(Y, "saved")
    if(length(out) != nsim)
      stop(paste("Internal error: wrong number of simulations generated:",
                 length(out), "!=", nsim))
  } else {
    # //////////////////////////////////////////////////
    # execute 'nsim' independent runs
    out <- list()
    # pre-digest arguments
    rmhinfolist <- rmh(rmodel, rstart, rcontr, preponly=TRUE, verbose=verbose)
    # go
    if(nsim > 0) {
      if(progress) {
        cat(paste("Generating", nsim, "simulated", 
                  ngettext(nsim, "pattern", "patterns"),
                  "..."))
        flush.console()
      }
      # call rmh
      # passing only arguments unrecognised by rmhcontrol
      known <- names(argh) %in% names(formals(rmhcontrol.default))
      fargs <- argh[!known]
      rmhargs <- append(list(InfoList=rmhinfolist, verbose=verbose), fargs)
      if(progress)
        pstate <- list()
      for(i in 1:nsim) {
        out[[i]] <- do.call(rmhEngine, rmhargs)
        if(progress) pstate <- progressreport(i, nsim, state=pstate)
      }
    }
  }
  if(nsim == 1 && drop) {
    out <- out[[1]]
  } else {
    out <- as.solist(out)
    if(nsim > 0)
      names(out) <- paste("Simulation", 1:nsim)
  }
  out <- timed(out, starttime=starttime)
  return(out)
}  
