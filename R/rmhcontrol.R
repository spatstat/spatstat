#
#
#   rmhcontrol.R
#
#   $Revision: 1.35 $  $Date: 2019/12/31 04:56:58 $
#
#

rmhcontrol <- function(...) {
  UseMethod("rmhcontrol")
}

rmhcontrol.rmhcontrol <- function(...) {
  argz <- list(...)
  if(length(argz) == 1)
    return(argz[[1]])
  stop("Arguments not understood")
}

rmhcontrol.list <- function(...) {
  argz <- list(...)
  nama <- names(argz)
  if(length(argz) == 1 && !any(nzchar(nama)))
    do.call(rmhcontrol.default, argz[[1]])
  else
    do.call.matched(rmhcontrol.default, argz)
}

rmhcontrol.default <- function(..., p=0.9, q=0.5, nrep=5e5,
                        expand=NULL, periodic=NULL, ptypes=NULL,
                        x.cond=NULL, fixall=FALSE, nverb=0,
                        nsave=NULL, nburn=nsave, track=FALSE,
                        pstage=c("block", "start"))
{
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    # allow rmhcontrol(NULL), otherwise flag an error
    if(!(nargh == 1 && is.null(argh[[1]])))
      stop(paste("Unrecognised arguments to rmhcontrol;",
                 "valid arguments are listed in help(rmhcontrol.default)"))
  }
  # impose default values
  if(missing(p)) p <- spatstat.options("rmh.p")
  if(missing(q)) q <- spatstat.options("rmh.q")
  if(missing(nrep)) nrep <- spatstat.options("rmh.nrep")
  # validate arguments
  if(!is.numeric(p) || length(p) != 1
     || p < 0 || p > 1)
    stop("p should be a number in [0,1]")
  if(!is.numeric(q) || length(q) != 1
     || q < 0 || q > 1)
    stop("q should be a number in [0,1]")
  if(!is.numeric(nrep) || length(nrep) != 1
     || nrep < 1)
    stop("nrep should be an integer >= 1")
  nrep <- as.integer(nrep)
  if(!is.numeric(nverb) || length(nverb) != 1
     || nverb < 0 || nverb > nrep)
    stop("nverb should be an integer <= nrep")
  nverb <- as.integer(nverb)
  if(!is.logical(fixall) || length(fixall) != 1)
    stop("fixall should be a logical value")
  if(!is.null(periodic) && (!is.logical(periodic) || length(periodic) != 1))
    stop(paste(sQuote("periodic"), "should be a logical value or NULL"))
  if(saving <- !is.null(nsave)) {
    nsave <- as.integer(as.vector(nsave))
    if(length(nsave) == 1L) {
      if(nsave <= 0)
        stop("nsave should be a positive integer")
      stopifnot(nsave < nrep)
    } else {
      stopifnot(all(nsave > 0))
      stopifnot(sum(nsave) <= nrep)
    }
    if(missing(nburn) || is.null(nburn)) {
      nburn <- min(nsave[1], nrep-sum(nsave))
    } else {
      check.1.integer(nburn)
      stopifnot(nburn + sum(nsave) <= nrep)
    }
  }
  stopifnot(is.logical(track))
  pstage <- match.arg(pstage) 

#################################################################
# Conditioning on point configuration
#
# condtype = "none": no conditioning
# condtype = "Palm": conditioning on the presence of specified points
# condtype = "window": conditioning on the configuration in a subwindow
#
  if(is.null(x.cond)) {
    condtype <- "none"
    n.cond <- NULL
  } else if(is.ppp(x.cond)) {
    condtype <- "window"
    n.cond <- x.cond$n
  } else if(is.data.frame(x.cond)) {
    if(ncol(x.cond) %in% c(2,3)) {
      condtype <- "Palm"
      n.cond <- nrow(x.cond)
    } else stop("Wrong number of columns in data frame x.cond")
  } else if(is.list(x.cond)) {
    if(length(x.cond) %in% c(2,3)) {
      x.cond <- as.data.frame(x.cond)
      condtype <- "Palm"
      n.cond <- nrow(x.cond)
    } else stop("Wrong number of components in list x.cond")
  } else stop("Unrecognised format for x.cond")

  if(condtype == "Palm" && n.cond == 0) {
    warning(paste("Ignored empty configuration x.cond;",
                  "conditional (Palm) simulation given an empty point pattern",
                  "is equivalent to unconditional simulation"), call.=FALSE)
    condtype <- "none"
    x.cond <- NULL
    n.cond <- NULL
  }
    
#################################################################
# Fixing the number of points?
#  
# fixcode = 1 <--> no conditioning
# fixcode = 2 <--> conditioning on n = number of points
# fixcode = 3 <--> conditioning on the number of points of each type.

  fixcode    <- 2 - (p<1) + fixall - fixall*(p<1)
  fixing <- switch(fixcode, "none", "n.total", "n.each.type")
  
# Warn about silly combination
  if(fixall && p < 1)
	warning("fixall = TRUE conflicts with p < 1. Ignored.", call.=FALSE)

###############################################################  
# `expand' determines expansion of the simulation window

  expand <- rmhexpand(expand)

# No expansion is permitted if we are conditioning on the
# number of points
  
  if(fixing != "none") {
    if(expand$force.exp)
      stop(paste("When conditioning on the number of points,",
                 "no expansion may be done."), call.=FALSE)
    # no expansion
    expand <- .no.expansion
  }

###################################################################
# return augmented list  
  out <- list(p=p, q=q, 
              nrep=nrep, nverb=nverb,
              expand=expand, 
              periodic=periodic, 
              ptypes=ptypes,
              fixall=fixall,
              fixcode=fixcode,
              fixing=fixing,
              condtype=condtype,
              x.cond=x.cond,
              saving=saving, nsave=nsave, nburn=nburn,
              track=track, pstage=pstage)
  class(out) <- c("rmhcontrol", class(out))
  return(out)
}

print.rmhcontrol <- function(x, ...) {
  verifyclass(x, "rmhcontrol")

  splat("Metropolis-Hastings algorithm control parameters")
  splat("Probability of shift proposal: p =", x$p)
  if(x$fixing == "none") {
    splat("Conditional probability of death proposal: q =", x$q)
    if(!is.null(x$ptypes)) {
      splat("Birth proposal probabilities for each type of point:")
      print(x$ptypes)
    }
  }
  switch(x$fixing,
         none={},
         n.total=splat("The total number of points is fixed"),
         n.each.type=splat("The number of points of each type is fixed"))
  switch(x$condtype,
         none={},
         window={
           splat("Conditional simulation given the",
                 "configuration in a subwindow")
           print(x$x.cond$window)
         },
         Palm={
           splat("Conditional simulation of Palm type")
         })
  splat("Number of M-H iterations: nrep =", x$nrep)
  if(x$saving) {
    nsave <- x$nsave
    len <- length(nsave)
    howmany <- if(len == 1L) nsave else
               if(len < 5L) commasep(nsave) else
               paste(paste(nsave[1:5], collapse=", "), "[...]")
    splat("After a burn-in of", x$nburn, "iterations,",
          "save point pattern after every", howmany, "iterations.")
  }
  pstage <- x$pstage %orifnull% "start"
  hdr <- "Generate random proposal points:"
  switch(pstage,
         start = splat(hdr, "at start of simulations."),
         block = splat(hdr,
                       "before each block of",
                       if(length(x$nsave) == 1L) x$nsave else "",
                       "iterations."))
  cat(paste("Track proposal type and acceptance/rejection?",
            if(x$track) "yes" else "no", "\n"))
  if(x$nverb > 0)
    cat(paste("Progress report every nverb=", x$nverb, "iterations\n"))
  else
    cat("No progress reports (nverb = 0).\n")

  # invoke print.rmhexpand
  print(x$expand)

  cat("Periodic edge correction? ")
  if(is.null(x$periodic)) cat("Not yet determined.\n") else 
  if(x$periodic) cat("Yes.\n") else cat("No.\n")
  #
  
  return(invisible(NULL))
}

default.rmhcontrol <- function(model, w=NULL) {
  # set default for 'expand'
  return(rmhcontrol(expand=default.expand(model, w=w)))
}

update.rmhcontrol <- function(object, ...) {
  do.call.matched(rmhcontrol.default,
                  resolve.defaults(list(...), as.list(object),
                                   .StripNull=TRUE))
}

rmhResolveControl <- function(control, model) {
  # adjust control information once the model is known
  stopifnot(inherits(control, "rmhcontrol"))
  # change *default* expansion rule to something appropriate for model
  # (applies only if expansion rule is undecided)
  control$expand <- change.default.expand(control$expand, default.expand(model))
  return(control)
}

