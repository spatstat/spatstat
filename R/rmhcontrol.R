#
#
#   rmhcontrol.R
#
#   $Revision: 1.25 $  $Date: 2014/10/24 00:22:30 $
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
    do.call("rmhcontrol.default", argz[[1]])
  else
    do.call.matched("rmhcontrol.default", argz)
}

rmhcontrol.default <- function(..., p=0.9, q=0.5, nrep=5e5,
                        expand=NULL, periodic=NULL, ptypes=NULL,
                        x.cond=NULL, fixall=FALSE, nverb=0,
                        nsave=NULL, nburn=nsave, track=FALSE)
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
    if(!is.numeric(nsave) || length(nsave) != 1
       || nsave < 0 || nsave >= nrep)
      stop("nsave should be an integer < nrep")
    if(is.null(nburn)) nburn <- min(nsave, nrep-nsave)
    if(!is.null(nburn)) stopifnot(nburn + nsave <= nrep)
  }
  stopifnot(is.logical(track))

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
                  "is equivalent to unconditional simulation"))
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
	warning("fixall = TRUE conflicts with p < 1. Ignored.\n")

###############################################################  
# `expand' determines expansion of the simulation window

  expand <- rmhexpand(expand)

# No expansion is permitted if we are conditioning on the
# number of points
  
  if(fixing != "none") {
    if(expand$force.exp)
      stop(paste("When conditioning on the number of points,",
                 "no expansion may be done.\n"), call.=FALSE)
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
              track=track)
  class(out) <- c("rmhcontrol", class(out))
  return(out)
}

print.rmhcontrol <- function(x, ...) {
  verifyclass(x, "rmhcontrol")

  cat("Metropolis-Hastings algorithm control parameters\n")
  cat(paste("Probability of shift proposal: p =", x$p, "\n"))
  if(x$fixing == "none") {
    cat(paste("Conditional probability of death proposal: q =", x$q, "\n"))
    if(!is.null(x$ptypes)) {
      cat("Birth proposal probabilities for each type of point:\n")
      print(x$ptypes)
    }
  }
  switch(x$fixing,
         none={},
         n.total=cat("The total number of points is fixed\n"),
         n.each.type=cat("The number of points of each type is fixed\n"))
  switch(x$condtype,
         none={},
         window={
           cat(paste("Conditional simulation given the",
                     "configuration in a subwindow\n"))
           print(x$x.cond$window)
         },
         Palm={
           cat("Conditional simulation of Palm type\n")
         })
  cat(paste("Number of M-H iterations: nrep =", x$nrep, "\n"))
  if(x$saving) 
    cat(paste("Save point pattern every", x$nsave, "iterations",
              "after a burn-in of", x$nburn, "iterations\n"))
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

default.rmhcontrol <- function(model) {
  # set default for 'expand'
  return(rmhcontrol(expand=default.expand(model)))
}

update.rmhcontrol <- function(object, ...) {
  do.call.matched("rmhcontrol.default",
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

