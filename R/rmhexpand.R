#
#   rmhexpand.R
#
#   Rules/data for expanding the simulation window in rmh
#
#   $Revision: 1.6 $  $Date: 2012/05/11 12:38:42 $
#

# Establish names and rules for each type of expansion

RmhExpandRule <- local({

  .RmhExpandTable <-
  list(area=list(descrip ="Area expansion factor",
                 minval = 1,
                 expands = function(x) { unname(x) > 1 }),
       length=list(descrip ="Length expansion factor",
                   minval = 1,
                   expands = function(x) { unname(x) > 1 }),
       distance=list(descrip="Expansion buffer distance",
                     minval = 0,
                     expands = function(x) { unname(x) > 0 }))
  
  RmhExpandRule <- function(nama) {
    if(length(nama) == 0) nama <- "area"
    if(length(nama) > 1)
      stop("Internal error: too many names in RmhExpandRule", call.=FALSE)
    if(!(nama %in% names(.RmhExpandTable)))
      stop(paste("Internal error: unrecognised expansion type",
                 sQuote(nama)),
           call.=FALSE)
    return(.RmhExpandTable[[nama]])
  }
  RmhExpandRule
})
  
                        
rmhexpand <- function(x=NULL, ..., area=NULL, length=NULL, distance=NULL) {
  trap.extra.arguments(..., .Context="In rmhexpand")
  # check for incompatibility
  n <- (!is.null(x)) + (!is.null(area)) +
       (!is.null(length)) + (!is.null(distance))
  if(n > 1) stop("Only one argument should be given")
  # absorb other arguments into 'x'
  if(is.null(x) && n > 0) {
      if(!is.null(area)) x <- c(area=area)
      if(!is.null(length)) x <- c(length=length)
      if(!is.null(distance)) x <- c(distance=distance)
  }
  if(is.null(x)) {
    # No expansion rule supplied.
    # Use spatstat default, indicating that the user did not choose it.
    force.exp <- force.noexp <- FALSE
    x <- spatstat.options("expand")
    x <- rmhexpand(x)$expand
  } else {
    # process x
    if(inherits(x, "rmhexpand"))
      return(x)
    if(is.owin(x)) {
      force.exp <- TRUE
      force.noexp <- FALSE
    } else {
      # expecting c(name=value) or list(name=value)
      if(is.list(x))
        x <- unlist(x)
      if(!is.numeric(x))
        stop(paste("Expansion argument must be either",
                   "a number, a window, or NULL.\n"))
      # x is numeric
      check.1.real(x, "In rmhexpand(x)")
      explain.ifnot(is.finite(x), "In rmhexpand(x)")
      # an unlabelled numeric value is interpreted as an area expansion factor
      if(!any(nzchar(names(x))))
        names(x) <- "area"
      # validate
      rule <- RmhExpandRule(names(x))
      if(x < rule$minval) {
        warning(paste(rule$descrip, "<", rule$minval,
                      "has been reset to", rule$minval),
                call.=FALSE)
        x[] <- rule$minval
      }
      force.exp <- rule$expands(x)
      force.noexp <- !force.exp
    }
  }
  result <- list(expand=x, force.exp=force.exp, force.noexp=force.noexp)
  class(result) <- "rmhexpand"
  return(result)
}

.no.expansion <- list(expand=c(area=1), force.exp=FALSE, force.noexp=TRUE)
class(.no.expansion) <- "rmhexpand"

print.rmhexpand <- function(x, ..., prefix=TRUE) {
  if(prefix) cat("Expand the simulation window? ")
  if(x$force.noexp) {
    cat("No.\n")
  } else {
    if(x$force.exp) cat("Yes:\n") else cat("Not determined. Default is:\n")

    y <- x$expand
    if(is.null(y)) {
      print(rmhexpand(spatstat.options("expand")), prefix=FALSE)
    } else if(is.numeric(y)) {
      descrip <- RmhExpandRule(names(y))$descrip
      cat(paste("\t", descrip, unname(y), "\n"))
    } else {
      print(y)
    }
  }
  return(invisible(NULL))
}

summary.rmhexpand <- function(object, ...) {
  decided <- with(object, force.exp || force.noexp)
  ex <- object$expand
  if(is.null(ex))
    ex <- rmhexpand(spatstat.options("expand"))$expand
  if(is.owin(ex)) {
    willexpand <- TRUE
    descrip <- "Window"
  } else if(is.numeric(ex)) {
    rule <- RmhExpandRule(names(ex))
    descrip    <- rule$descrip
    willexpand <- if(object$force.exp) TRUE else
                  if(object$force.noexp) FALSE else
                  (unname(ex) > rule$minval)
  } else stop("Internal error: unrecognised format in summary.rmhexpand",
              call.=FALSE)
              
  out <- list(rule.decided=decided,
              window.decided=decided && is.owin(ex), 
              expand=ex,
              descrip=descrip,
              willexpand=willexpand)
  class(out) <- "summary.rmhexpand"
  return(out)
}

print.summary.rmhexpand <- function(x, ...) {
  cat("Expansion rule\n")
  ex <- x$expand
  if(x$window.decided) {
    cat("Window is decided.\n")
    print(ex)
  } else {
    if(x$rule.decided) {
      cat("Rule is decided.\n")
    } else {
      cat("Rule is not decided.\nDefault is:\n")
    }
    if(!x$willexpand) {
      cat("No expansion\n")
    } else {
      if(is.numeric(ex)) cat(paste(x$descrip, ex, "\n")) else print(ex)
    }
  }
  return(invisible(NULL))
}

expand.owin <- function(W, ...) {
  ex <- list(...)
  if(length(ex) > 1) stop("Too many arguments")
  # get an rmhexpand object
  if(inherits(ex[[1]], "rmhexpand")) {
    ex <- ex[[1]]
  } else ex <- do.call("rmhexpand", ex)
  f <- ex$expand
  if(is.null(f)) return(W)
  if(is.owin(f)) return(f)
  if(!is.numeric(f)) stop("Format not understood")
  switch(names(f),
         area = {
           if(f == 1)
             return(W)
           bb <- bounding.box(W)
           xr <- bb$xrange
           yr <- bb$yrange
           fff <- (sqrt(f) - 1)/2
           Wexp <- grow.rectangle(bb, fff * diff(xr), fff * diff(yr))
         },
         length = {
           if(f == 1)
             return(W)
           bb <- bounding.box(W)
           xr <- bb$xrange
           yr <- bb$yrange
           fff <- (f - 1)/2
           Wexp <- grow.rectangle(bb, fff * diff(xr), fff * diff(yr))
         },
         distance = {
           if(f == 0)
             return(W)
           Wexp <- if(is.rectangle(W)) grow.rectangle(W, f) else dilation(W, f)
         },
         stop("Internal error: unrecognised type")
         )
  return(Wexp)
}

will.expand <- function(x) {
  stopifnot(inherits(x, "rmhexpand"))
  if(x$force.exp) return(TRUE)
  if(x$force.noexp) return(FALSE)
  return(summary(x)$willexpand)
}

is.expandable <- function(x) { UseMethod("is.expandable") }

change.default.expand <- function(x, newdefault) {
  stopifnot(inherits(x, "rmhexpand"))
  decided <- with(x, force.exp || force.noexp)
  if(!decided)
    x$expand <- rmhexpand(newdefault)$expand
  return(x)
}

  
