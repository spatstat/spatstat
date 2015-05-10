#
# Functions for extracting and setting the name of the unit of length
#
#   $Revision: 1.20 $   $Date: 2015/05/10 01:44:48 $
#
#

unitname <- function(x) {
  UseMethod("unitname")
}

unitname.owin <- function(x) {
  u <- as.units(x$units)
  return(u)
}

unitname.ppp <- function(x) {
  u <- as.units(x$window$units)
  return(u)
}

unitname.im <- function(x) {
  u <- as.units(x$units)
  return(u)
}

unitname.default <- function(x) {
  return(as.units(attr(x, "units")))
}

"unitname<-" <- function(x, value) {
  UseMethod("unitname<-")
}

"unitname<-.owin" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

"unitname<-.ppp" <- function(x, value) {
  w <- x$window
  unitname(w) <- value
  x$window <- w
  return(x)
}

"unitname<-.im" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

"unitname<-.default" <- function(x, value) {
  if(is.null(x)) return(x)
  attr(x, "units") <- as.units(value)
  return(x)
}


###  class 'units'

makeunits <- function(sing="unit", plur="units", mul = 1) {
  if(!is.character(sing))
    stop("In unit name, first entry should be a character string")
  if(!is.character(plur))
    stop("In unit name, second entry should be a character string")
  if(!is.numeric(mul)) {
    mul <- try(as.numeric(mul), silent=TRUE)
    if(inherits(mul, "try-error"))
      stop("In unit name, third entry should be a number")
  }
  if(length(mul) != 1 || mul <= 0)
    stop("In unit name, third entry should be a single positive number")
  u <- list(singular=sing, plural=plur, multiplier=mul)
  if(mul != 1 && (sing=="unit" || plur=="units"))
    stop(paste("A multiplier is not allowed",
               "if the unit does not have a specific name"))
  class(u) <- "units"
  return(u)
}
  
as.units <- function(s) {
  s <- as.list(s)
  n <- length(s)
  if(n > 3)
    stop(paste("Unit name should be a character string,",
               "or a vector/list of 2 character strings,",
               "or a list(character, character, numeric)"))
  
  out <- switch(n+1,
                makeunits(),
                makeunits(s[[1]], s[[1]]),
                makeunits(s[[1]], s[[2]]),
                makeunits(s[[1]], s[[2]], s[[3]]))
  return(out)
}

print.units <- function(x, ...) {
  mul <- x$multiplier
  if(mul == 1)
    cat(paste(x$singular, "/", x$plural, "\n"))
  else 
    cat(paste(mul, x$plural, "\n"))
  return(invisible(NULL))
}

as.character.units <- function(x, ...) {
  mul <- x$multiplier
  return(if(mul == 1) x$plural else paste(mul, x$plural))
}

summary.units <- function(object, ...) {
  x <- object
  scaled <- (x$multiplier != 1)
  named  <- (x$singular != "unit")
  vanilla <- !named && !scaled
  out <-
    if(vanilla) {
      list(legend = NULL,
           axis   = NULL, 
           explain = NULL,
           singular = "unit",
           plural   = "units")
    } else if(named & !scaled) {
      list(legend = paste("Unit of length: 1", x$singular),
           axis   = paren(x$plural, type=spatstat.options('units.paren')),
           explain = NULL,
           singular = x$singular,
           plural   = x$plural)
    } else {
      expanded <- paste(x$multiplier, x$plural)
      expla <- paren(paste("one unit =", expanded),
                     type=spatstat.options('units.paren'))
      list(legend = paste("Unit of length:", expanded),
           axis   = expla, 
           explain  = expla,
           singular = "unit",
           plural   = "units")
    }
  out <- append(out, list(scaled  = scaled,
                          named   = named,
                          vanilla = vanilla))
  class(out) <- "summary.units"
  return(out)
}

print.summary.units <- function(x, ...) {
  if(x$vanilla)
    cat("Unit of length (unnamed)\n")
  else
    cat(paste(x$legend, "\n"))
  invisible(NULL)
}

compatible.units <- function(A, B, ..., coerce=TRUE) {
  stopifnot(inherits(A, "units"))
  if(missing(B)) return(TRUE)
  stopifnot(inherits(B, "units"))
  # check for null units
  Anull <- summary(A)$vanilla
  Bnull <- summary(B)$vanilla
  # `coerce' determines whether `vanilla' units are compatible with other units
  coerce <- as.logical(coerce)
  # 
  agree <- if(!Anull && !Bnull) identical(all.equal(A,B), TRUE) else
           if(Anull && Bnull) TRUE else coerce 
  #
  if(!agree) return(FALSE)
  # A and B agree
  if(length(list(...)) == 0) return(TRUE)
  # recursion
  return(compatible.units(B, ...))
}
