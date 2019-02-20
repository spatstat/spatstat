#
# Functions for extracting and setting the name of the unit of length
#
#   $Revision: 1.29 $   $Date: 2019/02/20 03:34:50 $
#
#

unitname <- function(x) {
  UseMethod("unitname")
}

unitname.owin <- function(x) {
  u <- as.unitname(x$units)
  return(u)
}

unitname.ppp <- function(x) {
  u <- as.unitname(x$window$units)
  return(u)
}

unitname.im <- function(x) {
  u <- as.unitname(x$units)
  return(u)
}

unitname.default <- function(x) {
  return(as.unitname(attr(x, "units")))
}

"unitname<-" <- function(x, value) {
  UseMethod("unitname<-")
}

"unitname<-.owin" <- function(x, value) {
  x$units <- as.unitname(value)
  return(x)
}

"unitname<-.ppp" <- function(x, value) {
  w <- x$window
  unitname(w) <- value
  x$window <- w
  return(x)
}

"unitname<-.im" <- function(x, value) {
  x$units <- as.unitname(value)
  return(x)
}

"unitname<-.default" <- function(x, value) {
  if(is.null(x)) return(x)
  attr(x, "units") <- as.unitname(value)
  return(x)
}


###  class 'unitname'

makeunitname <- function(sing="unit", plur="units", mul = 1) {
  if(!is.character(sing))
    stop("In unit name, first entry should be a character string")
  if(!is.character(plur))
    stop("In unit name, second entry should be a character string")
  mul <- try(as.numeric(mul), silent=TRUE)
  if(inherits(mul, "try-error"))
    stop("In unit name, third entry should be a number")
  if(length(mul) != 1 || mul <= 0)
    stop("In unit name, third entry should be a single positive number")
  u <- list(singular=sing, plural=plur, multiplier=mul)
  if(mul != 1 && (sing=="unit" || plur=="units"))
    stop(paste("A multiplier is not allowed",
               "if the unit does not have a specific name"))
  class(u) <- "unitname"
  return(u)
}
  
as.unitname <- function(s) {
  if(inherits(s, "unitname")) return(s)
  s <- as.list(s)
  n <- length(s)
  if(n > 3)
    stop(paste("Unit name should be a character string,",
               "or a vector/list of 2 character strings,",
               "or a list(character, character, numeric)"))
  out <- switch(n+1,
                makeunitname(),
                makeunitname(s[[1]], s[[1]]),
                makeunitname(s[[1]], s[[2]]),
                makeunitname(s[[1]], s[[2]], s[[3]]))
  return(out)
}

print.unitname <- function(x, ...) {
  mul <- x$multiplier
  if(mul == 1)
    cat(paste(x$singular, "/", x$plural, "\n"))
  else 
    cat(paste(mul, x$plural, "\n"))
  return(invisible(NULL))
}

as.character.unitname <- function(x, ...) {
  mul <- x$multiplier
  return(if(mul == 1) x$plural else paste(mul, x$plural))
}

is.vanilla <- function(u) {
  u <- as.unitname(u)
  z <- (u$singular == "unit") && (u$multiplier == 1)
  return(z)
}

summary.unitname <- function(object, ...) {
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
  class(out) <- "summary.unitname"
  return(out)
}

print.summary.unitname <- function(x, ...) {
  if(x$vanilla)
    cat("Unit of length (unnamed)\n")
  else
    cat(paste(x$legend, "\n"))
  invisible(NULL)
}

compatible.unitname <- function(A, B, ..., coerce=TRUE) {
  A <- as.unitname(A)
  if(missing(B)) return(TRUE)
  B <- as.unitname(B)
  # check for null units
  Anull <- summary(A)$vanilla
  Bnull <- summary(B)$vanilla
  # `coerce' determines whether `vanilla' units are compatible with other units
  coerce <- as.logical(coerce)
  # 
  agree <- if(!Anull && !Bnull) isTRUE(all.equal(A,B)) else
           if(Anull && Bnull) TRUE else coerce 
  #
  if(!agree) return(FALSE)
  # A and B agree
  if(length(list(...)) == 0) return(TRUE)
  # recursion
  return(compatible.unitname(B, ...))
}

harmonize.unitname <-
harmonise.unitname <- function(..., coerce=TRUE, single=FALSE) {
  argh <- list(...)
  n <- length(argh)
  if(n == 0) return(NULL)
  u <- lapply(argh, as.unitname)
  if(n == 1) return(if(single) u[[1L]] else u)
  if(coerce) {
    #' vanilla units are compatible with another unit
    s <- lapply(u, summary)
    v <- sapply(s, getElement, name="vanilla")
    if(all(v))
      return(if(single) u[[1L]] else u)
    u <- u[!v]
  }
  z <- unique(u)
  if(length(z) > 1) stop("Unitnames are incompatible", call.=FALSE)
  if(single) return(z[[1]])
  z <- rep(z, n)
  names(z) <- names(argh)
  return(z)
}

# class 'numberwithunit':  numeric value(s) with unit of length

numberwithunit <- function(x, u) {
  u <- as.unitname(u)
  x <- as.numeric(x)
  unitname(x) <- u
  class(x) <- c(class(x), "numberwithunit")
  return(x)
}

"%unit%" <- function(x, u) {
  numberwithunit(x, u)
}

format.numberwithunit <- function(x, ..., collapse=" x ", modifier=NULL) {
  u <- summary(unitname(x))
  uname <- if(all(x == 1)) u$singular else u$plural
  y <- format(as.numeric(x), ...)
  z <- pasteN(paste(y, collapse=collapse), 
              modifier, uname, u$explain)
  return(z)
}

as.character.numberwithunit <- function(x, ...) {
  return(format(x))
}

print.numberwithunit <- function(x, ...) {
  cat(format(x, ...), fill=TRUE)
  return(invisible(NULL))
}


    
