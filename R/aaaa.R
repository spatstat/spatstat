#'
#'    aaaa.R
#'
#'   Code that must be read before the rest of the R code in spatstat
#' 
#'    $Revision: 1.3 $  $Date: 2014/12/07 11:27:51 $

#' ...................................................................
#'   intermaker:
#'   Class structure for functions like 'Strauss'
#'   so they print a nice description.
#'

intermaker <- function(f, blank) {
  # f is the creator function like 'Strauss'
  class(f) <- c("intermaker", class(f))
  # blank is the prototype interaction object
  attr(f, "b") <- blank[c("creator", "name", "par", "parnames")]
  return(f)
}

print.intermaker <- function(x, ...) {
  b <- attr(x, "b")
  argh <- names(formals(x))
  if(length(argh) > 0) {
    desc <- b$parnames
    names(desc) <- names(b$par)
    explain <- paste(", where",
                     commasep(paste(sQuote(argh), "is", desc[argh])))
  } else explain <- NULL
  blah <- paste0("Function ",
                 b$creator,
                 paren(paste(argh, collapse=", ")), 
                 ": creates the interpoint interaction of the ",
                 b$name,
                 explain)
  splat(blah)
  return(invisible(NULL))
}
