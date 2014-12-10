#'
#'    aaaa.R
#'
#'   Code that must be read before the rest of the R code in spatstat
#' 
#'    $Revision: 1.4 $  $Date: 2014/12/10 10:34:53 $

#' ...................................................................
#'   intermaker:
#'   Class structure for functions like 'Strauss'
#'   so they print a nice description.
#'

intermaker <- function(f, blank) {
  # f is the creator function like 'Strauss'
  class(f) <- c("intermaker", class(f))
  # blank is the prototype interaction object: extract some fields
  desired <- c("creator", "name", "par", "parnames", "pardesc")
  avail <- desired[desired %in% names(blank)]
  attr(f, "b") <- blank[avail]
  return(f)
}

print.intermaker <- function(x, ...) {
  b <- attr(x, "b")
  argh <- names(formals(x))
  explain <- NULL
  if(length(argh) > 0) {
    desc <- b$pardesc %orifnull% b$parnames
    namep <- names(b$par)
    if(length(desc) == length(namep) && all(argh %in% namep)) {
      names(desc) <- namep
      explain <- paste(", where",
                       commasep(paste(sQuote(argh), "is the", desc[argh])))
    }
  }
  blah <- paste0("Function ",
                 b$creator,
                 paren(paste(argh, collapse=", ")), 
                 ": creates the interpoint interaction of the ",
                 b$name,
                 explain)
  splat(blah)
  return(invisible(NULL))
}
