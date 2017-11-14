## scriptUtils.R
##       $Revision: 1.6 $ $Date: 2017/11/14 06:42:02 $

## slick way to use precomputed data
##    If the named file exists, it is loaded, giving access to the data.
##    Otherwise, 'expr' is evaluated, and all objects created
##    are saved in the designated file, for loading next time.

reload.or.compute <- function(filename, expr, 
                              objects=NULL,
                              destination=parent.frame(),
                              force=FALSE) {
  stopifnot(is.character(filename) && length(filename) == 1)
  if(force || !file.exists(filename)) {
    ## evaluate 'expr' in a fresh environment
    ee <- as.expression(substitute(expr))
    en <- new.env()
    local(eval(ee), envir=en)
    ## default is to save all objects that were created
    if(is.null(objects))
      objects <- ls(envir=en)
    ## save them in the designated file
    evalq(save(list=objects, file=filename, compress=TRUE), envir=en)
    ## assign them into the parent frame 
    for(i in seq_along(objects))
      assign(objects[i], get(objects[i], envir=en), envir=destination)
    result <- objects
  } else {
    result <- load(filename, envir=destination)
    if(!all(ok <- (objects %in% result))) {
      nbad <- sum(!ok)
      warning(paste(ngettext(nbad, "object", "objects"),
                    commasep(sQuote(objects[!ok])),
                    ngettext(nbad, "was", "were"),
                    "not present in data file", dQuote(filename)),
              call.=FALSE)
    }
  }
  return(invisible(result))
}
