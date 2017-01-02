##  terse.R
##
##  code to control terseness and layout of printed output
##
##  $Revision: 1.11 $  $Date: 2016/09/23 02:07:24 $
##


## paragraph break in long output e.g. ppm
parbreak <- function(terse = spatstat.options("terse")) {
  if(waxlyrical('space', terse)) cat("\n")
  return(invisible(NULL))
}

waxlyrical <- local({

  ##  Values of spatstat.options('terse'):
  ##        0    default
  ##        1    suppress obvious wastage e.g. 'gory details'
  ##        2    contract space between paragraphs in long output
  ##        3    suppress extras e.g. standard errors and CI 
  ##        4    suppress error messages eg failed to converge

  TerseCutoff <- list(gory=1,
                      space=2,
                      extras=3,
                      errors=4)

  waxlyrical <- function(type, terse = spatstat.options("terse")) {
    if(!(type %in% names(TerseCutoff)))
      stop(paste("Internal error: unrecognised permission request",
                 sQuote(type)),
           call.=TRUE)
    return(terse < TerseCutoff[[type]])
  }
  
  waxlyrical
  
})

ruletextline <- function(ch="-", n=getOption('width'),
                         terse=spatstat.options('terse')) {
  if(waxlyrical('space', terse)) {
    chn <- paste(rep(ch, n), collapse="")
    chn <- substr(chn, 1, n)
    cat(chn, fill=TRUE)
  }
  return(invisible(NULL))
}
  
