##  terse.R
##
##  code to control terseness and layout of printed output
##
##  $Revision: 1.5 $  $Date: 2014/10/24 00:22:30 $
##

## 'split cat'
## Replacement for 'cat(paste(...))' ensuring a multi-word output string
## doesn't extend over text margin

splat <- function(...) {
  s <- paste(...)
  ## split at newline characters, if present
  ss <- unlist(strsplit(s, "\n"))
  for(ssi in ss) 
    cat(unlist(strsplit(ssi, " ")), fill=TRUE)
  return(invisible(NULL))
}

pasteFormula <- function(f) {
  ## convert formula to a single string
  sf <- paste(format(f), collapse=" ")
  ## remove excessive blanks
  sf <- gsub("   ", " ", sf)
  sf <- gsub("  ", " ", sf)
  return(sf)
}

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

  
