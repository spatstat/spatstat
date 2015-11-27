##  terse.R
##
##  code to control terseness and layout of printed output
##
##  $Revision: 1.9 $  $Date: 2015/11/27 06:59:57 $
##

## 'split cat'
## Replacement for 'cat(paste(...))' ensuring a multi-word output string
## doesn't extend over text margin

splat <- function(..., indent=0) {
  s <- paste(...)
  ## split at newline characters, if present
  ss <- unlist(strsplit(s, "\n"))
  if(indent == 0) {
    for(ssi in ss) 
      cat(unlist(strsplit(ssi, " ")), fill=TRUE)
  } else {
    ispace <- paste(rep(" ", indent-1), collapse="")
    for(ssi in ss) 
      cat(ispace, unlist(strsplit(ssi, " ")), fill=TRUE)
  }
  return(invisible(NULL))
}

choptext <- local({

  choptext <- function(..., prefix="") {
    s <- paste(...)
    ## split at newline characters, if present
    lines <- unlist(strsplit(s, "\n"))
    ## cut into pieces that don't overreach width
    w <- getOption('width')
    lines <- sapply(lines, chopline, w=w, prefix=prefix)
    return(lines)
  }

  chopline <- function(s, w=getOption('width'), prefix="") {
    words <- unlist(strsplit(s, " "))
    nwords <- length(words)
    wordlengths <- nchar(words)
    lines <- character(0)
    prefixlength <- nchar(prefix)
    while(nwords > 0) {
      wordends <- cumsum(wordlengths + c(prefixlength, rep(1, nwords-1)))
      n <- which.max(wordends * (wordends <= w))
      if(n == 0) n <- 1
      lines <- c(lines, paste(words[1:n], collapse=" "))
      words <- words[-(1:n)]
      wordlengths <- wordlengths[-(1:n)]
      nwords <- nwords - n
      prefixlength <- 0
    }
    return(lines)
  }
                         
  choptext
})

exhibitStringList <- function(prefix, strings) {
  shortblurb <- paste(prefix, paste(strings, collapse=", "), "\n")
  if(nchar(shortblurb) < options("width")[[1]]) {
    cat(shortblurb)
  } else {
    cat(paste(prefix,"\n"))
    splat("  ", paste(strings, collapse=" "))
  }
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

ruletextline <- function(ch="-", n=getOption('width'),
                         terse=spatstat.options('terse')) {
  if(waxlyrical('space', terse)) {
    chn <- paste(rep(ch, n), collapse="")
    chn <- substr(chn, 1, n)
    cat(chn, fill=TRUE)
  }
  return(invisible(NULL))
}
  
