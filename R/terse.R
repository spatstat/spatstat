##  terse.R
##
##  code to control terseness and layout of printed output
##
##  $Revision: 1.11 $  $Date: 2016/09/23 02:07:24 $
##

## 'split cat'
## Replacement for 'cat(paste(...))' ensuring a multi-word output string
## doesn't extend over text margin

splat <- function(..., indent=0) {
  st <- pasteN(...) # removes NULL arguments without making a space
  ## split at newline characters, if present
  ss <- unlist(strsplit(st, "\n"))
  ## 
  if(is.numeric(indent)) {
    nindent <- indent
    indent <- paste(rep(" ", nindent), collapse="")
  } else if(is.character(indent)) {
    nindent <- nchar(indent)
  } else stop("indent should be character or numeric")
  w <- getOption('width')
  if(nindent >= w) {
    warning("indentation is more than the permissible text width: ignored")
    nindent <- 0
  }
  ##
  if(nindent == 0) {
    for(ssi in ss) 
      cat(unlist(strsplit(ssi, " ")), fill=TRUE)
  } else {
    wfill <- w - nindent
    for(ssi in ss) {
      vi <- choptextline(ssi, w=w, indent=indent)
      for(vij in vi) {
        cat(indent)
        cat(vij, fill=wfill)
      }
    }
  }
  return(invisible(NULL))
}

choptext <- function(..., prefix="", indent="") {
  s <- paste(...)
  ## split at newline characters, if present
  lines <- unlist(strsplit(s, "\n"))
  ## cut into pieces that don't overreach width
  w <- getOption('width')
  lines <- sapply(lines, choptextline, w=w, prefix=prefix, indent=indent)
  lines <- unname(as.vector(lines))
  return(lines)
}

choptextline <- function(st, w=getOption('width'), prefix="", indent="") {
  words <- unlist(strsplit(st, " "))
  nwords <- length(words)
  wordlengths <- nchar(words)
  lines <- character(0)
  prefixlength <- nchar(prefix)
  indentlength <- nchar(indent)
  while(nwords > 0) {
    inset <- prefixlength + indentlength
    wordends <- cumsum(wordlengths + c(inset, rep(1, nwords-1)))
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
  
