#'
#'    utiltext.R
#'
#'   Utilities for text output, etc
#'
#'   $Revision: 1.1 $ $Date: 2016/12/30 03:23:09 $
#'

# text magic

commasep <- function(x, join=" and ", flatten=TRUE) {
  px <- paste(x)
  nx <- length(px)
  if(nx <= 1) return(px)
  commas <- c(rep(", ", length(px)-2),
              join,
              "")
  out <- paste0(px, commas, collapse=if(flatten) "" else NULL)
  return(out)
}

paren <- function(x, type="(") {
  if(length(x) == 0) return(x)
  if(identical(type, "")) type <- "blank"
  switch(type,
         "(" = {
           out <- paste("(", x, ")", sep="")
         },
         "[" = {
           out <- paste("[", x, "]", sep="")
         },
         "{" = {
           out <- paste("{", x, "}", sep="")
         },
         blank = {
           out <- paste(x)
         },
         stop(paste("Unrecognised parenthesis type:", sQuote(type)))
         )
  out
}

unparen <- function(x) {
  x <- as.character(x)
  firstchar <- substr(x, 1, 1)
  n <- nchar(x)
  lastchar <- substr(x, n, n)
  enclosed <- n > 2 & (
                       (firstchar == "(" & lastchar == ")") |
                       (firstchar == "[" & lastchar == "]") |
                       (firstchar == "{" & lastchar == "}") )
  if(any(enclosed))
    x[enclosed] <- substr(x[enclosed], 2, n-1)
  return(x)
}

strsplitretain <- local({
  strsplitretain <- function(x, split=",") {
    ## split strings after occurrence of character b, but retain b
    y <- strsplit(x, split)
    lapply(y, addback, b=split)
  }
  addback <- function(x, b=",") {
    n <- length(x)
    if(n <= 1) x else c(paste0(x[-n], b), x[n])
  }    
  strsplitretain
})

truncline <- function(x, nc) {
  if(length(x) > 1)
    return(unlist(lapply(as.list(x), truncline, nc=nc)))
  ## split string into words
  y <- strsplit(x, " ", fixed=TRUE)[[1]]
  ## find max number of whole words that take up nc characters
  maxwords <- max(0, which(cumsum(nchar(y) + 1) <= nc+1))
  if(maxwords == length(y))
    return(x)
  ## truncation will occur.
  pad <- " [..]"
  nc <- nc - nchar(pad)
  maxwords <- max(0, which(cumsum(nchar(y) + 1) <= nc+1))
  z <- paste(y[seq_len(maxwords)], collapse=" ")
  d <- nc - nchar(z)
  if(d < 0)
    z <- substr(z, 1, nc)
  z <- paste0(z, pad)
  return(z)
}

is.blank <- function(s) {
  y <- strsplit(s, "")
  z <- lapply(y, "==", e2=" ")
  ans <- sapply(z, all)
  return(ans)
}

padtowidth <- local({

  blankstring <- function(n) paste(rep(" ", n), collapse="")

  padtowidth <- function(a, b, justify=c("left", "right", "centre")) {
    justify <- match.arg(justify)
    if(is.character(b)) b <- nchar(b) else stopifnot(is.numeric(b))
    extra <- pmax(0, b - nchar(a))
    rpad <- lpad <- ""
    switch(justify,
           left = {
             rpad <- sapply(extra, blankstring)
           },
           right = {
             lpad <- sapply(extra, blankstring)
           },
           centre = {
             lpad <- sapply(floor(extra/2), blankstring)
             rpad <- sapply(ceiling(extra/2), blankstring)
           })
    result <- paste0(lpad, a, rpad)
    return(result)
  }

  padtowidth
})

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


##  grammar, etc

  
ordinal <- function(k) paste(k, ordinalsuffix(k), sep="")

ordinalsuffix <- function(k) {
  last <- abs(k) %% 10
  lasttwo <- abs(k) %% 100
  isteen <- (lasttwo > 10 & lasttwo < 20)
  ending <- ifelse(isteen, "th",
                   ifelse(last == 1, "st",
                          ifelse(last == 2, "nd",
                                 ifelse(last == 3, "rd",
                                        "th"))))
  return(ending)
}

articlebeforenumber <- function(k) {
  k <- abs(k)
  if(k == 11) return("an")
  leading <- floor(k/10^floor(log10(k)))
  if(leading == 8) return("an")
  return("a")
}

numalign <- function(i, nmax, zero="0") {
  stopifnot(i <= nmax)
  nplaces <- as.integer(ceiling(log10(nmax+1)))
  out <- paste(rep(zero, nplaces), collapse="")
  istring <- paste(i)
  ilen <- nchar(istring)
  substr(out, nplaces-ilen+1, nplaces) <- istring
  return(out)
}

singlestring <- function(s, coll="") {
  s <- as.character(s)
  if(length(s) > 1)
    s <- paste(s, collapse=coll)
  return(s)
}

verbalogic <- function(x, op="and") {
  stopifnot(is.character(x))
  istrue <- (x == "TRUE")
  isfalse <- (x == "FALSE")
  isvariable <- !istrue & !isfalse
  y <- x[isvariable]
  switch(op,
         and={
           if(any(isfalse))
             return("FALSE")
           if(all(istrue))
             return("TRUE")
           return(paste(y, collapse=" and "))
         },
         or={
           if(all(isfalse))
             return("FALSE")
           if(any(istrue))
             return("TRUE")
           return(paste(y, collapse=" or "))
         },
         not={
           x[isfalse] <- "TRUE"
           x[istrue] <- "FALSE"
           x[isvariable] <- paste("not {", y, "}")
         },
         stop(paste("Unrecognised operation", sQuote(op))))
}

sensiblevarname <- function(guess, fallback, maxlen=12) {
  out <- if(is.character(guess) &&
            length(guess) == 1  &&
            make.names(guess) == guess) guess else fallback
  out <- substr(out, 1, maxlen)
  return(out)
}

## deparse() can sometimes be equivalent to dumping the whole object
short.deparse <- function(x, maxlen=60) {
  deparse(x,
          nlines=1,
          width.cutoff=maxlen,
          control="delayPromises")
}

## deparse() can produce multiple lines of text
flat.deparse <- function(x) {
  y <- paste(deparse(x), collapse=" ")
  y <- gsub("\n", " ", y)
  y <- gsub(" ", "", y)
  return(y)
}

good.names <- function(nama, defaults, suffices) {
  # ensure sensible, unique names 
  stopifnot(is.character(defaults))
  if(!missing(suffices))
    defaults <- paste(defaults, suffices, sep="")
  result <- nama
  if(is.null(result))
    result <- defaults
  else if(any(blank <- !nzchar(result)))
    result[blank] <- defaults[blank]
  if(anyDuplicated(result))
    result <- make.names(result, unique=TRUE)
  return(result)
}

nzpaste <- function(..., sep=" ", collapse=NULL) {
  # Paste only the non-empty strings
  v <- list(...)
  ok <- sapply(lapply(v, nzchar), any)
  do.call(paste, append(v[ok], list(sep=sep, collapse=collapse)))
}

pasteN <- function(...) {
  # remove NULL arguments then paste
  argh <- list(...)
  argh <- argh[!sapply(argh, is.null)]
  do.call(paste, argh)
}

substringcount <- function(x, y) {
  ## count occurrences of 'x' in 'y'
  yy <- paste0("a", y, "a")
  splot <- strsplit(yy, split=x, fixed=TRUE)
  nhits <- lengths(splot) - 1
  return(nhits)
}

is.parseable <- local({
  is.parseable <- function(x) sapply(x, canparse)

  canparse <- function(z) {
    !inherits(try(parse(text=z), silent=TRUE), "try-error")
  }

  is.parseable
})

make.parseable <- function(x) {
  if(all(is.parseable(x))) x else make.names(x)
}

# paste(expression(..)) seems to be broken

paste.expr <- function(x) {
  unlist(lapply(lapply(x, deparse),
                paste, collapse=""))
}

pasteFormula <- function(f) {
  ## convert formula to a single string
  sf <- paste(format(f), collapse=" ")
  ## remove excessive blanks
  sf <- gsub("   ", " ", sf)
  sf <- gsub("  ", " ", sf)
  return(sf)
}

cat.factor <- function (..., recursive=FALSE) {
  lll <- list(...)
  chk <- sapply(lll,is.factor)
  if(!all(chk))
    stop("First argument is a factor and at least one other argument is not.\n")
  lll <- lapply(lll,as.data.frame,nm="v1")
  return(do.call(rbind,lll)[,1])
}

prange <- function(x) {
  stopifnot(length(x) == 2)
  paren(paste(x, collapse=", "), "[")
}

#   gsub(".", replacement, x) but only when "." appears as a variable

gsubdot <- function(replacement, x) {
  x <- as.character(x)
  stopifnot(length(x) == 1)
  # find all positions of "." in x
  dotpos <- gregexpr("\\.", x)[[1]]
  if(all(dotpos == -1)) return(x)
  # find all positions of "." preceded or followed by alphanumeric
  dotbefore <- gregexpr("\\.[0-9A-Za-z]", x)[[1]]
  dotafter <- gregexpr("[0-9A-Za-z]\\.", x)[[1]] - 1
  # exclude them
  dotpos <- setdiff(dotpos, union(dotbefore, dotafter))
  #
  if(length(dotpos) == 0) return(x)
  lenrep <-length(replacement)
  while(length(dotpos) > 0) {
    dp <- dotpos[1]
    x <- paste0(substr(x, 0, dp-1), replacement, substr(x, dp+1, nchar(x)))
    dotpos <- dotpos[-1] + lenrep-1
  }
  return(x)
}



simplenumber <- local({

  iswhole <- function(x, tol=0) { abs(x %% 1) <= tol }

  iszero <- function(x, tol=0) { abs(x) <= tol }
  
  simplenumber <- function(x, unit = "", multiply="*",
                           tol=.Machine$double.eps) {
    ## Try to express x as a simple multiple or fraction
    if(length(x) > 1)
      return(sapply(as.list(x), simplenumber,
                    unit=unit, multiply=multiply, tol=tol))
    s <- if(x < 0) "-" else ""
    x <- abs(x)
    if(unit == "") {
      if(iswhole(x, tol)) return(paste0(s, round(x)))
      for(i in 1:12) {
        if(iswhole(i/x, tol)) return(paste0(s, i, "/", round(i/x)))
        if(iswhole(i*x, tol)) return(paste0(s, round(i*x), "/", i))
      }
    } else {
      if(iszero(x, tol)) return("0")
      if(iszero(x-1, tol)) return(paste0(s,unit))
      if(iswhole(x, tol)) return(paste0(s, round(x), multiply, unit))
      if(iswhole(1/x, tol)) return(paste0(s, unit, "/", round(1/x)))
      for(i in 2:12) {
        if(iswhole(i/x, tol))
          return(paste0(s, i, multiply, unit, "/", round(i/x)))
        if(iswhole(i*x, tol))
          return(paste0(s, round(i*x), multiply, unit, "/", i))
      }
    }
    return(NULL)
  }

  simplenumber
})


fontify <- function(x, font="italic") {
  if(!nzchar(font) || font == "plain")
    return(x)
  if(is.character(x))
    return(paste0(font, "(", x, ")"))
  if(is.expression(x)) {
    if((n <- length(x)) > 0) {
      for(i in 1:n) 
        x[[i]] <- fontify(x[[i]], font)
    }
    return(x)
  }
  if(is.language(x) || is.numeric(x)) 
    return(substitute(f(X), list(f=as.name(font), X=x)))
  if(all(sapply(x, is.language)))
    return(lapply(x, fontify))
  return(NULL)
}

romansort <- local({

  # sort character strings in order of Roman alphabet
  
  romansort <- function(x) {
    if(!is.character(x)) return(sort(x))
    x <- as.vector(x)
    ## convert each 'word' to a vector of single characters
    cc <- strsplit(x, "")
    ## find position of each character in Roman alphabet
    mm <- lapply(cc, match, table=c(letters, LETTERS))
    mmax <- max(unlist(mm), na.rm=TRUE)
    ## encode
    nn <- sapply(mm, powercode, base=mmax)
    ## find ordering
    oo <- order(nn, na.last=TRUE)
    return(x[oo])
  }

  powercode <- function(x, base) sum(x * base^rev((seq_len(length(x))-1)))

  romansort
})

variablesintext <- function(x) all.vars(as.expression(parse(text=x)))

## convert numeric matrix to character, and blank out lower sub-diagonal.
uptrimat <- function(x) {
  stopifnot(is.matrix(x))
  x[] <- as.character(x)
  x[row(x) > col(x)] <- ""
  return(noquote(x))
}

## convert lty codes to text 
lty2char <- function(i) {
  if(is.numeric(i)) c("blank", "solid", "dashed", "dotted",
                      "dotdash", "longdash", "twodash")[(i %% 7) + 1] else i
}



