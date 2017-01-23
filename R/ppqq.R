##
##  ppqq.R
##
## P-P and Q-Q versions of fv objects
##

PPversion <- local({

  PPversion <- function(f, theo="theo", columns=".") {
    if(!any(colnames(f) == theo))
      stop(paste(sQuote(theo), "is not the name of a column of f"))
    ## set up inverse theoretical function f_0: 'theo' |-> 'r'
    xname <- fvnames(f, ".x")
    df <- as.data.frame(f)
    theo.table <- df[,theo]
    x.table    <- df[,xname]
    invfun <- approxfun(x=theo.table, y=x.table, rule=1)
    ## evaluate f_0^{-1}(theo) for evenly-spaced grid of 'theo' values
    ra <- range(theo.table)
    theo.seq <- seq(from=ra[1], to=ra[2], length.out=nrow(df))
    x.vals <- invfun(theo.seq)
    ## convert f to a function and evaluate at these 'r' values
    ynames <- setdiff(fvnames(f, columns), theo)
    ff <- as.function(f, value=ynames)
    y.vals <- lapply(ynames, evalselected, x=x.vals, f=ff)
    ## build data frame
    all.vals <- append(list(theo=theo.seq), y.vals)
    names(all.vals) <- c(theo, ynames)
    DF <- as.data.frame(all.vals)
    ## set up fv object
    atr <- attributes(f)
    cnames <- colnames(f)
    i.theo <- match(theo,   cnames)
    i.yval <- match(ynames, cnames)
    ii <- c(i.theo, i.yval)
    old.best <- fvnames(f, ".y")
    best <- if(old.best %in% ynames) old.best else ynames[length(ynames)]
    result <- fv(DF,
                 argu = theo,
                 ylab = atr$ylab,
                 valu = best,
                 fmla = . ~ .x,
                 alim = ra,
                 labl = atr$labl[ii], 
                 desc = atr$desc[ii],
                 unitname = NULL,
                 fname = atr$fname,
                 yexp = atr$yexp)
    fvnames(result, ".") <- c(ynames, theo)
    return(result)
  }

  evalselected <- function(what, f, x){ f(x, what=what) } 

  PPversion
})


QQversion <- function(f, theo="theo", columns=".") {
  if(!any(colnames(f) == theo))
    stop(paste(sQuote(theo), "is not the name of a column of f"))
  ## extract relevant columns of data
  xname <- fvnames(f, ".x")
  ynames <- fvnames(f, columns)
  df <- as.data.frame(f)
  theo.table <- df[,theo]
  x.table    <- df[,xname]
  y.table    <- df[,ynames, drop=FALSE]
  ## set up inverse theoretical function f_0: 'theo' |-> 'r'
  invfun <- approxfun(x=theo.table, y=x.table, rule=1)
  ## apply f_0^{-1} to tabulated function values
  z.table <- as.data.frame(lapply(y.table, invfun))
  ## build data frame
  DF <- cbind(df[,xname,drop=FALSE], z.table)
  ## set up fv object
  atr <- attributes(f)
  cnames <- colnames(f)
  i.x <- match(xname,   cnames)
  i.y <- match(ynames, cnames)
  ii <- c(i.x, i.y)
  old.best <- fvnames(f, ".y")
  best <- if(old.best %in% ynames) old.best else ynames[length(ynames)]
  if(versionstring.spatstat() < package_version("1.38-2")) {
    fvl <- fvlabels(f, expand=TRUE)
    theo.string <- fvl[colnames(f) == theo]
  } else {
    theo.string <- fvlabels(f, expand=TRUE)[[theo]]
  }
  ## remove '(r)' from outer function
  theo.string <- sub(paren(xname), "", theo.string, fixed=TRUE)
  theo.expr <- parse(text=theo.string)
  theo.lang <- theo.expr[[1]]
  ylab <- substitute({{THEO}^{-1}}(FUN),
                     list(FUN=atr$ylab, THEO=theo.lang))
  yexp <- substitute({{THEO}^{-1}}(FUN),
                     list(FUN=atr$yexp, THEO=theo.lang))
  oldlabl <- atr$labl
  labl.iy <- sprintf("{{%s}^{-1}}(%s)",  theo.string, oldlabl[i.y])
  labl.ii <- c(oldlabl[i.x], labl.iy)
  result <- fv(DF,
               argu = atr$argu,
               ylab = ylab,
               valu = best,
               fmla = . ~ .x,
               alim = atr$alim,
               labl = labl.ii,
               desc = atr$desc[ii],
               unitname = NULL,
               fname = atr$fname,
               yexp = yexp)
  fvnames(result, ".") <- ynames
  unitname(result) <- unitname(f)
  return(result)
}
