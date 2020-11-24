#'
#'     pool.R
#'
#'     pool      Generic
#'     pool.fv
#'     pool.fasp
#' 
#'  $Revision: 1.6 $  $Date: 2020/11/24 01:37:59 $

pool <- function(...) {
  UseMethod("pool")
}

pool.fv <- local({

  Square <- function(A) { force(A); eval.fv(A^2, relabel=FALSE) }
  Add <- function(A,B){ force(A); force(B); eval.fv(A+B, relabel=FALSE) }
  Cmul <- function(A, f) { force(A); force(f); eval.fv(f * A, relabel=FALSE) }

  pool.fv <- function(..., weights=NULL, relabel=TRUE, variance=TRUE) {
    argh <- list(...)
    n <- narg <- length(argh)
    if(narg == 0) return(NULL)
    if(narg == 1) return(argh[[1]])
    ## validate 
    isfv <- unlist(lapply(argh, is.fv))
    if(!all(isfv))
      stop("All arguments must be fv objects")
    argh <- do.call(harmonise, append(argh, list(strict=TRUE)))
    template <- vanilla.fv(argh[[1]])
    ## compute products
    if(!is.null(weights)) {
      check.nvector(weights, narg, things="Functions")
      Y <- Map(Cmul, argh, weights)
      XY <- Map(Cmul, argh, weights^2)
      sumX <- sum(weights)
      sumX2 <- sum(weights^2)
    } else {
      ## default: weights=1
      Y <- XY <- argh
      sumX <- sumX2 <- narg
    }
    ## sum
    sumY <- Reduce(Add, Y)
    attributes(sumY) <- attributes(template)
    ## ratio-of-sums
    Ratio <- eval.fv(sumY/sumX, relabel=FALSE)
    if(variance) {
      ## variance calculation
      meanX <- sumX/n
      meanY <- eval.fv(sumY/n, relabel=FALSE)
      sumY2 <- Reduce(Add, lapply(Y, Square))
      varX   <- (sumX2 - n * meanX^2)/(n-1)
      varY   <- eval.fv((sumY2 - n * meanY^2)/(n-1), relabel=FALSE)
      sumXY <- Reduce(Add, XY)
      covXY <- eval.fv((sumXY - n * meanX * meanY)/(n-1), relabel=FALSE)
      ## variance by delta method
      relvar <- eval.fv(pmax.int(0, varY/meanY^2 + varX/meanX^2
                                 - 2 * covXY/(meanX * meanY)),
                        relabel=FALSE)
      Variance <- eval.fv(Ratio^2 * relvar/n,
                          relabel=FALSE)
      ## two sigma CI
      hiCI <- eval.fv(Ratio + 2 * sqrt(Variance), relabel=FALSE)
      loCI <- eval.fv(Ratio - 2 * sqrt(Variance), relabel=FALSE)
    }
    ## tweak labels of main estimate
    attributes(Ratio) <- attributes(template)
    if(relabel)
      Ratio <- prefixfv(Ratio,
                        tagprefix="pool",
                        descprefix="pooled ",
                        lablprefix="")
    if(!variance)
      return(Ratio)
    ## tweak labels of variance terms
    attributes(Variance) <- attributes(template)
    Variance <- prefixfv(Variance,
                         tagprefix="var",
                         descprefix="delta-method variance estimate of ",
                         lablprefix="bold(var)~")
    attributes(hiCI) <- attributes(loCI) <-  attributes(template)
    hiCI <- prefixfv(hiCI,
                     tagprefix="hi",
                     descprefix="upper limit of two-sigma CI based on ",
                     lablprefix="bold(hi)~")
    loCI <- prefixfv(loCI,
                     tagprefix="lo",
                     descprefix="lower limit of two-sigma CI based on ",
                     lablprefix="bold(lo)~")
    ## glue together
    result <- Reduce(bind.fv, list(Ratio, Variance, hiCI, loCI))
    ## don't plot variances, by default
    fvnames(result, ".") <- setdiff(fvnames(result, "."),
                                    fvnames(Variance, "."))
    return(result)
  }

  pool.fv
})


##

pool.fasp <- local({

  pool.fasp <- function(...) {
    Alist <- list(...)
    Yname <- short.deparse(sys.call())
    if(nchar(Yname) > 60) Yname <- paste(substr(Yname, 1L, 40L), "[..]")
    nA <-  length(Alist)
    if(nA == 0) return(NULL)
    ## validate....
    ## All arguments must be fasp objects
    notfasp <- !unlist(lapply(Alist, inherits, what="fasp"))
    if(any(notfasp)) {
      n <- sum(notfasp)
      why <- paste(ngettext(n, "Argument", "Arguments"),
                   commasep(which(notfasp)),
                   ngettext(n, "does not", "do not"),
                   "belong to the class",
                   dQuote("fasp"))
      stop(why)
    }
    ## All arguments must have envelopes
    notenv <- !unlist(lapply(Alist, has.env))
    if(any(notenv)) {
      n <- sum(notenv)
      why <- paste(ngettext(n, "Argument", "Arguments"),
                   commasep(which(notenv)),
                   ngettext(n, "does not", "do not"),
                   "contain envelope data")
      stop(why)
    }
  
    if(nA == 1L) return(Alist[[1L]])
  
    ## All arguments must have the same dimensions
    witches <- lapply(Alist, getElement, name="which")
    witch1 <- witches[[1L]]
    same <- unlist(lapply(witches, identical, y=witch1))
    if(!all(same))
      stop("Function arrays do not have the same array dimensions")
  
    ## OK.
    ## Pool envelopes at each position
    result <- Alist[[1L]]
    fns <- result$fns
    for(k in seq_along(fns)) {
      funks <- lapply(Alist, extractfun, k=k)
      fnk <- do.call(pool.envelope, funks)
      attr(fnk, "einfo")$Yname <- Yname
      fns[[k]] <- fnk
    }
    result$fns <- fns
    return(result)
  }

  has.env <- function(z) {
    all(unlist(lapply(z$fns, inherits, what="envelope")))
  }

  extractfun <- function(z, k) { z$fns[[k]] }
  
  pool.fasp
  
})

  
