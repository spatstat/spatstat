#'
#'     pool.R
#'
#'  $Revision: 1.1 $  $Date: 2015/01/25 03:32:22 $

pool <- function(...) {
  UseMethod("pool")
}

pool.fv <- local({

  Square <- function(A) { force(A); eval.fv(A^2) }
  Add <- function(A,B){ force(A); force(B); eval.fv(A+B) }
  Cmul <- function(A, f) { force(A); force(f); eval.fv(f * A) }

  pool.fv <- function(..., weights=NULL) {
    argh <- list(...)
    n <- narg <- length(argh)
    if(narg == 0) return(NULL)
    if(narg == 1) return(argh[[1]])
    ## validate 
    isfv <- unlist(lapply(argh, is.fv))
    if(!all(isfv))
      stop("All arguments must be fv objects")
    argh <- do.call(harmonise, argh)
    template <- vanilla.fv(argh[[1]])
    ## compute products
    if(!is.null(weights)) {
      check.nvector(weights, narg, "Functions")
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
    Ratio <- eval.fv(sumY/sumX)
    ## variance calculation
    meanX <- sumX/n
    meanY <- eval.fv(sumY/n)
    sumY2 <- Reduce(Add, lapply(Y, Square))
    varX   <- (sumX2 - n * meanX^2)/(n-1)
    varY   <- eval.fv((sumY2 - n * meanY^2)/(n-1))
    sumXY <- Reduce(Add, XY)
    covXY <- eval.fv((sumXY - n * meanX * meanY)/(n-1))
    ## variance by delta method
    relvar <- eval.fv(pmax.int(0, varY/meanY^2 + varX/meanX^2
                               - 2 * covXY/(meanX * meanY)))
    Variance <- eval.fv(Ratio^2 * relvar/n)
    ## two sigma CI
    hiCI <- eval.fv(Ratio + 2 * sqrt(Variance))
    loCI <- eval.fv(Ratio - 2 * sqrt(Variance))
    ## relabel
    attributes(Ratio) <- attributes(Variance) <- attributes(template)
    Ratio <- prefixfv(Ratio,
                      tagprefix="pool",
                      descprefix="pooled ",
                      lablprefix="")
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
    ##
    result <- Reduce(bind.fv, list(Ratio, Variance, hiCI, loCI))
    return(result)
  }

  pool.fv
})


  
