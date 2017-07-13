#'
#'     pool.R
#'
#'  $Revision: 1.5 $  $Date: 2017/06/05 10:31:58 $

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


  
