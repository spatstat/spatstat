# levelset.R
#
#  $Revision: 1.5 $  $Date: 2015/01/15 07:10:37 $
#
# level set of an image

levelset <- function(X, thresh, compare="<=") {
  # force X and thresh to be evaluated in this frame
  verifyclass(X, "im")
  thresh <- thresh
  switch(compare,
         "<"  = { A <- eval.im(X < thresh) },
         ">"  = { A <- eval.im(X > thresh) },
         "<=" = { A <- eval.im(X <= thresh) },
         ">=" = { A <- eval.im(X >= thresh) },
         "==" = { A <- eval.im(X == thresh) },
         "!=" = { A <- eval.im(X != thresh) },
         stop(paste("unrecognised comparison operator", sQuote(compare))))
  W <- as.owin(eval.im(ifelse1NA(A)))
  return(W)
}

# compute owin containing all pixels where image expression is TRUE

solutionset <- function(..., envir) {
  if(missing(envir))
    envir <- parent.frame()
  A <- try(eval.im(..., envir=envir), silent=TRUE)
  if(inherits(A, "try-error"))
    A <- try(eval(..., envir=envir), silent=TRUE)
  if(inherits(A, "try-error"))
    stop("Unable to evaluate expression")
  if(!is.im(A))
    stop("Evaluating the expression did not yield a pixel image")
  if(A$type != "logical")
    stop("Evaluating the expression did not yield a logical-valued image")
  W <- as.owin(eval.im(ifelse1NA(A)))
  return(W)
}


