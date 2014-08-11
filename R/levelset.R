# levelset.R
#
#  $Revision: 1.3 $  $Date: 2008/07/16 17:39:32 $
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
  W <- as.owin(eval.im(ifelse(A, 1, NA)))
  return(W)
}

# compute owin containing all pixels where image expression is TRUE

solutionset <- function(..., envir) {
  if(missing(envir))
    envir <- parent.frame()
  A <- eval.im(..., envir=envir)
  if(A$type != "logical")
    stop("Evaluating the expression did not yield a logical-valued image")
  W <- as.owin(eval.im(ifelse(A, 1, NA)))
  return(W)
}


