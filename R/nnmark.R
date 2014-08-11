#
# nnmark.R
#
# $Revision: 1.3 $ $Date: 2013/07/04 09:17:14 $

nnmark <- function(X, ..., k=1, at=c("pixels", "points")) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  at <- match.arg(at)
  mX <- marks(X)
  switch(at,
         pixels = {
           Y <- nnmap(X, k=k, what="which", ...)
           switch(markformat(X),
                  vector={
                    result <- eval.im(mX[Y])
                  },
                  dataframe = {
                    result <- as.listof(lapply(mX, function(z) eval.im(z[Y])))
                  },
                  stop("Marks must be a vector or dataframe"))
         },
         points = {
           Y <- nnwhich(X, k=k)
           switch(markformat(X),
                  vector={
                    result <- eval.im(mX[Y])
                  },
                  dataframe = {
                    result <- mX[Y,, drop=FALSE]
                  },
                  stop("Marks must be a vector or dataframe"))
         })
  return(result)
}



  
