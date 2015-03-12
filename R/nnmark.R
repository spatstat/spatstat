#
# nnmark.R
#
# $Revision: 1.5 $ $Date: 2014/09/05 06:08:49 $

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
                    result <- as.solist(lapply(mX, function(z) eval.im(z[Y])))
                  },
                  stop("Marks must be a vector or dataframe"))
         },
         points = {
           Y <- nnwhich(X, k=k)
           switch(markformat(X),
                  vector={
                    result <- mX[Y]
                  },
                  dataframe = {
                    result <- mX[Y,, drop=FALSE]
                    row.names(result) <- NULL
                  },
                  stop("Marks must be a vector or dataframe"))
         })
  return(result)
}



  
