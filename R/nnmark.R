#
# nnmark.R
#
# $Revision: 1.6 $ $Date: 2015/10/21 09:06:57 $

nnmark <- local({

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
                      result <- solapply(mX, lookedup, indeximage=Y)
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

  lookedup <- function(xvals, indeximage) eval.im(xvals[indeximage])

  nnmark
})




  
