\name{spatstat-deprecated}
\alias{as.psp.owin}
\alias{which.max.im}
\alias{circumradius}
\alias{circumradius.owin}
\alias{circumradius.ppp}
\alias{circumradius.linnet}
\title{Deprecated spatstat functions}
\description{
  Deprecated spatstat functions.
}
\usage{
\method{as.psp}{owin}(x, \dots, window=NULL, 
       check=spatstat.options("checksegments"), fatal=TRUE)
which.max.im(x)
circumradius(x, \dots)
\method{circumradius}{owin}(x, \dots)
\method{circumradius}{ppp}(x, \dots)
\method{circumradius}{linnet}(x, \dots)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat} package.
  
  \code{as.psp.owin} has been replaced by \code{\link{edges}}.
  
  \code{which.max.im(x)} is replaced by
  \code{\link{im.apply}(x, which.max)}.

  \code{circumradius} is replaced by the more appropriately named
  \code{boundingradius}.
}
\keyword{internal}
