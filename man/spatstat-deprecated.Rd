\name{spatstat-deprecated}
\alias{as.psp.owin}
\alias{bounding.box}
\alias{clf.test}
\alias{conspire}
\alias{eval.hyper}
\alias{ksmooth.ppp}
\alias{smooth.ppp}
\alias{smooth.fv}
\alias{smooth.msr}
\alias{ks.test.ppm}
\alias{mpl}
\alias{rtoro}
\alias{superimposePSP}
\alias{which.max.im}
\title{Deprecated spatstat functions}
\description{
  Deprecated spatstat functions.
}
\usage{
\method{as.psp}{owin}(x, \dots, window=NULL, 
       check=spatstat.options("checksegments"), fatal=TRUE)
bounding.box(\dots)
clf.test(\dots)
conspire(\dots)
ksmooth.ppp(x, sigma, \dots, edge=TRUE)
smooth.ppp(X, ..., weights = rep(1, npoints(X)), at="pixels")
smooth.fv(x, which = "*", ...,
          method=c("smooth.spline", "loess"),
          xinterval=NULL)
smooth.msr(X, ...)
ks.test.ppm(\dots)
mpl(Q, trend, interaction, data, correction, rbord, use.gam)
rtoro(X, which=NULL, radius=NULL, width=NULL, height=NULL) 
eval.hyper(e, h, simplify=TRUE, ee=NULL)
superimposePSP(\dots, W=NULL, check=TRUE)
which.max.im(x)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat} package.
  
  \code{as.psp.owin} has been replaced by \code{\link{edges}}.
  
  \code{bounding.box} has been replaced by \code{\link{boundingbox}}.
  
  \code{clf.test} has been renamed \code{\link{dclf.test}}.

  \code{conspire} has been replaced by \code{\link{plot.fv}}.

  \code{ksmooth.ppp} has been replaced by \code{\link{density.ppp}}.

  \code{smooth.ppp} has been replaced by \code{\link{Smooth.ppp}}.

  \code{smooth.fv} has been replaced by \code{\link{Smooth.fv}}.

  \code{smooth.msr} has been replaced by \code{\link{Smooth.msr}}.

  \code{mpl} has been replaced by \code{\link{ppm}}.

  \code{ks.test.ppm} has been replaced by
  \code{\link{cdf.test}}.
  
  \code{rtoro} has been replaced by \code{\link{rshift}}.

  \code{eval.hyper} has been replaced by \code{\link{with.hyperframe}}.
  
  \code{superimposePSP} has been replaced by \code{\link{superimpose.psp}}
  which is a method for the generic function \code{\link{superimpose}}.

  \code{which.max.im(x)} is replaced by \code{\link{im.apply}(x, which.max)}.
}
\keyword{internal}
