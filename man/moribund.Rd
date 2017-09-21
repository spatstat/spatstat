\name{moribund} 
\alias{kstest}
\alias{kstest.ppp}
\alias{kstest.ppm}
\alias{kstest.lpp}
\alias{kstest.lppm}
\alias{kstest.slrm}
\alias{plot.kstest}
\alias{bermantest}
\alias{bermantest.ppm}
\alias{bermantest.ppp}
\alias{bermantest.lppm}
\alias{bermantest.lpp}
\title{Outdated Functions} 
\description{
  These outdated functions are retained only for
  compatibility; they will soon be marked as Deprecated.
}
\usage{
kstest(\dots)
kstest.ppp(\dots)
kstest.ppm(\dots)
kstest.lpp(\dots)
kstest.lppm(\dots)
kstest.slrm(\dots)
\method{plot}{kstest}(x, \dots)

bermantest(\dots)
bermantest.ppp(\dots)
bermantest.ppm(\dots)
bermantest.lpp(\dots)
bermantest.lppm(\dots)
}
\arguments{
  \item{x}{
    An object of class \code{"kstest"} or \code{"cdftest"}.
  }
  \item{\dots}{
    Arguments passed to other functions.
  }
}
\details{
  These functions will be Deprecated in future releases of \pkg{spatstat}.

  The \code{kstest} functions have been superseded by
  \code{\link{cdf.test}}.

  The \code{bermantest} functions have been superseded by
  \code{\link{berman.test}}.
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\seealso{
  \code{\link{cdf.test}},
  \code{\link{berman.test}},
  \code{\link{plot.cdftest}}
}
\keyword{spatial}

