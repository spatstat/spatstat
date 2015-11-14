\name{as.linnet.linim}
\alias{as.linnet.lpp}
\alias{as.linnet.linim}
\alias{as.linnet.linfun}
\title{
  Extract Linear Network from Data on a Linear Network
}
\description{
  Given some kind of data on a linear network,
  the command \code{as.linnet} extracts the linear network itself.
}
\usage{
 \method{as.linnet}{linim}(X, \dots)

 \method{as.linnet}{linfun}(X, \dots)

 \method{as.linnet}{lpp}(X, \dots, fatal=TRUE, sparse)
}
\arguments{
  \item{X}{
    Data on a linear network.
    A point pattern (class \code{"lpp"}),
    pixel image (class \code{"linim"}) or function (class \code{"linfun"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{fatal}{
    Logical value indicating whether data in the wrong format
    should lead to an error (\code{fatal=TRUE}) or a warning
    (\code{fatal=FALSE}).
  }
  \item{sparse}{
    Logical value indicating whether to use a sparse matrix
    representation, as explained in \code{\link{linnet}}.
    Default is to keep the same representation as in \code{X}.
  }
}
\details{
  These are methods for the generic \code{\link{as.linnet}}
  for various classes.
  
  The network on which the data are defined is extracted.
}
\value{
  A linear network (object of class \code{"linnet"}).
}
\author{
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link{linnet}}, 
  \code{\link{methods.linnet}}.
}
\examples{
  # make some data
  xcoord <- linfun(function(x,y,seg,tp) { x }, simplenet)
  as.linnet(xcoord)
  X <- as.linim(xcoord)
  as.linnet(X)
}
\keyword{spatial}
\keyword{manip}
