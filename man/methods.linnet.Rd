\name{methods.linnet}
\alias{methods.linnet} %DoNotExport
\Rdversion{1.1}
\alias{as.linnet}
\alias{as.linnet.linnet}
\alias{as.linnet.lpp}
\alias{as.owin.linnet}
\alias{as.psp.linnet}
\alias{print.linnet}
\alias{summary.linnet}
\alias{unitname.linnet}
\alias{unitname<-.linnet}
\alias{vertexdegree}
\alias{vertices.linnet}
\alias{volume.linnet}
\title{
  Methods for Linear Networks
}
\description{
  These are methods for the class \code{"linnet"} of linear networks.
}
\usage{
as.linnet(X, ...)

\method{as.linnet}{linnet}(X, ..., sparse)

\method{as.linnet}{lpp}(X, ..., fatal=TRUE, sparse)

\method{as.owin}{linnet}(W, ...)

\method{as.psp}{linnet}(x, ..., fatal=TRUE)

\method{print}{linnet}(x, ...)

\method{summary}{linnet}(object, ...)

\method{unitname}{linnet}(x)

\method{unitname}{linnet}(x) <- value

vertexdegree(x)

\method{vertices}{linnet}(w)

\method{volume}{linnet}(x)
}
\arguments{
  \item{x,X,object,w,W}{
    An object of class \code{"linnet"} representing 
    a linear network.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{value}{
    A valid name for the unit of length for \code{x}.
    See \code{\link{unitname}}.
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
  The function \code{as.linnet} is generic.
  It converts data from some other format
  into an object of class \code{"linnet"}.
  The method \code{as.linnet.lpp} extracts the linear network
  information from an \code{lpp} object.

  The other functions are methods for the generic commands
  \code{\link{as.owin}},
  \code{\link{as.psp}},
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link{unitname}},
  \code{\link{unitname<-}},
  \code{\link{vertices}}
  and \code{\link{volume}}
  for the class \code{"linnet"}.

  The method \code{as.owin.linnet} extracts the window containing
  the linear network, and returns it as an object of class \code{"owin"}.
  
  The method \code{as.psp.linnet} extracts the
  lines of the linear network as a line segment pattern (object of class
  \code{"psp"}).

  The method \code{vertices.linnet} extracts the vertices (nodes)
  of the linear network. The function \code{vertexdegree} calculates
  the topological degree of each vertex (the number of lines
  emanating from that vertex) and returns these values as an integer vector.
}
\value{
  For \code{as.linnet} the value is an object of class \code{"linnet"}.
  For other functions, see the help file for the corresponding
  generic function.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{linnet}}
}
\examples{
  data(simplenet)
  simplenet
  summary(simplenet)
  volume(simplenet)
  unitname(simplenet) <- c("cubit", "cubits")
}
\keyword{spatial}
\keyword{methods}
