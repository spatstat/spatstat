\name{as.data.frame.hyperframe}
\alias{as.data.frame.hyperframe}
\title{Coerce Hyperframe to Data Frame}
\description{
  Converts a hyperframe to a data frame.
}
\usage{
\method{as.data.frame}{hyperframe}(x, row.names = NULL,
                                  optional = FALSE, ..., 
                                  discard=TRUE, warn=TRUE)
}
\arguments{
  \item{x}{Hyperframe (object of class \code{"hyperframe"}).}
  \item{row.names}{Optional character vector of row names.}
  \item{optional}{Argument passed to \code{\link{as.data.frame}}
    controlling what happens to row names.}
  \item{\dots}{Ignored.}
  \item{discard}{Logical. Whether to discard columns of the hyperframe
    that do not contain atomic data. See Details.
  }
  \item{warn}{Logical. Whether to issue a warning when columns are discarded.}
}
\details{
  This is a method for the generic function \code{\link{as.data.frame}}
  for the class of hyperframes (see \code{\link{hyperframe}}.

  If \code{discard=TRUE}, any columns of the hyperframe that
  do not contain atomic data will be removed (and a warning will
  be issued if \code{warn=TRUE}). 
  If \code{discard=FALSE}, then such columns are converted to
  strings indicating what class of data
  they originally contained.
}
\value{
  A data frame.
}
\examples{
  h <- hyperframe(X=1:3, Y=letters[1:3], f=list(sin, cos, tan))
  as.data.frame(h, discard=TRUE, warn=FALSE)
  as.data.frame(h, discard=FALSE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
