\name{cbind.hyperframe}
\alias{cbind.hyperframe}
\alias{rbind.hyperframe}
\title{
  Combine Hyperframes by Rows or by Columns
}
\description{
  Methods for \code{cbind} and \code{rbind} for hyperframes.
}
\usage{
\method{cbind}{hyperframe}(...)
\method{rbind}{hyperframe}(...)
}
\arguments{
  \item{\dots}{
    Any number of hyperframes (objects of class \code{\link{hyperframe}}).
  }
}
\details{
  These are methods for \code{\link{cbind}}
  and \code{\link{rbind}} for hyperframes.
  
  Note that \emph{all} the arguments must be hyperframes (because of
  the peculiar dispatch rules of \code{\link{cbind}} and
  \code{\link{rbind}}).

  To combine a hyperframe with a data frame, one should either 
  convert the data frame to a hyperframe using
  \code{\link{as.hyperframe}}, or explicitly invoke the
  function \code{cbind.hyperframe} or \code{rbind.hyperframe}.
  
  In other words: if \code{h} is a hyperframe and \code{d} is a data frame,
  the result of \code{cbind(h,d)} will be the same as
  \code{cbind(as.data.frame(h), d)}, so that all hypercolumns
  of \code{h} will be deleted (and a warning will be issued).
  To combine \code{h} with \code{d}
  so that all columns of \code{h} are retained,
  type either \code{cbind(h, as.hyperframe(d))} or
  \code{cbind.hyperframe(h,d)}.
}
\value{
  Another hyperframe.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{hyperframe}},
  \code{\link{as.hyperframe}}
}
\examples{
  lambda <- runif(5, min=10, max=30)
  X <- lapply(as.list(lambda), function(x) { rpoispp(x) })
  h <- hyperframe(lambda=lambda, X=X)
  g <- hyperframe(id=letters[1:5], Y=rev(X))
  gh <- cbind(h, g)
  hh <- rbind(h, h)
}
\keyword{spatial}
\keyword{manip}
