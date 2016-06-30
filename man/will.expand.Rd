\name{will.expand}
\alias{will.expand}
\title{
  Test Expansion Rule
}
\description{
  Determines whether an expansion rule will
  actually expand the window or not.
}
\usage{
will.expand(x)
}
\arguments{
  \item{x}{
    Expansion rule.
    An object of class \code{"rmhexpand"}.
  }
}
\details{
  An object of class \code{"rmhexpand"} describes a rule for
  expanding a simulation window. See \code{\link{rmhexpand}}
  for details.

  One possible expansion rule is to do nothing, i.e. not to expand
  the window.
  
  This command inspects the expansion rule \code{x}
  and determines whether it will or will not actually expand the window.
  It returns \code{TRUE} if the window will be expanded.
}
\value{
  Logical value.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{rmhexpand}},
  \code{\link{expand.owin}}
}
\examples{
  x <- rmhexpand(distance=0.2)
  y <- rmhexpand(area=1)
  will.expand(x)
  will.expand(y)
}
\keyword{spatial}
\keyword{manip}
