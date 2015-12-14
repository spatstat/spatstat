\name{as.layered}
\alias{as.layered}
\alias{as.layered.default}
\alias{as.layered.ppp}
\alias{as.layered.splitppp}
\alias{as.layered.solist}
\alias{as.layered.listof}
\alias{as.layered.msr}
\title{Convert Data To Layered Object}
\description{
  Converts spatial data into a layered object.
}
\usage{
 as.layered(X)

 \method{as.layered}{default}(X)

 \method{as.layered}{ppp}(X)

 \method{as.layered}{splitppp}(X)

 \method{as.layered}{solist}(X)

 \method{as.layered}{listof}(X)

 \method{as.layered}{msr}(X)
}
\arguments{
  \item{X}{
    Some kind of spatial data.
  }
}
\value{
  An object of class \code{"layered"} (see \code{\link{layered}}).
}
\details{
  This function converts the object \code{X}
  into an object of class \code{"layered"}.

  The argument \code{X} should contain some kind of spatial data
  such as a point pattern, window, or pixel image.

  If \code{X} is a simple object then it will be converted into
  a \code{layered} object containing only one layer which is equivalent
  to \code{X}.

  If \code{X} can be interpreted as consisting of
  multiple layers of data, then the result will be a \code{layered}
  object consisting of these separate layers of data.
  \itemize{
    \item if \code{X} is a list of class \code{"listof"} or \code{"solist"},
    then \code{as.layered(X)} consists of several layers,
    one for each entry in the list \code{X};
    \item if \code{X} is a multitype point pattern,
    then \code{as.layered(X)} consists of several layers,
    each containing the sub-pattern consisting of points of one type;
    \item if \code{X} is a vector-valued measure,
    then \code{as.layered(X)} consists of several layers,
    each containing a scalar-valued measure.
  }
}
\seealso{
  \code{\link{layered}},
  \code{\link{split.ppp}}
}
\examples{
   as.layered(cells)
   as.layered(amacrine)

   P <- rpoispp(100)
   fit <- ppm(P ~ x+y)
   rs <- residuals(fit, type="score")
   as.layered(rs)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{manip}
 
 
