\name{spatdim}
\alias{spatdim}
\title{Spatial Dimension of a Dataset}
\description{
  Extracts the spatial dimension of an object in the
  \pkg{spatstat} package.
}
\usage{spatdim(X)}
\arguments{
  \item{X}{Object belonging to any class defined in the
    \pkg{spatstat} package.}
}
\value{
  An integer, or \code{NA}.
}
\details{
  This function returns the number of spatial coordinate dimensions
  of the dataset \code{X}. The results for some of the more common
  types of objects are as follows:
  \tabular{ll}{
    \bold{object class} \tab \bold{dimension} \cr
    \code{"ppp"} \tab 2 \cr
    \code{"lpp"} \tab 2 \cr
    \code{"pp3"} \tab 3 \cr
    \code{"ppx"} \tab number of \emph{spatial} dimensions \cr
    \code{"owin"} \tab 2 \cr
    \code{"psp"} \tab 2 \cr
    \code{"ppm"} \tab 2
  }
  Note that time dimensions are not counted.

  If \code{X} is not a recognised spatial object, the result is \code{NA}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
  spatdim(lansing)
}





