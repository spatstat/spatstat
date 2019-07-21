\name{spatdim}
\alias{spatdim}
\title{Spatial Dimension of a Dataset}
\description{
  Extracts the spatial dimension of an object in the
  \pkg{spatstat} package.
}
\usage{spatdim(X, intrinsic=FALSE)}
\arguments{
  \item{X}{
    Object belonging to any class defined in the
    \pkg{spatstat} package.
  }
  \item{intrinsic}{
    Logical value indicating whether to return the
    number of intrinsic dimensions. See Details.
  }
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

  Some spatial objects are lower-dimensional subsets of the
  space in which they live. This lower number of dimensions
  is returned if \code{intrinsic=TRUE}. For example, a dataset
  on a linear network
  (an object \code{X} of class \code{"linnet", "lpp", "linim", "linfun"}
  or \code{"lintess"}) returns \code{spatdim(X) = 2} but
  \code{spatdim(X, intrinsic=TRUE) = 1}.
  
  If \code{X} is not a recognised spatial object, the result is \code{NA}.
}
\author{
  \spatstatAuthors.
}
\examples{
  spatdim(lansing)
  A <- osteo$pts[[1]]
  spatdim(A)
  spatdim(domain(A))
  spatdim(chicago)
  spatdim(chicago, intrinsic=TRUE)
}





