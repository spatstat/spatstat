\name{laslett}
\alias{laslett}
\title{
  Laslett's Transform
}
\description{
  Apply Laslett's Transform to a spatial region,
  returning the original and transformed regions,
  and the lower tangent points before and after transformation.
}
\usage{
laslett(X, \dots, verbose = FALSE, plotit = TRUE, discretise = FALSE)
}
\arguments{
  \item{X}{
    Spatial region to be transformed.
    A window (object of class \code{"owin"}) or a logical-valued pixel
    image (object of class \code{"im"}).
  }
  \item{\dots}{
    Graphics arguments to control the plot, passed to
    \code{\link{plot.laslett}} when \code{plotit=TRUE}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
  \item{plotit}{
    Logical value indicating whether to plot the result.
  }
  \item{discretise}{
    Logical value indicating whether polygonal windows should first be
    converted to pixel masks before the Laslett transform is
    computed. This should be set to \code{TRUE} for very complicated
    polygons. 
  }
}
\details{
  To Be Provided.
}
\value{
  A list, which also belongs to the class \code{"laslett"}
  so that it can immediately be printed and plotted.
  
  The list elements are:
  \describe{
    \item{oldX:}{the original dataset \code{X};}
    \item{TanOld:}{a point pattern, whose window is \code{Frame(X)},
      containing the lower tangent points of \code{X};}
    \item{TanNew:}{a point pattern, whose window is the Laslett transform
      of \code{Frame(X)}, and which contains the Laslett-transformed
      positions of the tangent points;}
    \item{Rect:}{a rectangular window, which is the largest rectangle
      lying inside the transformed set;}
    \item{df:}{a data frame giving the locations of the tangent points
      before and after transformation.}
  }
}
\references{
  To Be Provided.
}
\author{
  Kassel Hingee and Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}.
}
\seealso{
  \code{\link{plot.laslett}}
}
\examples{
a <- laslett(heather$coarse)
b <- laslett(discs(runifpoint(30) \%mark\% 0.15))
}
\keyword{spatial}
\keyword{manip}

