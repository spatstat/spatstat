\name{rectcontact}
\alias{rectcontact}
\title{
  Contact Distribution Function using Rectangular Structuring Element
}
\description{
  Computes an estimate of the contact distribution function
  of a set, using a rectangular structuring element.
}
\usage{
rectcontact(X, \dots, asp = 1, npasses=4,
            eps = NULL, r = NULL, breaks = NULL, correction = c("rs", "km"))
}
\arguments{
  \item{X}{
    Logical-valued image. The \code{TRUE} values in the image
    determine the spatial region whose contact distribution function
    should be estimated. 
  }
  \item{\dots}{
    Ignored.
  }
  \item{asp}{
    Aspect ratio for the rectangular metric. A single positive number.
    See \code{\link{rectdistmap}} for explanation.
  }
  \item{npasses}{
    Number of passes to perform in the distance algorithm. A positive integer.
    See \code{\link{rectdistmap}} for explanation.
  }
  \item{eps}{
    Pixel size, if the image should be converted to a finer grid.
  }
  \item{r}{
    Optional vector of distance values. Do Not Use This.
  }
  \item{breaks}{
    Do Not Use This.
  }
  \item{correction}{
    Character vector specifying the edge correction.
  }
}
\details{
  To be written.
}
\value{
  Object of class \code{"fv"}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{Hest}}
}
\examples{
  ## make an image which is TRUE/FALSE inside/outside the letter R
  V <- letterR
  Frame(V) <- grow.rectangle(Frame(V), 0.5)
  Z <- as.im(V, value=TRUE, na.replace=FALSE)
  ## analyse
  plot(rectcontact(Z))
}
\keyword{spatial}
\keyword{nonparametric}

