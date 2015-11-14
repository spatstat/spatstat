\name{layout.boxes}
\alias{layout.boxes}
\title{
  Generate a Row or Column Arrangement of Rectangles.
}
\description{
  A simple utility to generate a row or column of boxes (rectangles)
  for use in point-and-click panels.
}
\usage{
layout.boxes(B, n, horizontal = FALSE, aspect = 0.5, usefrac = 0.9)
}
\arguments{
  \item{B}{
    Bounding rectangle for the boxes.
    An object of class \code{"owin"}.
  }
  \item{n}{
    Integer. The number of boxes.
  }
  \item{horizontal}{
    Logical. If \code{TRUE}, arrange the boxes in a horizontal row.
    If \code{FALSE} (the default), arrange them in a vertical column.
  }
  \item{aspect}{
    Aspect ratio (height/width) of each box.
  }
  \item{usefrac}{
    Number between 0 and 1. The 
    fraction of height or width of \code{B} that should be
    occupied by boxes.
  }
}
\details{
  This simple utility generates a list of boxes (rectangles)
  inside the bounding box \code{B} arranged in a regular
  row or column. It is useful for generating the
  positions of the panel buttons in the function \code{\link{simplepanel}}.
}
\value{
  A list of rectangles.
}
\examples{
  B <- owin(c(0,10),c(0,1))
  boxes <- layout.boxes(B, 5, horizontal=TRUE)
  plot(B, main="", col="blue")
  niets <- lapply(boxes, plot, add=TRUE, col="grey")
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{simplepanel}}
}
\keyword{utilities}


