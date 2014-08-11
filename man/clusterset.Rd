\name{clusterset}
\alias{clusterset}
\title{
  Allard-Fraley Estimator of Cluster Feature
}
\description{
  Detect high-density features in a spatial point pattern
  using the (unrestricted) Allard-Fraley estimator.
}
\usage{
  clusterset(X, result=c("marks", "domain"),
            \dots, verbose=TRUE,
            fast=FALSE,
            exact=!fast)
}
\arguments{
  \item{X}{
    A dimensional spatial point pattern (object of class
    \code{"ppp"}).
  }
  \item{result}{
    Character string specifying the type of result. See Details.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
  \item{fast}{
    Logical. If \code{FALSE} (the default), the Dirichlet tile areas
    will be computed exactly using polygonal geometry, so that the
    optimal choice of tiles will be computed exactly.
    If \code{TRUE}, the Dirichlet tile areas
    will be approximated using pixel counting, so the optimal
    choice will be approximate.
  }
  \item{exact}{
    Logical. If \code{TRUE}, the Allard-Fraley estimator
    of the domain will be computed exactly using polygonal geometry.
    If \code{FALSE}, the Allard-Fraley estimator of the domain
    will be approximated by a binary pixel mask.
    The default is initially set to \code{FALSE}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}} to control the
    pixel resolution if \code{exact=FALSE}.
  }
}
\details{
  Allard and Fraley (1997) developed a technique for recognising
  features of high density in a spatial point pattern in the presence of
  random clutter.

  This algorithm computes the \emph{unrestricted} Allard-Fraley estimator.
  The Dirichlet (Voronoi) tessellation of the point pattern \code{X} is
  computed. The smallest \code{m} Dirichlet cells are selected,
  where the number \code{m} is determined by a maximum likelihood
  criterion.
  \itemize{
    \item 
    If \code{fast=FALSE} (the default), the areas of the tiles
    of the Dirichlet tessellation will be computed exactly
    using polygonal geometry. This ensures that the optimal selection of
    tiles is computed exactly.
    \item
    If \code{fast=TRUE}, the Dirichlet tile areas
    will be approximated by counting pixels.
    This is faster, and is usually correct (depending on the pixel
    resolution, which is controlled by the arguments \code{\dots}).
  }

  If \code{result="marks"} the result is the point pattern \code{X}
  with a vector of marks labelling each point with a value \code{yes} or
  \code{no} depending on whether the corresponding Dirichlet cell is
  selected by the Allard-Fraley estimator. In other words each point of
  \code{X} is labelled as either a cluster point or a non-cluster point.

  If \code{result="domain"}, the result is the Allard-Fraley estimator
  of the cluster feature set, which is the union of all the
  selected Dirichlet cells, represented as a window (object of class
  \code{"owin"}).
  \itemize{
    \item
    If \code{exact=TRUE}, the Allard-Fraley set estimator
    will be computed exactly using polygonal geometry.
    The result is a polygonal window.
    \item
    If \code{exact=FALSE}, the Allard-Fraley set estimator
    will be approximated by a binary pixel mask.
    This is faster than the exact computation.
    The result is a binary mask.
  }
  When \pkg{spatstat} is initialised, the default is
  \code{exact=FALSE}.
}
\value{
  If \code{result="marks"}, a multitype point pattern (object of class
  \code{"ppp"}).

  If  \code{result="domain"}, a window (object of class
  \code{"owin"}). 
}
\references{
  Allard, D. and Fraley, C. (1997) 
  Nonparametric maximum likelihood estimation of features in
  spatial point processes using Voronoi tessellation.
  \emph{Journal of the American Statistical Association}
  \bold{92}, 1485--1493.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{nnclean}}, 
  \code{\link{sharpen}}
}
\examples{
  opa <- par(mfrow=c(1,2))
  W <- grow.rectangle(as.rectangle(letterR), 1)
  X <- superimpose(runifpoint(300, letterR),
                   runifpoint(50, W), W=W)
  plot(W, main="clusterset(X)")
  plot(clusterset(X, fast=TRUE), add=TRUE, chars=c("o", "+"), cols=1:2)
  plot(letterR, add=TRUE)
  plot(W, main="clusterset(X, 'd')")
  plot(clusterset(X, "d", exact=FALSE), add=TRUE)
  plot(letterR, add=TRUE)
  par(opa)
}
\keyword{spatial}
\keyword{classif}
