\name{lineardisc}
\alias{lineardisc}
\alias{countends}
\title{
  Compute Disc of Given Radius in Linear Network
}
\description{
  Computes the \sQuote{disc} of given radius and centre
  in a linear network.
}
\usage{
  lineardisc(L, x = locator(1), r, plotit = TRUE,
             cols=c("blue", "red","green"))

  countends(L, x = locator(1), r, toler=NULL)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{x}{
    Location of centre of disc.
    Either a point pattern (object of class \code{"ppp"})
    containing exactly 1 point, or a numeric vector of length 2.
  }
  \item{r}{
    Radius of disc.
  }
  \item{plotit}{
    Logical. Whether to plot the disc.
  }
  \item{cols}{
    Colours for plotting the disc. A numeric or character vector of
    length 3 specifying the colours of the disc centre, disc lines and
    disc endpoints respectively.
  }
  \item{toler}{
    Optional. Distance threshold for \code{countends}. See Details.
    There is a sensible default.
  }
}
\details{
  The \sQuote{disc} \eqn{B(u,r)} of centre \eqn{x} and radius \eqn{r}
  in a linear network \eqn{L} is the set of all points
  \eqn{u} in \eqn{L} such that the shortest path distance from \eqn{x}
  to \eqn{u} is less than or equal to \eqn{r}. This is a union of line
  segments contained in \eqn{L}.

  The \emph{relative boundary} of the disc \eqn{B(u,r)}
  is the set of points \eqn{v} such that the shortest path distance from
  \eqn{x} to \eqn{u} is \emph{equal} to \eqn{r}.
  
  The function \code{lineardisc} computes the
  disc of radius \eqn{r} and its relative boundary,
  optionally plots them, and returns them.
  The faster function \code{countends} simply counts the number of
  points in the relative boundary.

  The optional threshold \code{toler} is used to suppress numerical
  errors in \code{countends}.
  If the distance from \eqn{u} to a network vertex \eqn{v}
  is between \code{r-toler} and \code{r+toler}, the vertex
  will be treated as lying on the relative boundary.    
}
\value{
  The value of \code{lineardisc} is a list with two entries:
  \item{lines }{Line segment pattern (object of class \code{"psp"})
    representing the interior disc}
  \item{endpoints}{Point pattern (object of class \code{"ppp"})
    representing the relative boundary of the disc.
  }
  The value of \code{countends} is an integer giving the number of
  points in the relative boundary.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian
}
\seealso{
  \code{\link{linnet}}
}
\references{
Ang, Q.W. (2010)
\emph{Statistical methodology for events on a network}.
Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
Ang, Q.W., Baddeley, A. and Nair, G. (2012)
Geometrically corrected second-order analysis of 
events on a linear network, with applications to
ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
}
\examples{
    # letter 'A' 
    v <- ppp(x=(-2):2, y=3*c(0,1,2,1,0), c(-3,3), c(-1,7))
    edg <- cbind(1:4, 2:5)
    edg <- rbind(edg, c(2,4))
    letterA <- linnet(v, edges=edg)

   lineardisc(letterA, c(0,3), 1.6)
   # count the endpoints
   countends(letterA, c(0,3), 1.6)
   # cross-check (slower)
   en <- lineardisc(letterA, c(0,3), 1.6, plotit=FALSE)$endpoints
   npoints(en)
}
\keyword{spatial}
