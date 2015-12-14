\name{deltametric}
\Rdversion{1.1}
\alias{deltametric}
\title{
  Delta Metric
}
\description{
  Computes the discrepancy between two sets \eqn{A} and \eqn{B}
  according to Baddeley's delta-metric.
}
\usage{
deltametric(A, B, p = 2, c = Inf, ...)
}
\arguments{
  \item{A,B}{
    The two sets which will be compared.
    Windows (objects of class \code{"owin"}),
    point patterns (objects of class \code{"ppp"})
    or line segment patterns (objects of class \code{"psp"}).
  }
  \item{p}{
    Index of the \eqn{L^p} metric.
    Either a positive numeric value, or \code{Inf}.
  }
  \item{c}{
    Distance threshold. 
    Either a positive numeric value, or \code{Inf}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution of the distance maps computed by \code{\link{distmap}}.
  }
}
\details{
  Baddeley (1992a, 1992b) defined a distance
  between two sets \eqn{A} and \eqn{B} contained in a space \eqn{W} by 
  \deqn{
    \Delta(A,B) = \left[
    \frac 1 {|W|}
    \int_W
    \left| \min(c, d(x,A)) - \min(c, d(x,B)) \right|^p \, {\rm d}x
    \right]^{1/p}
  }{
    Delta(A,B) = [ (1/|W|) * integral of |min(c, d(x,A))-min(c, d(x,B))|^p dx ]^(1/p)
  }
  where \eqn{c \ge 0}{c >= 0} is a distance threshold parameter,
  \eqn{0 < p \le \infty}{0 < p <= Inf} is the exponent parameter,
  and \eqn{d(x,A)} denotes the 
  shortest distance from a point \eqn{x} to the set \eqn{A}.
  Also \code{|W|} denotes the area or volume of the containing space \eqn{W}.

  This is defined so that it is a \emph{metric}, i.e.
  \itemize{
    \item \eqn{\Delta(A,B)=0}{Delta(A,B)=0} if and only if \eqn{A=B}
    \item \eqn{\Delta(A,B)=\Delta(B,A)}{Delta(A,B)=Delta(B,A)}
    \item \eqn{\Delta(A,C) \le \Delta(A,B) + \Delta(B,C)}{Delta(A,C) <=
      Delta(A,B) + Delta(B,C)}
  }
  It is topologically equivalent to the Hausdorff metric
  (Baddeley, 1992a) but has better stability properties
  in practical applications (Baddeley, 1992b).

  If \eqn{p=\infty}{p=Inf} and \eqn{c=\infty}{c=Inf} the Delta metric
  is equal to the Hausdorff metric.

  The algorithm uses \code{\link{distmap}} to compute the distance maps
  \eqn{d(x,A)} and \eqn{d(x,B)}, then approximates the integral
  numerically.
  The accuracy of the computation depends on the pixel resolution
  which is controlled through the extra arguments \code{\dots} passed
  to \code{\link{as.mask}}.
}
\value{
  A numeric value.
}
\references{
  Baddeley, A.J. (1992a)
  Errors in binary images and an \eqn{L^p} version of the Hausdorff metric.
  \emph{Nieuw Archief voor Wiskunde} \bold{10}, 157--183.

  Baddeley, A.J. (1992b)
  An error metric for binary images.
  In W. Foerstner and S. Ruwiedel (eds)
  \emph{Robust Computer Vision}. Karlsruhe: Wichmann.
  Pages 59--78.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{distmap}}
}
\examples{
  X <- runifpoint(20)
  Y <- runifpoint(10)
  deltametric(X, Y, p=1,c=0.1)
}
\keyword{spatial}
\keyword{math}
