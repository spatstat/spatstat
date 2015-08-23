\name{fourierbasis}
\alias{fourierbasis}
\title{Fourier Basis Functions}
\description{Evaluates the Fourier basis functions
  on a \eqn{d}-dimensional box
  with \eqn{d}-dimensional frequencies \eqn{k_i} at the
  \eqn{d}-dimensional coordinates \eqn{x_j}.
}
\usage{
  fourierbasis(x, k, win = boxx(rep(list(0:1), ncol(k))))
}
\arguments{
  \item{x}{
    Coordinates. 
    A \code{data.frame} or matrix with
    \eqn{m} rows and \eqn{d} columns giving
    the \eqn{d}-dimensional coordinates.
  }
  \item{k}{Frequencies.
    A \code{data.frame} or matrix with \eqn{n} rows and \eqn{d} columns
    giving the frequencies of the Fourier-functions.
  }
  \item{win}{
    window (of class \code{"owin"}, \code{"box3"} or \code{"boxx"})
    giving the \eqn{d}-dimensional box domain of the Fourier functions.
  }
}
\details{
  The result is an \eqn{n} by \eqn{m} matrix where the \eqn{(i,j)}'th
  entry is the \eqn{d}-dimensional Fourier basis function with
  frequency \eqn{k_i} evaluated at the point \eqn{x_j}, i.e.,
  \deqn{
    \frac{1}{|W|}
    \exp(2\pi i <k_i,x_j>/|W|)
  }{
    1/|W| * exp(2*pi*i*k_i*x_j/|W|)
  }
  where \eqn{<\cdot,\cdot>}{*} is the \eqn{d}-dimensional inner product
  and \eqn{|W|} is the volume of the 
  domain (window/box). Note that the algorithm does not check whether
  the coordinates given in \code{x} are contained in the given box.
  Actually the box is only 
  used to determine the volume of the domain for normalization.
}
\value{An \code{n} by \code{m} matrix of complex values.}
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
## 27 rows of three dimensional Fourier frequencies:
k <- expand.grid(-1:1,-1:1, -1:1)
## Two random points in the three dimensional unit box:
x <- rbind(runif(3),runif(3))
## 27 by 2 resulting matrix:
v <- fourierbasis(x, k)
head(v)
}
