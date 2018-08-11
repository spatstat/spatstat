\name{fourierbasis}
\alias{fourierbasis}
\alias{fourierbasisraw}
\title{Fourier Basis Functions}
\description{Evaluates the Fourier basis functions
  on a \eqn{d}-dimensional box
  with \eqn{d}-dimensional frequencies \eqn{k_i} at the
  \eqn{d}-dimensional coordinates \eqn{x_j}.
}
\usage{
  fourierbasis(x, k, win = boxx(rep(list(0:1), ncol(k))))
  fourierbasisraw(x, k, boxlengths)
}
\arguments{
  \item{x}{
    Coordinates. 
    A \code{data.frame} or matrix with
    \eqn{n} rows and \eqn{d} columns giving
    the \eqn{d}-dimensional coordinates.
  }
  \item{k}{Frequencies.
    A \code{data.frame} or matrix with \eqn{m} rows and \eqn{d} columns
    giving the frequencies of the Fourier-functions.
  }
  \item{win}{
    window (of class \code{"owin"}, \code{"box3"} or \code{"boxx"})
    giving the \eqn{d}-dimensional box domain of the Fourier functions.
  }
  \item{boxlengths}{
    numeric giving the side lengths of the box domain of the Fourier functions.
  }
}
\details{
  The result is an \eqn{m} by \eqn{n} matrix where the \eqn{(i,j)}'th
  entry is the \eqn{d}-dimensional Fourier basis function with
  frequency \eqn{k_i} evaluated at the point \eqn{x_j}, i.e.,
  \deqn{
    \frac{1}{\sqrt{|W|}}
    \exp(2\pi i \sum{l=1}^d k_{i,l} x_{j,l}/L_l)
  }{
    1/sqrt(|W|) * exp(2*pi*i*(k_{i,1}*x_{j,1}/L_1 + ... + k_{i,d}*x_{j,d}/L_d))
  }
  where \eqn{L_l}, \eqn{l=1,...,d} are the box side lengths
  and \eqn{|W|} is the volume of the 
  domain (window/box). Note that the algorithm does not check whether
  the coordinates given in \code{x} are contained in the given box.
  Actually the box is only used to determine the side lengths and volume of the
  domain for normalization.
  
  The stripped down faster version \code{fourierbasisraw} doesn't do checking or
  conversion of arguments and requires \code{x} and \code{k} to be matrices.
}
\value{An \code{m} by \code{n} matrix of complex values.}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
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
