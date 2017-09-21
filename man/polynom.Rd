\name{polynom}
\alias{polynom}
\title{
  Polynomial in One or Two Variables
}
\description{
  This function is used to represent a polynomial term in a model
  formula.
  It computes the homogeneous terms in the polynomial of degree \code{n}
  in one variable \code{x} or two variables \code{x,y}.
}
\usage{
  polynom(x, \dots)
}
\arguments{
  \item{x}{
    A numerical vector.
  }
  \item{\dots}{
    Either a single integer \code{n} specifying the degree of the polynomial,
    or two arguments \code{y,n} giving another vector of data \code{y}
    and the degree of the polynomial.
  }
}
\details{
  This function is typically used inside a model formula
  in order to specify the most general possible polynomial
  of order \code{n} involving one numerical variable \code{x}
  or two numerical variables \code{x,y}.

  It is equivalent to \code{\link[stats]{poly}(, raw=TRUE)}.

  If only one numerical vector argument \code{x} is given, 
  the function computes the vectors \code{x^k} for
  \code{k = 1, 2, \dots, n}. These vectors are combined into a matrix
  with \code{n} columns.
  
  If two numerical vector arguments \code{x,y} are given,
  the function computes the vectors \code{x^k * y^m} for
  \code{k >= 0} and \code{m >= 0} satisfying
  \code{0 < k + m <= n}. These vectors are combined into a matrix
  with one column for each homogeneous term.
}
\value{
   A numeric matrix, with rows corresponding to the entries of \code{x},
   and columns corresponding to the terms in the polynomial.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[stats]{poly}},
  \code{\link{harmonic}}
}
\examples{
   x <- 1:4
   y <- 10 * (0:3)
   polynom(x, 3)
   polynom(x, y, 3)
}
\keyword{arith}
