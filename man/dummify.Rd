\name{dummify}
\alias{dummify}
\title{
  Convert Data to Numeric Values by Constructing Dummy Variables
}
\description{
  Converts data of any kind to numeric values.
  A factor is expanded to a set of dummy variables. 
}
\usage{
dummify(x)
}
\arguments{
  \item{x}{
    Vector, factor, matrix or data frame to be converted.
  }
}
\details{
  This function converts data (such as a factor) to numeric values
  in order that the user may calculate, for example, 
  the mean, variance, covariance and correlation of the data.

  If \code{x} is a numeric vector or integer vector, it is returned
  unchanged.

  If \code{x} is a logical vector, it is converted to a 0-1 matrix with
  2 columns. The first column contains a 1 if the logical value is
  \code{FALSE}, and the second column contains a 1 if the logical
  value is \code{TRUE}.

  If \code{x} is a complex vector, it is converted to a matrix with 2
  columns, containing the real and imaginary parts.

  If \code{x} is a factor, the result is a matrix of 0-1 dummy
  variables. The matrix has one column for each possible level of the
  factor. The \code{(i,j)} entry is 
  equal to 1 when the \code{i}th factor value equals the
  \code{j}th level, and is equal to 0 otherwise.

  If \code{x} is a matrix or data frame, the appropriate conversion is
  applied to each column of \code{x}.

  Note that, unlike \code{\link[stats]{model.matrix}}, this command converts a
  factor into a full set of dummy variables (one column for each level of
  the factor).
}
\value{
  A numeric matrix.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\examples{
   chara <- sample(letters[1:3], 8, replace=TRUE)
   logi <- (runif(8) < 0.3)
   comp <- round(4*runif(8) + 3*runif(8) * 1i, 1)
   nume <- 8:1 + 0.1
   df <- data.frame(nume, chara, logi, comp)
   df
   dummify(df)
}
\keyword{math}
