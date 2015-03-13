\name{studpermu.test}
\alias{studpermu.test}
\title{
  Studentised Permutation Test
}
\description{
  Ute Hahn's studentised permutation test
}
\usage{
   studpermu.test(X, formula, summaryfunction = Kest,
       \dots, rinterval = NULL, nperm = 1000,
        use.Tbar = FALSE, minpoints = 20, rsteps = 128,
        r = NULL, arguments.in.data = FALSE)
}
\arguments{
  \item{X}{
    Data. Either a \code{hyperframe} or a list of lists of point patterns.
  }
  \item{formula}{
    Formula describing the grouping, when \code{X} is a hyperframe.
    The left side of the formula identifies which column of \code{X}
    contains the point patterns.
    The right side identifies the grouping factor. 
    If the formula is missing, the grouping variable is taken to be the
    first column of \code{X} that contains a factor, and the point
    patterns are taken from the first column that contains point patterns.
  }
  \item{summaryfunction}{
    Summary function applicable to point patterns.
  }
  \item{\dots}{
    Additional arguments passed to \code{summaryfunction}.
  }
  \item{rinterval}{
    Interval of distance values \eqn{r} over which the
    summary function should be evaluated and over which the test
    statistic will be integrated. If \code{NULL}, the default
    range of the summary statistic is used (taking the intersection
    of these ranges over all patterns).
  }
  \item{nperm}{
    Number of random permutations for the test.
  }
  \item{use.Tbar}{
    Logical value indicating choice of test statistic.
    If \code{TRUE}, use the alternative test statistic,
    which is appropriate for summary functions with
    roughly constant variance, such as \eqn{K(r)/r} or \eqn{L(r)}.
  }
  \item{minpoints}{
    Minimum permissible number of points in a point pattern
    for inclusion in the test calculation.
  }
  \item{rsteps}{
    Number of discretisation steps in the \code{rinterval}.
  }
  \item{r}{
    Optional vector of distance values as the argument for
    \code{summaryfunction}. Should not usually be given.
    There is a sensible default.
  }
  \item{arguments.in.data}{
    Logical. If \code{TRUE}, individual extra arguments to
    \code{summaryfunction} will be taken from \code{X}
    (which must be a hyperframe). This assumes that 
    the first argument of \code{summaryfunction} is the
    point pattern dataset.
  }
}
\details{
  Studentized permutation test for grouped point patterns.

  Interpreted version, random permutations only.
  
  A group needs to contain at least two point patterns with at least
  \code{minpoints} points in each pattern.

  Returns an object of class \code{"htest"} and \code{"studpermutest"},
  that can be plotted. The plot shows the summary functions for the
  groups (and the group means if requested).
}
\value{
  Object of class \code{"studpermutest"}.
}
\author{
  Ute Hahn.

  Modified by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\examples{
  42
}
\keyword{spatial}
\keyword{htest}
