\name{studpermu.test}
\alias{studpermu.test}
\title{
  Studentised Permutation Test
}
\description{
  Perform a studentised permutation test for a difference between
  groups of point patterns.
}
\usage{
   studpermu.test(X, formula, summaryfunction = Kest,
       \dots, rinterval = NULL, nperm = 999,
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
  This function performs the studentized permutation test
  of Hahn (2012) for a difference between groups of point patterns.

  The first argument \code{X} should be either
  \describe{
    \item{a list of lists of point patterns.}{
      Each element of \code{X} will be interpreted as a group of
      point patterns, assumed to be replicates of the same point process.
    }
    \item{a hyperframe:}{
      One column of the hyperframe should contain point patterns,
      and another column should contain a factor indicating the
      grouping. The argument \code{formula} should be a formula in the
      \R language specifying the grouping: it should be of the form
      \code{P ~ G} where \code{P} is the name of the column of point
      patterns, and \code{G} is the name of the factor.
    }
  }
  A group needs to contain at least two point patterns with at least
  \code{minpoints} points in each pattern.

  The function returns an object of class \code{"htest"}
  and \code{"studpermutest"} that can be printed and plotted.
  The printout shows the test result and \eqn{p}-value.
  The plot shows the summary functions for the
  groups (and the group means if requested).
}
\value{
  Object of class \code{"studpermutest"}.
}
\references{
  Hahn, U. (2012) 
  A studentized permutation test for the comparison of
  spatial point patterns.
  \emph{Journal of the American Statistical Association}
  \bold{107} (498), 754--764.
}
\author{
  Ute Hahn.

  Modified for \code{spatstat} by
\adrian


\rolf

and \ege

}
\examples{
  np <- if(interactive()) 99 else 19
  testpyramidal <- studpermu.test(pyramidal, Neurons ~ group, nperm=np)
  testpyramidal
}
\keyword{spatial}
\keyword{htest}
