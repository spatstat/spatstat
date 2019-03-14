\name{nestsplit}
\alias{nestsplit}
\title{
  Nested Split
}
\description{
  Applies two splitting operations to a point pattern,
  producing a list of lists of patterns.
}
\usage{
  nestsplit(X, \dots)
}
\arguments{
  \item{X}{
    Point pattern to be split.
    Object of class \code{"ppp"}.
  }
  \item{\dots}{
    Data determining the splitting factors or splitting regions.
    See Details.
  }
}
\details{
  This function splits the point pattern \code{X} into several
  sub-patterns using \code{\link{split.ppp}}, then splits each of the
  sub-patterns into sub-sub-patterns using \code{\link{split.ppp}}
  again.
  The result is a hyperframe containing the sub-sub-patterns
  and two factors indicating the grouping.

  The arguments \code{\dots} determine the two splitting factors
  or splitting regions. Each argument may be:
  \itemize{
    \item a factor (of length equal to the number of points in \code{X})
    \item the name of a column of marks of \code{X} (provided this
    column contains factor values)
    \item a tessellation (class \code{"tess"})
    \item a pixel image (class \code{"im"}) with factor values
    \item a window (class \code{"owin"})
    \item identified by name (in the form \code{name=value})
    as one of the formal arguments of \code{\link{quadrats}} or
    \code{\link{tess}}
  }
  The arguments will be processed to yield a list of two
  splitting factors/tessellations. The splits will be applied
  to \code{X} consecutively to produce the sub-sub-patterns.
}
\value{
  A hyperframe with three columns. The first column contains the
  sub-sub-patterns. The second and third columns are factors
  which identify the grouping according to the two splitting factors.
}
\author{
  Original idea by Ute Hahn.
  Code by \spatstatAuthors.
}
\seealso{
  \code{\link{split.ppp}},
  \code{\link{quantess}}
}
\examples{
  # factor and tessellation
  Nft <- nestsplit(amacrine, marks(amacrine), quadrats(amacrine, 3, 1))
  Ntf <- nestsplit(amacrine, quadrats(amacrine, 3, 1), marks(amacrine))
  Ntf

  # two factors
  big <- with(marks(betacells), area > 300)
  Nff <- nestsplit(betacells, "type", factor(big))

  # two tessellations
  Tx <- quantess(redwood, "x", 4)
  Td <- dirichlet(runifpoint(5, Window(redwood)))
  Ntt <- nestsplit(redwood, Td, Tx)
  Ntt2 <- nestsplit(redwood, Td, ny=3)
}
\keyword{spatial}
\keyword{manip}
