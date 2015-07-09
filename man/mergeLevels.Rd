\name{mergeLevels}
\alias{mergeLevels}
\title{
  Merge Levels of a Factor
}
\description{
  Specified levels of the factor will be merged into a single level.
}
\usage{
mergeLevels(.f, \dots)
}
\arguments{
  \item{.f}{
    A factor (or a factor-valued pixel image
    or a point pattern with factor-valued marks).
  }
  \item{\dots}{
    List of \code{name=value} pairs, where \code{name} is the
    new merged level, and \code{value} is the vector of old
    levels that will be merged.
  }
}
\details{
  This utility function takes a factor \code{.f}
  and merges specified levels of the factor.

  The grouping is specified by the arguments \code{\dots}
  which must each be given in the form \code{new=old}, where
  \code{new} is the name for the new merged level,
  and \code{old} is a character vector containing the old levels
  that are to be merged.

  The result is a new factor (or factor-valued object),
  in which the levels listed in \code{old}
  have been replaced by a single level \code{new}.

  An argument of the form \code{name=character(0)} or \code{name=NULL}
  is interpreted to mean that all other levels of the old factor
  should be mapped to \code{name}.
}
\value{
  Another factor of the same length as \code{.f}
  (or object of the same kind as \code{.f}).
}
\section{Tips for manipulating factor levels}{
  To remove unused levels from a factor \code{f},
  just type \code{f <- factor(f)}.

  To change the ordering of levels in a factor,
  use \code{\link[base]{factor}(f, levels=l)} or
  \code{\link[stats]{relevel}(f, ref)}.
}
\seealso{
  \code{\link[base]{factor}},
  \code{\link[stats]{relevel}}
}
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
   likert <- c("Strongly Agree", "Agree", "Neutral",
               "Disagree", "Strongly Disagree")
   answers <- factor(sample(likert, 15, replace=TRUE), levels=likert)
   answers
   mergeLevels(answers, Positive=c("Strongly Agree", "Agree"),
                        Negative=c("Strongly Disagree", "Disagree"))
}
\keyword{manip}
\keyword{spatial}
