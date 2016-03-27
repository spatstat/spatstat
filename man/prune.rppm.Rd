\name{prune.rppm}
\alias{prune.rppm}
\title{
  Prune a Recursively Partitioned Point Process Model
}
\description{
  Given a model which has been fitted to point pattern data
  by recursive partitioning, apply pruning to reduce the
  complexity of the partition tree.
}
\usage{
\method{prune}{rppm}(tree, \dots) 
}
\arguments{
  \item{tree}{
    Fitted point process model of class \code{"rppm"}
    produced by the function \code{\link{rppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[rpart]{prune.rpart}}
    to control the pruning procedure.
  }
}
\details{
  This is a method for the generic function \code{\link[rpart]{prune}}
  for the class \code{"rppm"}. An object of this class is a
  point process model, fitted to point pattern data by
  recursive partitioning, by the function \code{\link{rppm}}.
  
  The recursive partition tree will be pruned using
  \code{\link[rpart]{prune.rpart}}. The result is another
  object of class \code{"rppm"}.
}
\value{
  Object of class \code{"rppm"}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{rppm}},
  \code{\link{plot.rppm}},
  \code{\link{predict.rppm}}.
}
\examples{
  # Murchison gold data
  mur <- solapply(murchison, rescale, s=1000, unitname="km")
  mur$dfault <- distfun(mur$faults)
  fit <- rppm(gold ~ dfault + greenstone, data=mur)
  fit
  prune(fit, cp=0.1)
}
\keyword{spatial}
\keyword{models}
