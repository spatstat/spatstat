\name{rdpp}
\alias{rdpp}
\title{Simulation of a Determinantal Point Process}
\description{
  Generates simulated realisations from a determinantal point process.
}
\usage{
  rdpp(eig, index, basis = "fourierbasis",
       window = boxx(rep(list(0:1), ncol(index))),
       reject_max = 10000, progress = 0, debug = FALSE, \dots)
}
\arguments{
  \item{eig}{
    vector of values between 0 and 1 specifying the non-zero
    eigenvalues for the process.
  }
  \item{index}{
    \code{data.frame} or \code{matrix} (or something acceptable to
    \code{\link{as.matrix}}) specifying indices of the basis
    functions.
  }
  \item{basis}{character string giving the name of the basis.}
  \item{window}{
    window (of class \code{"owin"}, \code{"box3"} or \code{"boxx"})
    giving the domain of the point process.
  }
  \item{reject_max}{
    integer giving the maximal number of trials for rejection sampling.
  }
  \item{progress}{
    integer giving the interval for making a progress report. The value
    zero turns reporting off.
  }
  \item{debug}{
    logical value indicating whether debug informationb
    should be outputted.
  }
  \item{\dots}{Ignored.}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\examples{
index <- expand.grid(-2:2,-2:2)
eig <- exp(-rowSums(index^2))
X <- rdpp(eig, index)
X
## To simulate a det. projection p. p. with the given indices set eig=1:
XX <- rdpp(1, index)
XX
}
\keyword{datagen}
\keyword{spatial}
\keyword{models}
