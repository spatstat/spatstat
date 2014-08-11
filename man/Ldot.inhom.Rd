\name{Ldot.inhom}
\alias{Ldot.inhom}
\title{
  Inhomogeneous Multitype L Dot Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the dot \eqn{L} function.
}
\usage{
Ldot.inhom(X, i, \dots)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{L} function
    \eqn{L_{i\bullet}(r)}{Li.(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{Kdot.inhom}}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{L_{i\bullet}(r)}{Li.(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{L_{i\bullet}(r)}{Li.(r)}
    for a marked Poisson process, identical to \eqn{r}.
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L_{i\bullet}(r)}{Li.(r)}
  obtained by the edge corrections named.
}
\details{
  This a generalisation of the function \code{\link{Ldot}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Linhom}}.

  All the arguments are passed to \code{\link{Kdot.inhom}}, which 
  estimates the inhomogeneous multitype K function
  \eqn{K_{i\bullet}(r)}{Ki.(r)} for the point pattern.
  The resulting values are then
  transformed by taking \eqn{L(r) = \sqrt{K(r)/\pi}}{L(r) = sqrt(K(r)/pi)}.
}
\references{
  Moller, J. and Waagepetersen, R.
  Statistical Inference and Simulation for Spatial Point Processes
  Chapman and Hall/CRC
  Boca Raton, 2003.
}
\section{Warnings}{
  The argument \code{i} is interpreted as
  a level of the factor \code{X$marks}. It is converted to a character
  string if it is not already a character string.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.
}
\seealso{
 \code{\link{Ldot}},
 \code{\link{Linhom}},
 \code{\link{Kdot.inhom}},
 \code{\link{Lcross.inhom}}.
}
\examples{
    # Lansing Woods data
    data(lansing)
    lansing <- lansing[seq(1,lansing$n, by=10)]
    ma <- split(lansing)$maple
    lg <- unmark(lansing)

    # Estimate intensities by nonparametric smoothing
    lambdaM <- density.ppp(ma, sigma=0.15, at="points")
    lambdadot <- density.ppp(lg, sigma=0.15, at="points")
    L <- Ldot.inhom(lansing, "maple", lambdaI=lambdaM,
                                      lambdadot=lambdadot)


    # synthetic example: type A points have intensity 50,
    #                    type B points have intensity 50 + 100 * x
    lamB <- as.im(function(x,y){50 + 100 * x}, owin())
    lamdot <- as.im(function(x,y) { 100 + 100 * x}, owin())
    X <- superimpose(A=runifpoispp(50), B=rpoispp(lamB))
    L <- Ldot.inhom(X, "B",  lambdaI=lamB,     lambdadot=lamdot)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

