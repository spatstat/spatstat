\name{Lcross.inhom}
\alias{Lcross.inhom}
\title{
  Inhomogeneous Cross Type L Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the cross-type \eqn{L} function.
}
\usage{
Lcross.inhom(X, i, j, \dots)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{L} function
    \eqn{L_{ij}(r)}{Lij(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be
    converted to a character string).
    Defaults to the second level of \code{marks(X)}.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{Kcross.inhom}}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{L_{ij}(r)}{Lij(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{L_{ij}(r)}{Lij(r)}
    for a marked Poisson process, identically equal to \code{r}
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L_{ij}(r)}{Lij(r)}
  obtained by the edge corrections named.
}
\details{
  This is a generalisation of the function \code{\link{Lcross}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Linhom}}.

  All the arguments are passed to \code{\link{Kcross.inhom}}, which 
  estimates the inhomogeneous multitype K function
  \eqn{K_{ij}(r)}{Kij(r)} for the point pattern.
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
  The arguments \code{i} and \code{j} are always interpreted as
  levels of the factor \code{X$marks}. They are converted to character
  strings if they are not already character strings.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.
}
\seealso{
 \code{\link{Lcross}},
 \code{\link{Linhom}},
 \code{\link{Kcross.inhom}}
}
\examples{
    # Lansing Woods data
    data(lansing)
    lansing <- lansing[seq(1,lansing$n, by=10)]
    ma <- split(lansing)$maple
    wh <- split(lansing)$whiteoak

    # method (1): estimate intensities by nonparametric smoothing
    lambdaM <- density.ppp(ma, sigma=0.15, at="points")
    lambdaW <- density.ppp(wh, sigma=0.15, at="points")
    L <- Lcross.inhom(lansing, "whiteoak", "maple", lambdaW, lambdaM)

    # method (2): fit parametric intensity model
    fit <- ppm(lansing, ~marks * polynom(x,y,2))
    # evaluate fitted intensities at data points
    # (these are the intensities of the sub-processes of each type)
    inten <- fitted(fit, dataonly=TRUE)
    # split according to types of points
    lambda <- split(inten, lansing$marks)
    L <- Lcross.inhom(lansing, "whiteoak", "maple",
              lambda$whiteoak, lambda$maple)
    
    # synthetic example: type A points have intensity 50,
    #                    type B points have intensity 100 * x
    lamB <- as.im(function(x,y){50 + 100 * x}, owin())
    X <- superimpose(A=runifpoispp(50), B=rpoispp(lamB))
    L <- Lcross.inhom(X, "A", "B",
        lambdaI=as.im(50, X$window), lambdaJ=lamB)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
