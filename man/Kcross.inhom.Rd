\name{Kcross.inhom}
\alias{Kcross.inhom}
\title{
  Inhomogeneous Cross K Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the cross \eqn{K} function,
  which counts the expected number of points of type \eqn{j}
  within a given distance of a point of type \eqn{i},
  adjusted for spatially varying intensity.
}
\usage{
Kcross.inhom(X, i, j, lambdaI=NULL, lambdaJ=NULL, \dots,  r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate"),
         sigma=NULL, varcov=NULL,
         lambdaIJ=NULL,
         lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{K} function
    \eqn{K_{ij}(r)}{Kij(r)} will be computed.
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
  \item{lambdaI}{
    Optional.
    Values of the the estimated intensity of the sub-process of
    points of type \code{i}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the type \code{i} points in \code{X},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{lambdaJ}{
    Optional.
    Values of the the estimated intensity of the sub-process of
    points of type \code{j}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the type \code{j} points in \code{X},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location. 
  }
  \item{r}{
      Optional. Numeric vector giving the values of the argument \eqn{r}
      at which the cross K function
      \eqn{K_{ij}(r)}{Kij(r)} should be evaluated.
      There is a sensible default.
      First-time users are strongly advised not to specify this argument.
      See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
	This argument is for advanced use only.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"} ,\code{"translate"},
    \code{"translation"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{\dots}{
    Ignored.
  }
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel,
    used in computing leave-one-out kernel estimates of
    \code{lambdaI}, \code{lambdaJ} if they are omitted.
  }
  \item{varcov}{
    Variance-covariance matrix of anisotropic Gaussian kernel,
    used in computing leave-one-out kernel estimates of
    \code{lambdaI}, \code{lambdaJ} if they are omitted.
    Incompatible with \code{sigma}.
  }
  \item{lambdaIJ}{
    Optional. A matrix containing estimates of the
    product of the intensities \code{lambdaI} and \code{lambdaJ}
    for each pair of points of types \code{i} and \code{j} respectively.
  }
  \item{lambdaX}{
    Optional. Values of the intensity for all points of \code{X}.
    Either a pixel image (object of class \code{"im"}),
    a numeric vector containing the intensity values
    at each of the points in \code{X},
    a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"}),
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
    If present, this argument overrides both \code{lambdaI} and
    \code{lambdaJ}.
  }
  \item{update}{
    Logical value indicating what to do when
    \code{lambdaI}, \code{lambdaJ} or \code{lambdaX}
    is a fitted point process model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
    If \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.ppm}} or \code{\link{update.kppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{K_{ij}(r)}{Kij(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{K_{ij}(r)}{Kij(r)}
    for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K_{ij}(r)}{Kij(r)}
  obtained by the edge corrections named.
}
\details{
  This is a generalisation of the function \code{\link{Kcross}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Kinhom}}.

  The inhomogeneous cross-type \eqn{K} function is described by
  \ifelse{latex}{\out{M\o ller}}{Moller}
  and Waagepetersen (2003, pages 48-49 and 51-53).
  
  Briefly, given a multitype point process, suppose the sub-process
  of points of type \eqn{j} has intensity function
  \eqn{\lambda_j(u)}{lambda[j](u)} at spatial locations \eqn{u}.
  Suppose we place a mass of \eqn{1/\lambda_j(\zeta)}{1/lambda[j](z)}
  at each point \eqn{\zeta}{z} of type \eqn{j}. Then the expected total
  mass per unit area is 1. The
  inhomogeneous ``cross-type'' \eqn{K} function 
  \eqn{K_{ij}^{\mbox{inhom}}(r)}{K[ij]inhom(r)} equals the expected
  total mass within a radius \eqn{r} of a point of the process
  of type \eqn{i}. 
  
  If the process of type \eqn{i} points
  were independent of the process of type \eqn{j} points,
  then \eqn{K_{ij}^{\mbox{inhom}}(r)}{K[ij]inhom(r)}
  would equal \eqn{\pi r^2}{pi * r^2}.
  Deviations between the empirical \eqn{K_{ij}}{Kij} curve
  and the theoretical curve \eqn{\pi r^2}{pi * r^2} 
  suggest dependence between the points of types \eqn{i} and \eqn{j}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern, and the mark vector
  \code{X$marks} must be a factor.

  The arguments \code{i} and \code{j} will be interpreted as
  levels of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the number 3,
  \bold{not} the 3rd smallest level).
  If \code{i} and \code{j} are missing, they default to the first
  and second level of the marks factor, respectively.

  The argument \code{lambdaI} supplies the values
  of the intensity of the sub-process of points of type \code{i}.
  It may be either
  \describe{
    \item{a pixel image}{(object of class \code{"im"}) which
      gives the values of the type \code{i} intensity
      at all locations in the window containing \code{X};
    }
    \item{a numeric vector}{containing the values of the
      type \code{i} intensity evaluated only
      at the data points of type \code{i}. The length of this vector
      must equal the number of type \code{i} points in \code{X}.
    }
    \item{a function}{
      which can be evaluated to give values of the intensity at
      any locations.
    }
     \item{a fitted point process model}{
      (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
      whose fitted \emph{trend} can be used as the fitted intensity.
      (If \code{update=TRUE} the model will first be refitted to the
      data \code{X} before the trend is computed.)
    }
    \item{omitted:}{
      if \code{lambdaI} is omitted then it will be estimated
      using a leave-one-out kernel smoother. 
    }
  }
  If \code{lambdaI} is omitted, then it will be estimated using
  a `leave-one-out' kernel smoother,
  as described in Baddeley, \Moller 
  and Waagepetersen (2000).  The estimate of \code{lambdaI} for a given
  point is computed by removing the point from the
  point pattern, applying kernel smoothing to the remaining points using
  \code{\link{density.ppp}}, and evaluating the smoothed intensity
  at the point in question. The smoothing kernel bandwidth is controlled
  by the arguments \code{sigma} and \code{varcov}, which are passed to
  \code{\link{density.ppp}} along with any extra arguments.

  Similarly \code{lambdaJ} should contain
  estimated values of the intensity of the sub-process of points of
  type \code{j}. It may be either a pixel image, a function,
  a numeric vector, or omitted.
  
  Alternatively if the argument \code{lambdaX} is given, then it specifies
  the intensity values for all points of \code{X}, and the
  arguments \code{lambdaI}, \code{lambdaJ} will be ignored.
  
  The optional argument \code{lambdaIJ} is for advanced use only.
  It is a matrix containing estimated
  values of the products of these two intensities for each pair of
  data points of types \code{i} and \code{j} respectively.
  
  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{ij}(r)}{Kij(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must exceed the radius of the
  largest disc contained in the window.

  The argument \code{correction} chooses the edge correction
  as explained e.g. in \code{\link{Kest}}.

  The pair correlation function can also be applied to the
  result of \code{Kcross.inhom}; see \code{\link{pcf}}.
}
\references{
  Baddeley, A.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.

  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R.
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
 \code{\link{Kcross}},
 \code{\link{Kinhom}},
 \code{\link{Kdot.inhom}},
 \code{\link{Kmulti.inhom}},
 \code{\link{pcf}}
}
\examples{
    # Lansing Woods data
    woods <- lansing
    \testonly{woods <- woods[seq(1,npoints(woods), by=10)]}
    ma <- split(woods)$maple
    wh <- split(woods)$whiteoak

    # method (1): estimate intensities by nonparametric smoothing
    lambdaM <- density.ppp(ma, sigma=0.15, at="points")
    lambdaW <- density.ppp(wh, sigma=0.15, at="points")
    K <- Kcross.inhom(woods, "whiteoak", "maple", lambdaW, lambdaM)

    # method (2): leave-one-out
    K <- Kcross.inhom(woods, "whiteoak", "maple", sigma=0.15)

    # method (3): fit parametric intensity model
    fit <- ppm(woods ~marks * polynom(x,y,2))
    # alternative (a): use fitted model as 'lambda' argument
    K <- Kcross.inhom(woods, "whiteoak", "maple",
                      lambdaI=fit, lambdaJ=fit, update=FALSE)
    K <- Kcross.inhom(woods, "whiteoak", "maple",
                      lambdaX=fit, update=FALSE)
    # alternative (b): evaluate fitted intensities at data points
    # (these are the intensities of the sub-processes of each type)
    inten <- fitted(fit, dataonly=TRUE)
    # split according to types of points
    lambda <- split(inten, marks(woods))
    K <- Kcross.inhom(woods, "whiteoak", "maple",
              lambda$whiteoak, lambda$maple)
    
    # synthetic example: type A points have intensity 50,
    #                    type B points have intensity 100 * x
    lamB <- as.im(function(x,y){50 + 100 * x}, owin())
    X <- superimpose(A=runifpoispp(50), B=rpoispp(lamB))
    K <- Kcross.inhom(X, "A", "B",
        lambdaI=as.im(50, Window(X)), lambdaJ=lamB)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}
