\name{Kmulti}
\alias{Kmulti}
\title{
Marked K-Function
}
\description{
For a marked point pattern, 
estimate the multitype \eqn{K} function
which counts the expected number of points of subset \eqn{J}
within a given distance from a typical point in subset \code{I}.
}
\usage{
Kmulti(X, I, J, r=NULL, breaks=NULL, correction, \dots, ratio=FALSE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype \eqn{K} function
    \eqn{K_{IJ}(r)}{KIJ(r)} will be computed.
    It must be a marked point pattern.
    See under Details.
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured. See Details.
  }
  \item{J}{Subset index specifying the points in \code{X} to which
    distances are measured. See Details.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the multitype \eqn{K} function
    \eqn{K_{IJ}(r)}{KIJ(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{\dots}{Ignored.}
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\value{
An object of class \code{"fv"} (see \code{\link{fv.object}}).

Essentially a data frame containing numeric columns 
\item{r}{the values of the argument \eqn{r} 
at which the function \eqn{K_{IJ}(r)}{KIJ(r)} has been  estimated
}
\item{theo}{the theoretical value of  \eqn{K_{IJ}(r)}{KIJ(r)}
for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
}
together with a column or columns named 
\code{"border"}, \code{"bord.modif"},
\code{"iso"} and/or \code{"trans"},
according to the selected edge corrections. These columns contain
estimates of the function \eqn{K_{IJ}(r)}{KIJ(r)}
obtained by the edge corrections named.

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{K(r)}. 
}
\details{
The function \code{Kmulti}
generalises \code{\link{Kest}} (for unmarked point
patterns) and \code{\link{Kdot}} and \code{\link{Kcross}} (for
multitype point patterns) to arbitrary marked point patterns.

Suppose \eqn{X_I}{X[I]}, \eqn{X_J}{X[J]} are subsets, possibly
overlapping, of a marked point process.
The multitype \eqn{K} function 
is defined so that
\eqn{\lambda_J K_{IJ}(r)}{lambda[J] KIJ(r)} equals the expected number of
additional random points of \eqn{X_J}{X[J]} 
within a distance \eqn{r} of a
typical point of  \eqn{X_I}{X[I]}.
Here \eqn{\lambda_J}{lambda[J]}
is the intensity of \eqn{X_J}{X[J]} 
i.e. the expected number of points of \eqn{X_J}{X[J]} per unit area.
The function \eqn{K_{IJ}}{KIJ} is determined by the 
second order moment properties of \eqn{X}.

The argument \code{X} must be a point pattern (object of class
\code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.

The arguments \code{I} and \code{J} specify two subsets of the
point pattern. They may be any type of subset indices, for example,
logical vectors of length equal to \code{npoints(X)},
or integer vectors with entries in the range 1 to
\code{npoints(X)}, or negative integer vectors.

Alternatively, \code{I} and \code{J} may be \bold{functions}
that will be applied to the point pattern \code{X} to obtain
index vectors. If \code{I} is a function, then evaluating
\code{I(X)} should yield a valid subset index. This option
is useful when generating simulation envelopes using
\code{\link{envelope}}.

The argument \code{r} is the vector of values for the
distance \eqn{r} at which \eqn{K_{IJ}(r)}{KIJ(r)} should be evaluated. 
It is also used to determine the breakpoints
(in the sense of \code{\link{hist}})
for the computation of histograms of distances.

First-time users would be strongly advised not to specify \code{r}.
However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, 
and \code{max(r)} must be larger than the radius of the largest disc 
contained in the window. 

This algorithm assumes that \code{X} can be treated
as a realisation of a stationary (spatially homogeneous) 
random spatial point process in the plane, observed through
a bounded window.
The window (which is specified in \code{X} as \code{X$window})
may have arbitrary shape.

Biases due to edge effects are
treated in the same manner as in \code{\link{Kest}}.
The edge corrections implemented here are
\describe{
\item{border}{the border method or
``reduced sample'' estimator (see Ripley, 1988). This is
the least efficient (statistically) and the fastest to compute.
It can be computed for a window of arbitrary shape.
}
\item{isotropic/Ripley}{Ripley's isotropic correction
(see Ripley, 1988; Ohser, 1983).
This is currently implemented only for rectangular and polygonal windows.
}
\item{translate}{Translation correction (Ohser, 1983).
Implemented for all window geometries.
}
  }

  The pair correlation function \code{\link{pcf}} can also be applied to the
  result of \code{Kmulti}.

}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Diggle, P. J. (1986).
  Displaced amacrine cells in the retina of a
  rabbit : analysis of a bivariate spatial point pattern. 
  \emph{J. Neurosci. Meth.} \bold{18}, 115--125.
 
  Harkness, R.D and Isham, V. (1983)
  A bivariate spatial point pattern of ants' nests.
  \emph{Applied Statistics} \bold{32}, 293--303
 
  Lotwick, H. W. and Silverman, B. W. (1982).
  Methods for analysing spatial processes of several types of points.
  \emph{J. Royal Statist. Soc. Ser. B} \bold{44}, 406--413.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.

  Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
  Indices of dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

}
\section{Warnings}{
  The function \eqn{K_{IJ}}{KIJ} is not necessarily differentiable.

  The border correction (reduced sample) estimator of
  \eqn{K_{IJ}}{KIJ} used here is pointwise approximately 
  unbiased, but need not be a nondecreasing function of \eqn{r},
  while the true  \eqn{K_{IJ}}{KIJ} must be nondecreasing.
}
\seealso{
 \code{\link{Kcross}},
 \code{\link{Kdot}},
 \code{\link{Kest}},
 \code{\link{pcf}}
}
\examples{
     # Longleaf Pine data: marks represent diameter
    \testonly{
        longleaf <- longleaf[seq(1,npoints(longleaf), by=50), ]
    }
    K <- Kmulti(longleaf, longleaf$marks <= 15, longleaf$marks >= 25)
    plot(K)
    # functions determining subsets
    f1 <- function(X) { marks(X) <= 15 }
    f2 <- function(X) { marks(X) >= 15 }
    K <- Kmulti(longleaf, f1, f2)
   \testonly{
        rm(longleaf)
    }
 
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}



