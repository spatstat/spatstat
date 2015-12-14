\name{Jest}
\alias{Jest}
\title{Estimate the J-function}
\description{
  Estimates the summary function \eqn{J(r)} for a point pattern in a 
  window of arbitrary shape.
}
\usage{
  Jest(X, ..., eps=NULL, r=NULL, breaks=NULL, correction=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{J(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{Ignored.}
  \item{eps}{
    the resolution of the discrete approximation to Euclidean distance
    (see below). There is a sensible default.
  }
  \item{r}{vector of values for the argument \eqn{r} at which \eqn{J(r)} 
    should be evaluated. There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \code{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{correction}{
    Optional. Character string specifying the choice of edge
    correction(s) in \code{\link{Fest}} and \code{\link{Gest}}.
    See Details.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{J} has been  estimated}
  \item{rs}{the ``reduced sample'' or ``border correction''
              estimator of \eqn{J(r)} computed from
	      the border-corrected estimates of \eqn{F} and \eqn{G} }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{J(r)} computed from
    the Kaplan-Meier estimates of \eqn{F} and \eqn{G} }
  \item{han}{the Hanisch-style estimator of \eqn{J(r)} computed from
    the Hanisch estimate of \eqn{G} and the Chiu-Stoyan estimate of
    \eqn{F}
  }
  \item{un}{the uncorrected estimate of \eqn{J(r)}
             computed from the uncorrected estimates of \eqn{F} and
	     \eqn{G}
  }
  \item{theo}{the theoretical value of \eqn{J(r)}
    for a stationary Poisson process: identically equal to \eqn{1}
  }
  The data frame also has \bold{attributes}
  \item{F}{
    the output of \code{\link{Fest}} for this point pattern,
    containing three estimates of the empty space function \eqn{F(r)}
    and an estimate of its hazard function
  }
  \item{G}{
    the output of \code{\link{Gest}} for this point pattern,
    containing three estimates of the nearest neighbour distance distribution
    function \eqn{G(r)} and an estimate of its hazard function
  }
}
\note{
  Sizeable amounts of memory may be needed during the calculation.
}
\details{
  The \eqn{J} function (Van Lieshout and Baddeley, 1996)
  of a stationary point process is defined as 
  \deqn{J(r) = \frac{1-G(r)}{1-F(r)} }{ %
    J(r) = (1-G(r))/(1-F(r))}
  where \eqn{G(r)} is the nearest neighbour distance distribution
  function of the point process (see \code{\link{Gest}}) 
  and \eqn{F(r)} is its empty space function (see \code{\link{Fest}}).

  For a completely random (uniform Poisson) point process,
  the \eqn{J}-function is identically equal to \eqn{1}. 
  Deviations \eqn{J(r) < 1} or \eqn{J(r) > 1}
  typically indicate spatial clustering or spatial regularity, respectively.
  The \eqn{J}-function is one of the few characteristics that can be
  computed explicitly for a wide range of point processes. 
  See Van Lieshout and Baddeley (1996), Baddeley et al (2000),
  Thonnes and Van Lieshout (1999)  for further information.

  An estimate of \eqn{J} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern. The estimate of \eqn{J(r)} is compared against the 
  constant function \eqn{1}.
  Deviations \eqn{J(r) < 1} or \eqn{J(r) > 1}
  may suggest spatial clustering or spatial regularity, respectively.

  This algorithm estimates the \eqn{J}-function
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window. 
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape. 

  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised by
  \code{\link{as.ppp}()}. 

  The functions \code{\link{Fest}} and \code{\link{Gest}} are called to 
  compute estimates of \eqn{F(r)} and \eqn{G(r)} respectively.
  These estimates are then combined by simply taking the ratio
  \eqn{J(r) = (1-G(r))/(1-F(r))}.

  In fact several different estimates are computed
  using different edge corrections (Baddeley, 1998).

  The Kaplan-Meier estimate (returned as \code{km}) is the ratio 
  \code{J = (1-G)/(1-F)} of the Kaplan-Meier estimates of
  \eqn{1-F} and \eqn{1-G} computed by
  \code{\link{Fest}} and \code{\link{Gest}} respectively.
  This is computed if \code{correction=NULL} or if \code{correction}
  includes \code{"km"}.

  The Hanisch-style estimate (returned as \code{han}) is the ratio 
  \code{J = (1-G)/(1-F)} where \code{F} is the Chiu-Stoyan estimate of
  \eqn{F} and \code{G} is the Hanisch estimate of \eqn{G}.
  This is computed if \code{correction=NULL} or if \code{correction}
  includes \code{"cs"} or \code{"han"}.

  The reduced-sample or border corrected estimate
  (returned as \code{rs}) is
  the same ratio \code{J = (1-G)/(1-F)}
  of the border corrected estimates. 
  This is computed if \code{correction=NULL} or if \code{correction}
  includes \code{"rs"} or \code{"border"}.

  These edge-corrected estimators are slightly biased for \eqn{J}, 
  since they are ratios
  of approximately unbiased estimators.
  The logarithm of the
  Kaplan-Meier estimate is exactly unbiased for \eqn{\log J}{log J}.

  The uncorrected estimate (returned as \code{un}
  and computed only if \code{correction} includes \code{"none"})
  is the ratio \code{J = (1-G)/(1-F)}
  of the uncorrected (``raw'') estimates of the survival functions
  of \eqn{F} and \eqn{G},
  which are the empirical distribution functions of the 
  empty space distances \code{Fest(X,\dots)$raw}
  and of the nearest neighbour distances 
  \code{Gest(X,\dots)$raw}. The uncorrected estimates
  of \eqn{F} and \eqn{G} are severely biased.
  However the uncorrected estimate of \eqn{J}
  is approximately unbiased (if the process is close to Poisson);
  it is insensitive to edge effects, and should be used when
  edge effects are severe (see Baddeley et al, 2000).
  
  The algorithm for \code{\link{Fest}}
  uses two discrete approximations which are controlled
  by the parameter \code{eps} and by the spacing of values of \code{r}
  respectively. See \code{\link{Fest}} for details.
  First-time users are strongly advised not to specify these arguments.

  Note that the value returned by \code{Jest} includes 
  the output of \code{\link{Fest}} and \code{\link{Gest}}
  as attributes (see the last example below).
  If the user is intending to compute the \code{F,G} and \code{J}
  functions for the point pattern, it is only necessary to
  call \code{Jest}.
}
\references{
  Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37--78.
  
  Baddeley, A.J. and Gill, R.D. 
    The empty space hazard of a spatial pattern.
    Research Report 1994/3, Department of Mathematics,
    University of Western Australia, May 1994.

  Baddeley, A.J. and Gill, R.D.
     Kaplan-Meier estimators of interpoint distance
		distributions for spatial point processes.
     \emph{Annals of Statistics} \bold{25} (1997) 263--292.

  Baddeley, A., Kerscher, M., Schladitz, K. and Scott, B.T.
  Estimating the \emph{J} function without edge correction.
  \emph{Statistica Neerlandica} \bold{54} (2000) 315--328.

  Borgefors, G.
     Distance transformations in digital images.
     \emph{Computer Vision, Graphics and Image Processing}
     \bold{34} (1986) 344--371.

  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.

  Thonnes, E. and Van Lieshout, M.N.M,
  A comparative study on the power of Van Lieshout and Baddeley's J-function.
  \emph{Biometrical Journal} \bold{41} (1999) 721--734.

  Van Lieshout, M.N.M. and Baddeley, A.J.
  A nonparametric measure of spatial interaction in point patterns.
  \emph{Statistica Neerlandica} \bold{50} (1996) 344--361.
}
\seealso{
  \code{\link{Jinhom}},
  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Kest}},
  \code{\link{km.rs}},
  \code{\link{reduced.sample}},
  \code{\link{kaplan.meier}}
}
\examples{
   data(cells)
   J <- Jest(cells, 0.01)
   plot(J, main="cells data")
   # values are far above J = 1, indicating regular pattern

   data(redwood)
   J <- Jest(redwood, 0.01, legendpos="center")
   plot(J, main="redwood data")
   # values are below J = 1, indicating clustered pattern
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}

