\name{Iest}
\alias{Iest}
\title{Estimate the I-function}
\description{
  Estimates the summary function \eqn{I(r)} for a multitype point pattern.
}
\usage{
  Iest(X, ..., eps=NULL, r=NULL, breaks=NULL, correction=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{I(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{Ignored.}
  \item{eps}{
    the resolution of the discrete approximation to Euclidean distance
    (see below). There is a sensible default.
  }
  \item{r}{Optional. Numeric vector of values for the argument \eqn{r}
  	at which \eqn{I(r)} 
    should be evaluated. There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \code{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{correction}{
    Optional. Vector of character strings specifying the edge correction(s)
    to be used by \code{\link{Jest}}.
  }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{I} has been  estimated}
  \item{rs}{the ``reduced sample'' or ``border correction''
              estimator of \eqn{I(r)} computed from
	      the border-corrected estimates of \eqn{J} functions}
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{I(r)} computed from
    the Kaplan-Meier estimates of \eqn{J} functions}
  \item{han}{the Hanisch-style estimator of \eqn{I(r)} computed from
    the Hanisch-style estimates of \eqn{J} functions}
  \item{un}{the uncorrected estimate of \eqn{I(r)}
             computed from the uncorrected estimates of \eqn{J}
  }
  \item{theo}{the theoretical value of \eqn{I(r)}
    for a stationary Poisson process: identically equal to \eqn{0}
  }
}
\note{
  Sizeable amounts of memory may be needed during the calculation.
}
\details{
  The \eqn{I} function 
  summarises the dependence between types in a multitype point process
  (Van Lieshout and Baddeley, 1999)
  It is based on the concept of the \eqn{J} function for an
  unmarked point process (Van Lieshout and Baddeley, 1996).
  See \code{\link{Jest}} for information about the \eqn{J} function.
  
  The \eqn{I} function is defined as 
  \deqn{ %
    I(r) = \sum_{i=1}^m p_i J_{ii}(r) %
    - J_{\bullet\bullet}(r)}{ %
    I(r) = (sum p[i] Jii(r)) - J(r)
  }
  where \eqn{J_{\bullet\bullet}}{J} is the \eqn{J} function for
  the entire point process ignoring the marks, while 
  \eqn{J_{ii}}{Jii} is the \eqn{J} function for the
  process consisting of points of type \eqn{i} only,
  and \eqn{p_i}{p[i]} is the proportion of points which are of type \eqn{i}.

  The \eqn{I} function is designed to measure dependence between
  points of different types, even if the points are
  not Poisson. Let \eqn{X} be a stationary multitype point process,
  and write \eqn{X_i}{X[i]} for the process of points of type \eqn{i}.
  If the processes \eqn{X_i}{X[i]} are independent of each other,
  then the \eqn{I}-function is identically equal to \eqn{0}. 
  Deviations \eqn{I(r) < 1} or \eqn{I(r) > 1}
  typically indicate negative and positive association, respectively,
  between types.
  See Van Lieshout and Baddeley (1999)
  for further information.

  An estimate of \eqn{I} derived from a multitype spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern. The estimate of \eqn{I(r)} is compared against the 
  constant function \eqn{0}.
  Deviations \eqn{I(r) < 1} or \eqn{I(r) > 1}
  may suggest negative and positive association, respectively.

  This algorithm estimates the \eqn{I}-function
  from the multitype point pattern \code{X}.
  It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial marked point process in the plane, observed through
  a bounded window. 

  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised by
  \code{\link{as.ppp}()}. It must be a multitype point pattern
  (it must have a \code{marks} vector which is a \code{factor}).

  The function \code{\link{Jest}} is called to 
  compute estimates of the \eqn{J} functions in the formula above.
  In fact three different estimates are computed
  using different edge corrections. See \code{\link{Jest}} for
  information.
}
\references{
  Van Lieshout, M.N.M. and Baddeley, A.J. (1996)
  A nonparametric measure of spatial interaction in point patterns.
  \emph{Statistica Neerlandica} \bold{50}, 344--361.

  Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
  Indices of dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

}
\seealso{
  \code{\link{Jest}}
}
\examples{
   data(amacrine)
   Ic <- Iest(amacrine)
   plot(Ic, main="Amacrine Cells data")
   # values are below I= 0, suggesting negative association
   # between 'on' and 'off' cells.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}

