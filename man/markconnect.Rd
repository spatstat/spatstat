\name{markconnect}
\alias{markconnect}
\title{
  Mark Connection Function
}
\description{
  Estimate the marked connection function
  of a multitype point pattern.
}
\usage{
markconnect(X, i, j, r=NULL,
         correction=c("isotropic", "Ripley", "translate"),
         method="density", \dots, normalise=FALSE)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. 
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
  }
  \item{j}{Number or character string identifying the type (mark value)
    of the points in \code{X} to which distances are measured.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the mark connection function \eqn{p_{ij}(r)}{p[ij](r)}
    should be evaluated. There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{markcorr}},
    or passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
  \item{normalise}{
    If \code{TRUE}, normalise the pair connection function by
    dividing it by \eqn{p_i p_j}{p[i]*p[j]}, the estimated probability
    that randomly-selected points will have marks \eqn{i} and \eqn{j}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the mark connection function \eqn{p_{ij}(r)}{p[i,j](r)}
    has been  estimated
  }
  \item{theo}{the theoretical value of \eqn{p_{ij}(r)}{p[i,j](r)}
    when the marks attached to different points are independent
  }
  together with a column or columns named 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{p_{ij}(r)}{p[i,j](r)}
  obtained by the edge corrections named.
}
\details{
  The mark connection function \eqn{p_{ij}(r)}{p[i,j](r)}
  of a multitype point process \eqn{X}
  is a measure of the dependence between the types of two 
  points of the process a distance \eqn{r} apart.

  Informally \eqn{p_{ij}(r)}{p[i,j](r)} is defined
  as the conditional probability,
  given that there is a point of the process at a location \eqn{u}
  and another point of the process at a location \eqn{v}
  separated by a distance \eqn{||u-v|| = r}, that the first point
  is of type \eqn{i} and the second point is of type \eqn{j}.
  See Stoyan and Stoyan (1994).

  If the marks attached to the points of \code{X} are independent
  and identically distributed, then
  \eqn{p_{ij}(r) \equiv p_i p_j}{p[i,j](r) =  p[i]p[j]} where
  \eqn{p_i}{p[i]} denotes the probability that a point is of type
  \eqn{i}. Values larger than this,
  \eqn{p_{ij}(r) > p_i p_j}{p[i,j](r) >  p[i]p[j]},
  indicate positive association between the two types,
  while smaller values indicate negative association.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a multitype point pattern (a marked point pattern
  with factor-valued marks).

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{p_{ij}(r)}{p[i,j](r)} is estimated.
  There is a sensible default.

  This algorithm assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.

  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kest}}.
  The edge corrections implemented here are
  \describe{
    \item{isotropic/Ripley}{Ripley's isotropic correction
      (see Ripley, 1988; Ohser, 1983).
      This is implemented only for rectangular and polygonal windows
      (not for binary masks) and is slow for complicated polygons.
    }
    \item{translate}{Translation correction (Ohser, 1983).
      Implemented for all window geometries.
    }
    \item{none}{No edge correction.}
  }
  The option \code{correction="none"} should only be used if the
  number of data points is extremely large (otherwise an edge correction
  is needed to correct bias).

  Note that the estimator assumes the process is stationary (spatially
  homogeneous). 

  The mark connection function is estimated using density estimation
  techniques. The user can choose between
  \describe{
    \item{\code{"density"}}{
      which uses the standard kernel
      density estimation routine \code{\link{density}}, and
      works only for evenly-spaced \code{r} values;
    }
    \item{\code{"loess"}}{
      which uses the function \code{loess} in the
      package \pkg{modreg};
    }
    \item{\code{"sm"}}{
      which uses the function \code{sm.density} in the
      package \pkg{sm} and is extremely slow;
    }
    \item{\code{"smrep"}}{
      which uses the function \code{sm.density} in the
      package \pkg{sm} and is relatively fast, but may require manual
      control of the smoothing parameter \code{hmult}.
    }
  }
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  Multitype pair correlation \code{\link{pcfcross}}
  and multitype K-functions \code{\link{Kcross}}, \code{\link{Kdot}}.

  Use \code{\link{alltypes}} to compute the mark connection functions
  between all pairs of types.

  Mark correlation \code{\link{markcorr}} and 
  mark variogram \code{\link{markvario}}
  for numeric-valued marks.

}
\examples{
    # Hughes' amacrine data
    # Cells marked as 'on'/'off'
    data(amacrine)
    M <- markconnect(amacrine, "on", "off")
    plot(M)

    # Compute for all pairs of types at once
    plot(alltypes(amacrine, markconnect))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}


