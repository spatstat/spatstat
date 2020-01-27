\name{markcorr}
\alias{markcorr}
\title{
  Mark Correlation Function
}
\description{
  Estimate the marked correlation function
  of a marked point pattern.
}
\usage{
markcorr(X, f = function(m1, m2) { m1 * m2}, r=NULL,
         correction=c("isotropic", "Ripley", "translate"),
         method="density", \dots, weights=NULL,
         f1=NULL, normalise=TRUE, fargs=NULL, internal=NULL)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. 
  }
  \item{f}{Optional. Test function \eqn{f} used in the definition of the
    mark correlation function. An \R function with at least two
    arguments. There is a sensible default.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the mark correlation function 
    \eqn{k_f(r)}{k[f](r)} should be evaluated.
    There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"}, \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  
  \item{\dots}{
    Arguments passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
  \item{weights}{
    Optional. Numeric weights for each data point in \code{X}.
    A numeric vector, a pixel image, or a \code{function(x,y)}.
    Alternatively, an \code{expression} to be evaluated to yield the
    weights; the expression may involve the variables
    \code{x,y,marks} representing the coordinates and marks of\code{X}. 
  }
  \item{f1}{
    An alternative to \code{f}. If this argument is given,
    then \eqn{f} is assumed to take the form
    \eqn{f(u,v)=f_1(u)f_1(v)}{f(u,v)=f1(u) * f1(v)}.
  }
  \item{normalise}{
    If \code{normalise=FALSE},
    compute only the numerator of the expression for the
    mark correlation.
  }
  \item{fargs}{
    Optional. A list of extra arguments to be passed to the function
    \code{f} or \code{f1}.
  }
  \item{internal}{Do not use this argument.}
}
\value{
  A function value table (object of class \code{"fv"})
  or a list of function value tables, one for each column of marks.
  
  An object of class \code{"fv"} (see \code{\link{fv.object}})
  is essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the mark correlation function \eqn{k_f(r)}{k[f](r)}
    has been  estimated
  }
  \item{theo}{the theoretical value of \eqn{k_f(r)}{k[f](r)}
    when the marks attached to different points are independent,
    namely 1
  }
  together with a column or columns named 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the mark correlation function \eqn{k_f(r)}{k[f](r)}
  obtained by the edge corrections named.
}
\details{
  By default, this command calculates an estimate of
  Stoyan's mark correlation \eqn{k_{mm}(r)}{k[mm](r)}
  for the point pattern.

  Alternatively if the argument \code{f} or \code{f1} is given, then it
  calculates Stoyan's generalised mark correlation \eqn{k_f(r)}{k[f](r)}
  with test function \eqn{f}.

  Theoretical definitions are as follows
  (see Stoyan and Stoyan (1994, p. 262)):
  \itemize{
    \item
    For a point process \eqn{X} with numeric marks,
    Stoyan's mark correlation function \eqn{k_{mm}(r)}{k[mm](r)},
    is
    \deqn{
      k_{mm}(r) = \frac{E_{0u}[M(0) M(u)]}{E[M,M']}
    }{
      k[mm](r) = E[0u](M(0) * M(u))/E(M * M')
    }
    where \eqn{E_{0u}}{E[0u]} denotes the conditional expectation
    given that there are points of the process at the locations
    \eqn{0} and \eqn{u} separated by a distance \eqn{r},
    and where \eqn{M(0),M(u)} denote the marks attached to these
    two points. On the denominator, \eqn{M,M'} are random marks
    drawn independently from the marginal distribution of marks,
    and \eqn{E} is the usual expectation.
    \item
    For a multitype point process \eqn{X}, the mark correlation is 
    \deqn{
      k_{mm}(r) = \frac{P_{0u}[M(0) M(u)]}{P[M = M']}
    }{
      k[mm](r) = P[0u](M(0) = M(u))/P(M = M')
    }
    where \eqn{P} and \eqn{P_{0u}}{P[0u]} denote the
    probability and conditional probability.
    \item 
    The \emph{generalised} mark correlation function \eqn{k_f(r)}{k[f](r)}
    of a marked point process \eqn{X}, with test function \eqn{f},
    is
    \deqn{
      k_f(r) = \frac{E_{0u}[f(M(0),M(u))]}{E[f(M,M')]}
    }{
      k[f](r) = E[0u](f(M(0),M(u))]/E(f(M,M'))
    }
  }

  The test function \eqn{f} is any function
  \eqn{f(m_1,m_2)}{f(m1,m2)}
  with two arguments which are possible marks of the pattern,
  and which returns a nonnegative real value.
  Common choices of \eqn{f} are:
  for continuous nonnegative real-valued marks,
  \deqn{f(m_1,m_2) = m_1 m_2}{f(m1,m2)= m1 * m2}
  for discrete marks (multitype point patterns),
  \deqn{f(m_1,m_2) = 1(m_1 = m_2)}{f(m1,m2)= (m1 == m2)}
  and for marks taking values in \eqn{[0,2\pi)}{[0,2 * pi)},
  \deqn{f(m_1,m_2) = \sin(m_1 - m_2)}{f(m1,m2) = sin(m1-m2)}.
  
  Note that \eqn{k_f(r)}{k[f](r)} is not a ``correlation''
  in the usual statistical sense. It can take any 
  nonnegative real value. The value 1 suggests ``lack of correlation'':
  if the marks attached to the points of \code{X} are independent
  and identically distributed, then
  \eqn{k_f(r) \equiv 1}{k[f](r) =  1}.
  The interpretation of values larger or smaller than 1 depends
  on the choice of function \eqn{f}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern.

  The argument \code{f} determines the function to be applied to
  pairs of marks. It has a sensible default, which depends on the
  kind of marks in \code{X}. If the marks
  are numeric values, then \code{f <- function(m1, m2) { m1 * m2}}
  computes the product of two marks.
  If the marks are a factor (i.e. if \code{X} is a multitype point
  pattern) then \code{f <- function(m1, m2) { m1 == m2}} yields
  the value 1 when the two marks are equal, and 0 when they are unequal.
  These are the conventional definitions for numerical
  marks and multitype points respectively.

  The argument \code{f} may be specified by the user.
  It must be an \R function, accepting two arguments \code{m1}
  and \code{m2} which are vectors of equal length containing mark
  values (of the same type as the marks of \code{X}).
  (It may also take additional arguments, passed through \code{fargs}).
  It must return a vector of numeric
  values of the same length as \code{m1} and \code{m2}.
  The values must be non-negative, and \code{NA} values are not permitted.

  Alternatively the user may specify the argument \code{f1}
  instead of \code{f}. This indicates that the test function \eqn{f}
  should take the form \eqn{f(u,v)=f_1(u)f_1(v)}{f(u,v)=f1(u) * f1(v)}
  where \eqn{f_1(u)}{f1(u)} is given by the argument \code{f1}.
  The argument \code{f1} should be an \R function with at least one
  argument.
  (It may also take additional arguments, passed through \code{fargs}).
  
  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{k_f(r)}{k[f](r)} is estimated.

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
      (not for binary masks).
    }
    \item{translate}{Translation correction (Ohser, 1983).
      Implemented for all window geometries, but slow for
      complex windows. 
    }
  }
  Note that the estimator assumes the process is stationary (spatially
  homogeneous). 

  The numerator and denominator of the mark correlation function
  (in the expression above) are estimated using density estimation
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
  If \code{normalise=FALSE} then the algorithm will compute
  only the numerator
  \deqn{
    c_f(r) = E_{0u} f(M(0),M(u))
  }{
    c[f](r) = E[0u] f(M(0),M(u))
  }
  of the expression for the mark correlation function.
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  Mark variogram \code{\link{markvario}} for numeric marks.
  
  Mark connection function \code{\link{markconnect}} and 
  multitype K-functions \code{\link{Kcross}}, \code{\link{Kdot}}
  for factor-valued marks.

  Mark cross-correlation function \code{\link{markcrosscorr}}
  for point patterns with several columns of marks.
  
  \code{\link{Kmark}} to estimate a cumulative function
  related to the mark correlation function.
}
\examples{
    # CONTINUOUS-VALUED MARKS:
    # (1) Spruces
    # marks represent tree diameter
    # mark correlation function
    ms <- markcorr(spruces)
    plot(ms)

    # (2) simulated data with independent marks
    X <- rpoispp(100)
    X <- X \%mark\% runif(npoints(X))
    \dontrun{
    Xc <- markcorr(X)
    plot(Xc)
    }
    
    # MULTITYPE DATA:
    # Hughes' amacrine data
    # Cells marked as 'on'/'off'
    # (3) Kernel density estimate with Epanecnikov kernel
    # (as proposed by Stoyan & Stoyan)
    M <- markcorr(amacrine, function(m1,m2) {m1==m2},
                  correction="translate", method="density",
                  kernel="epanechnikov")
    plot(M)
    # Note: kernel="epanechnikov" comes from help(density)

    # (4) Same again with explicit control over bandwidth
    \dontrun{
    M <- markcorr(amacrine, 
                  correction="translate", method="density",
                  kernel="epanechnikov", bw=0.02)
    # see help(density) for correct interpretation of 'bw'
    }

   \testonly{
    niets <- markcorr(amacrine, function(m1,m2){m1 == m2}, method="loess")
    if(require(sm))
      niets <- markcorr(X, correction="isotropic", method="smrep", hmult=2)
    }

   # weighted mark correlation
   Y <- subset(betacells, select=type)
   a <- marks(betacells)$area
   v <- markcorr(Y, weights=a)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}


