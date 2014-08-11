\name{nncorr}
\alias{nncorr}
\alias{nnmean}
\alias{nnvario}
\title{Nearest-Neighbour Correlation Indices of Marked Point Pattern}
\description{
  Computes nearest-neighbour correlation indices of a marked point
  pattern, including the nearest-neighbour mark product index
  (default case of \code{nncorr}),
  the nearest-neighbour mark index (\code{nnmean}),
  and the nearest-neighbour variogram index (\code{nnvario}).
}
\usage{
     nncorr(X,
            f = function(m1, m2) { m1 * m2 },
            \dots,
            use = "all.obs", method = c("pearson", "kendall", "spearman"),
            denominator=NULL)
     nnmean(X)
     nnvario(X)
}
\arguments{
  \item{X}{
    The observed point pattern.
    An object of class \code{"ppp"}.
  }
  \item{f}{
    Function \eqn{f} used in the definition of the
    nearest neighbour correlation. There is a sensible default
    that depends on the type of marks of \code{X}.
  }
  \item{\dots}{
    Extra arguments passed to \code{f}.
  }
  \item{use,method}{
    Arguments passed to the standard correlation function \code{\link{cor}}.
  }
  \item{denominator}{
    Internal use only.
  }
}
\details{
  The nearest neighbour correlation index \eqn{\bar n_f}{nbar}
  of a marked point process \eqn{X}
  is a number measuring the dependence between the mark of a typical point
  and the mark of its nearest neighbour. 

  The command \code{nncorr} computes the nearest neighbour correlation index
  based on any test function \code{f} provided by the user.
  The default behaviour of \code{nncorr} is to compute the
  nearest neighbour mark product index.
  The commands \code{nnmean} and \code{nnvario} are
  convenient abbreviations for other special choices of \code{f}.

  In the default case, \code{nncorr(X)} computes three different
  versions of the nearest-neighbour correlation index:
  the unnormalised, normalised, and classical correlations.
  \describe{
    \item{unnormalised:}{
      The \bold{unnormalised} nearest neighbour correlation (Stoyan and Stoyan,
      1994, section 14.7) is defined as
      \deqn{\bar n_f = E[f(M, M^\ast)]}{nbar[f] = E[f(M, M*)]}
      where \eqn{E[]} denotes mean value,
      \eqn{M} is the mark attached to a
      typical point of the point process, and \eqn{M^\ast}{M*} is the mark
      attached to its nearest neighbour (i.e. the nearest other point of the
      point process).
  
      Here \eqn{f} is any function
      \eqn{f(m_1,m_2)}{f(m1,m2)}
      with two arguments which are possible marks of the pattern,
      and which returns a nonnegative real value.
      Common choices of \eqn{f} are:
      for continuous real-valued marks,
      \deqn{f(m_1,m_2) = m_1 m_2}{f(m1,m2)= m1 * m2}
      for discrete marks (multitype point patterns),
      \deqn{f(m_1,m_2) = 1(m_1 = m_2)}{f(m1,m2)= (m1 == m2)}
      and for marks taking values in \eqn{[0,2\pi)}{[0,2 * pi)},
      \deqn{f(m_1,m_2) = \sin(m_1 - m_2)}{f(m1,m2) = sin(m1-m2).}
      For example, in the second case, the unnormalised nearest neighbour
      correlation \eqn{\bar n_f}{nbar[f]} equals the proportion of
      points in the pattern which have the same mark as their nearest
      neighbour.

      Note that \eqn{\bar n_f}{nbar[f]} is not a ``correlation''
      in the usual statistical sense. It can take values greater than 1.
    }
    \item{normalised:}{
      We can define a \bold{normalised} nearest neighbour correlation
      by 
      \deqn{\bar m_f = \frac{E[f(M,M^\ast)]}{E[f(M,M')]}}{mbar[f] = E[f(M,M*)]/E[f(M,M')]}
      where again \eqn{M} is the
      mark attached to a typical point, \eqn{M^\ast}{M*} is the mark
      attached to its nearest neighbour, and \eqn{M'} is an independent
      copy of \eqn{M} with the same distribution.
      This normalisation is also not a ``correlation''
      in the usual statistical sense, but is normalised so that 
      the value 1 suggests ``lack of correlation'':
      if the marks attached to the points of \code{X} are independent
      and identically distributed, then
      \eqn{\bar m_f = 1}{mbar[f] =  1}.
      The interpretation of values larger or smaller than 1 depends
      on the choice of function \eqn{f}.
    }
    \item{classical:}{
      Finally if the marks of \code{X} are real numbers,
      we can also compute the
      \bold{classical} correlation, that is, the correlation coefficient
      of the two random variables \eqn{M} and \eqn{M^\ast}{M*}.
      The classical correlation has a value between \eqn{-1} and \eqn{1}.
      Values close to \eqn{-1} or \eqn{1} indicate strong dependence between
      the marks.
    }
  }

  In the default case where \code{f} is not given,
  \code{nncorr(X)} computes
  \itemize{
    \item
    If the marks of \code{X} are real numbers, 
    the unnormalised and normalised
    versions of the nearest-neighbour product index
    \eqn{E[M \, M^\ast]}{E[M * M*]},
    and the classical correlation
    between \eqn{M} and \eqn{M^\ast}{M*}.
    \item
    If the marks of \code{X} are factor valued,
    the unnormalised and normalised
    versions of the nearest-neighbour equality index
    \eqn{P[M = M^\ast]}{P[M = M*]}.
  }

  The wrapper functions \code{nnmean} and \code{nnvario}
  compute the correlation indices for two special choices of the
  function \eqn{f(m_1,m_2)}{f(m1,m2)}.
  \itemize{
    \item
    \code{nnmean} computes the correlation indices for 
    \eqn{f(m_1,m_2) = m_1}{f(m1,m2) = m1}. The unnormalised index
    is simply the mean value of the mark of the neighbour of a typical point,
    \eqn{E[M^\ast]}{E[M*]}, while the normalised index is
    \eqn{E[M^\ast]/E[M]}{E[M*]/E[M]}, the ratio of the mean mark of the
    neighbour of a typical point to the mean mark of a typical point.
    \item 
    \code{nnvario} computes the correlation indices for 
    \eqn{f(m_1,m_2) = (1/2) (m_1-m_2)^2}{f(m1,m2) = (1/2) * (m1-m2)^2}.
  }

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) and must be a marked point pattern.
  (The marks may be a data frame, containing several columns of mark variables;
  each column is treated separately.)

  If the argument \code{f} is given, it
  must be a function, accepting two arguments \code{m1}
  and \code{m2} which are vectors of equal length containing mark
  values (of the same type as the marks of \code{X}).
  It must return a vector of numeric
  values of the same length as \code{m1} and \code{m2}.
  The values must be non-negative.

  The arguments \code{use} and \code{method} control
  the calculation of the classical correlation using \code{\link{cor}},
  as explained in the help file for \code{\link{cor}}.

  Other arguments may be passed to \code{f} through the \code{...}
  argument.
  
  This algorithm assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.
  Biases due to edge effects are
  treated using the \sQuote{border method} edge correction.
}
\value{
  Labelled vector of length 2 or 3
  containing the unnormalised and normalised
  nearest neighbour correlations, and the classical correlation
  if appropriate.
  Alternatively a matrix with 2 or 3 rows, containing this information
  for each mark variable.
}
\examples{
  data(finpines)
  nncorr(finpines)
  # heights of neighbouring trees are slightly negatively correlated

  data(amacrine)
  nncorr(amacrine)
  # neighbouring cells are usually of different type
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
