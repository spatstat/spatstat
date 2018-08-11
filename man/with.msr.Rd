\name{with.msr}
\alias{with.msr}
\title{Evaluate Expression Involving Components of a Measure}
\description{
  An expression involving the names of components of a measure
  is evaluated.
}
\usage{
\method{with}{msr}(data, expr, \dots)
}
\arguments{
  \item{data}{
    A measure (object of class \code{"msr"}).
  }
  \item{expr}{
    An expression to be evaluated.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function
  \code{\link[base]{with}} for the class \code{"msr"}.
  The argument \code{data} should be an object of class \code{"msr"}
  representing a measure (a function which assigns a value to each
  subset of two-dimensional space).

  This function can be used to extract the components of the measure,
  or to perform more complicated manipulations of the components.

  The argument \code{expr} should be an un-evaluated expression
  in the \R language. The expression may involve any of the variable
  names listed below with their corresponding meanings. 
  \tabular{ll}{
    \code{qlocations} \tab (point pattern) all quadrature locations  \cr
    \code{qweights} \tab (numeric) all quadrature weights  \cr
    \code{density} \tab (numeric) density value at each quadrature point  \cr
    \code{discrete} \tab (numeric) discrete mass at each quadrature point  \cr
    \code{continuous} \tab (numeric) increment of continuous component  \cr
    \code{increment} \tab (numeric) increment of measure  \cr
    \code{is.atom} \tab (logical) whether quadrature point is an atom  \cr
    \code{atoms} \tab (point pattern) locations of atoms  \cr
    \code{atommass} \tab (numeric) massess of atoms
  }
  The measure is the sum of discrete and continuous components.
  The discrete component assigns non-zero mass to several points called atoms.
  The continuous component has a density which should be integrated
  over a region to determine the value for that region.

  An object of class \code{"msr"} approximates the continuous component
  by a sum over quadrature points. The quadrature points are chosen
  so that they include the atoms of the measure. In the list above,
  we have \code{increment = continuous + discrete},
  \code{continuous = density * qweights},
  \code{is.atom = (discrete > 0)},
  \code{atoms = qlocations[is.atom]} and
  \code{atommass = discrete[is.atom]}.
}
\value{
  The result of evaluating the expression could be
  an object of any kind.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{msr}},
  \code{\link{split.msr}},
  \code{\link{measureContinuous}}, \code{\link{measurePositive}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   rp <- residuals(fit, type="pearson")

   with(rp, atoms)
   with(rp, qlocations \%mark\% continuous)
}
\keyword{spatial}
\keyword{manip}

