\name{Math.imlist}
\alias{Math.imlist}
\alias{Ops.imlist}
\alias{Complex.imlist}
\alias{Summary.imlist}
\title{S3 Group Generic methods for List of Images}
\description{
  These are group generic methods for the class \code{"imlist"}
  of lists of images. These methods 
  allows the usual mathematical functions and operators to be applied
  directly to lists of images. See Details for a list of implemented functions.
}
\usage{
## S3 methods for group generics have prototypes:
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm = TRUE)}
%NAMESPACE S3method("Math", "imlist")
%NAMESPACE S3method("Ops", "imlist")
%NAMESPACE S3method("Complex", "imlist")
%NAMESPACE S3method("Summary", "imlist")
}
\arguments{
  \item{x, z, e1, e2}{objects of class \code{"imlist"}.}
  \item{\dots}{further arguments passed to methods.}
  \item{na.rm}{logical: should missing values be removed?}
}
\details{
  Below is a list of mathematical functions and operators which are
  defined for lists of images. Not all functions will make sense for all
  types of images. For example, none of the functions in the \code{"Math"} group
  make sense for character-valued images. Note that the \code{"Ops"}
  group methods are implemented using \code{\link{eval.im}}, which tries
  to harmonise images via \code{\link{harmonise.im}} if they aren't
  compatible to begin with.

  \enumerate{
    \item Group \code{"Math"}:
    \itemize{
      \item
      \code{abs}, \code{sign}, \code{sqrt},\cr
      \code{floor}, \code{ceiling}, \code{trunc},\cr
      \code{round}, \code{signif}

      \item
      \code{exp}, \code{log},  \code{expm1}, \code{log1p},\cr
      \code{cos}, \code{sin}, \code{tan},\cr
      \code{cospi}, \code{sinpi}, \code{tanpi},\cr
      \code{acos}, \code{asin}, \code{atan}

      \code{cosh}, \code{sinh}, \code{tanh},\cr
      \code{acosh}, \code{asinh}, \code{atanh}

      \item
      \code{lgamma}, \code{gamma}, \code{digamma}, \code{trigamma}
      \item \code{cumsum}, \code{cumprod}, \code{cummax}, \code{cummin}
    }

    \item Group \code{"Ops"}:
    \itemize{
      \item
      \code{"+"}, \code{"-"}, \code{"*"}, \code{"/"},
      \code{"^"}, \code{"\%\%"}, \code{"\%/\%"}

      \item \code{"&"}, \code{"|"}, \code{"!"}

      \item \code{"=="}, \code{"!="},
      \code{"<"}, \code{"<="}, \code{">="}, \code{">"}
    }

    \item Group \code{"Summary"}:
    \itemize{
      \item \code{all}, \code{any}
      \item \code{sum}, \code{prod}
      \item \code{min}, \code{max}
      \item \code{range}
    }

    \item Group \code{"Complex"}:
    \itemize{
      \item \code{Arg}, \code{Conj}, \code{Im}, \code{Mod}, \code{Re}
    }
  }
}
\value{
  The result of \code{"Math"}, \code{"Ops"} and \code{"Complex"} group
  operations is another list of images.
  The result of \code{"Summary"} group operations is a numeric vector of
  length 1 or 2.
}
\seealso{
  \code{\link{Math.im}} or \code{\link{eval.im}}
  for evaluating expressions involving images.
}
\examples{
  a <- Smooth(finpines, 2)
  log(a)/2 - sqrt(a)
  range(a)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
