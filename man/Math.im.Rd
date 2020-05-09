\name{Math.im}
\alias{Math.im}
\alias{Ops.im}
\alias{Complex.im}
\alias{Summary.im}
\title{S3 Group Generic methods for images}
\description{
  These are group generic methods for images of class \code{"im"}, which
  allows for usual mathematical functions and operators to be applied
  directly to images. See Details for a list of implemented functions.
}
\usage{
## S3 methods for group generics have prototypes:
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm=FALSE, drop=TRUE)}
%NAMESPACE S3method("Math", "im")
%NAMESPACE S3method("Ops", "im")
%NAMESPACE S3method("Complex", "im")
%NAMESPACE S3method("Summary", "im")
}
\arguments{
  \item{x, z, e1, e2}{objects of class \code{"im"}.}
  \item{\dots}{further arguments passed to methods.}
  \item{na.rm,drop}{
    Logical values specifying whether missing values
    should be removed. This will happen if either
    \code{na.rm=TRUE} or \code{drop=TRUE}. See Details.
  }
}

\details{
  Below is a list of mathematical functions and operators which are
  defined for images. Not all functions will make sense for all
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

  For the \code{Summary} group, the generic has an argument
  \code{na.rm=FALSE}, but for pixel images it makes
  sense to set \code{na.rm=TRUE} so that pixels outside the domain of
  the image are ignored. To enable this, we added the argument
  \code{drop}. Pixel values that are \code{NA} are removed
  if \code{drop=TRUE} or if \code{na.rm=TRUE}.

  For the \code{Ops} group, one of the arguments is permitted to be
  a single atomic value instead of an image.
}
\seealso{
  \code{\link{eval.im}} for evaluating expressions involving images.
}
\examples{
  ## Convert gradient values to angle of inclination:
  V <- atan(bei.extra$grad) * 180/pi
  ## Make logical image which is TRUE when heat equals 'Moderate': 
  A <- (gorillas.extra$heat == "Moderate")
  ## Summary:
  any(A)
  ## Complex:
  Z <- exp(1 + V * 1i)
  Z
  Re(Z)
}
\author{
  \spatstatAuthors and Kassel Hingee.
}
\keyword{spatial}
\keyword{methods}
