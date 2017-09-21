\name{Math.linim}
\alias{Math.linim}
\alias{Ops.linim}
\alias{Summary.linim}
\alias{Complex.linim}
\title{S3 Group Generic Methods for Images on a Linear Network}
\description{
  These are group generic methods for images of class \code{"linim"}, which
  allows for usual mathematical functions and operators to be applied
  directly to pixel images on a linear network.
  See Details for a list of implemented functions.
}
\usage{
## S3 methods for group generics have prototypes:
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm = FALSE)}
%NAMESPACE S3method("Math", "linim")
%NAMESPACE S3method("Ops", "linim")
%NAMESPACE S3method("Complex", "linim")
%NAMESPACE S3method("Summary", "linim")
}
\arguments{
  \item{x, z, e1, e2}{objects of class \code{"linim"}.}
  \item{\dots}{further arguments passed to methods.}
  \item{na.rm}{logical: should missing values be removed?}
}
\details{
  An object of class \code{"linim"} represents a pixel image on
  a linear network. See \code{\link{linim}}.
  
  Below is a list of mathematical functions and operators which are
  defined for these images. Not all functions will make sense for all
  types of images. For example, none of the functions in the \code{"Math"} group
  make sense for character-valued images. Note that the \code{"Ops"}
  group methods are implemented using \code{\link{eval.linim}}.

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
\seealso{
  \code{\link{eval.linim}} for evaluating expressions involving images.
}
\examples{
  fx <- function(x,y,seg,tp) { (x - y)^2 }
  fL <- linfun(fx, simplenet)
  Z <- as.linim(fL)
  A <- Z+2
  A <- -Z
  A <- sqrt(Z)
  A <- !(Z > 0.1)
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{methods}
