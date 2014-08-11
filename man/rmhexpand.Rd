\name{rmhexpand}
\alias{rmhexpand}
\title{
  Specify Simulation Window or Expansion Rule
}
\description{
  Specify a spatial domain in which point process simulations
  will be performed. Alternatively, specify a rule which will be
  used to determine the simulation window.
}
\usage{
rmhexpand(x = NULL, ..., area = NULL, length = NULL, distance = NULL)
}
\arguments{
  \item{x}{
    Any kind of data determining the simulation window or the
    expansion rule.
    A window (object of class \code{"owin"}) specifying the
    simulation window, a numerical value specifying an expansion
    factor or expansion distance, a list containing
    one numerical value, an object of class \code{"rmhexpand"},
    or \code{NULL}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{area}{
    Area expansion factor. 
    Incompatible with other arguments.
  }
  \item{length}{
    Length expansion factor.
    Incompatible with other arguments.
    }
  \item{distance}{
    Expansion distance (buffer width).
    Incompatible with other arguments.
  }
}
\details{
  In the Metropolis-Hastings algorithm \code{\link{rmh}}
  for simulating spatial point processes, 
  simulations are usually carried out on a spatial domain that is
  larger than the original window of the point process model,
  then subsequently clipped to the original window.

  The command \code{rmhexpand} can be used to specify the simulation window,
  or to specify a rule which will later be used to determine the
  simulation window from data.

  The arguments are all incompatible: at most one of them should be
  given.

  If the first argument \code{x} is given, it may be any of the
  following:
  \itemize{
    \item
    a window (object of class \code{"owin"}) specifying the
    simulation window.
    \item
    an object of class \code{"rmhexpand"} specifying the
    expansion rule.
    \item
    a single numerical value, without attributes.
    This will be interpreted as the value of the argument \code{area}.
    \item
    either \code{c(area=v)} or \code{list(area=v)},
    where \code{v} is a single numeric value.
    This will be interpreted as the value of the argument \code{area}.
    \item
    either \code{c(length=v)} or \code{list(length=v)},
    where \code{v} is a single numeric value.
    This will be interpreted as the value of the argument \code{length}.
    \item
    either \code{c(distance=v)} or \code{list(distance=v)},
    where \code{v} is a single numeric value.
    This will be interpreted as the value of the argument \code{distance}.
    \item 
    \code{NULL}, meaning that the expansion rule is not yet
    determined. 
  }

  If one of the arguments \code{area}, \code{length} or \code{distance}
  is given, then the simulation window is determined from the original
  data window as follows.
  \describe{
    \item{area}{
      The bounding box of the original data window will be extracted,
      and the simulation window will be a scalar dilation of this rectangle.
      The argument \code{area} should be a numerical value, greater than
      or equal to 1. It specifies the area expansion factor, i.e. the ratio
      of the area of the simulation window to the area of the
      original point process window's bounding box.
    }
    \item{length}{
      The bounding box of the original data window will be extracted,
      and the simulation window will be a scalar dilation of this rectangle.
      The argument \code{length} should be a numerical value, greater than
      or equal to 1. It specifies the length expansion factor, i.e. the ratio
      of the width (height) of the simulation window to the width
      (height) of the original point process window's bounding box.
    }
    \item{distance}{
      The argument \code{distance} should be a numerical value, greater than
      or equal to 0. It specifies the width of a buffer region around the
      original data window. 
      If the original data window is a rectangle, then this window
      is extended by a margin of width equal to \code{distance}
      around all sides of the original rectangle. The result is a
      rectangle.
      If the original data window is not a rectangle, then
      morphological dilation is applied using
      \code{\link{dilation.owin}}
      so that a margin or buffer of width equal to \code{distance}
      is created around all sides of the original window. The result
      is a non-rectangular window, typically of a different shape.
    }
  }
}
\section{Undetermined expansion}{
  If \code{expand=NULL}, this is interpreted to mean that the
  expansion rule is \dQuote{not yet decided}. Expansion will be decided
  later, by the simulation algorithm \code{\link{rmh}}.
  If the model cannot be expanded (for example if the covariate data
  in the model are not available on a larger domain) then expansion
  will not occur. If the model can be expanded, then
  if the point process model has a finite interaction range \code{r},
  the default is \code{rmhexpand(distance=2*r)}, and
  otherwise \code{rmhexpand(area=2)}.
}
\value{
  An object of class \code{"rmhexpand"} specifying the
  expansion rule. There is a \code{print} method for this class.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{expand.owin}} to apply the rule to a window.

  \code{\link{will.expand}} to test whether expansion will occur.

  \code{\link{rmh}}, \code{\link{rmhcontrol}} for background details.
}
\examples{
  rmhexpand()
  rmhexpand(2)
  rmhexpand(1)
  rmhexpand(length=1.5)
  rmhexpand(distance=0.1)
  rmhexpand(letterR)
}
\keyword{spatial}
\keyword{datagen}
