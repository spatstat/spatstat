\name{bind.fv}
\alias{bind.fv}
\alias{cbind.fv}
\title{
  Combine Function Value Tables
}
\description{
  Advanced Use Only.
  Combine objects of class \code{"fv"},
  or glue extra columns of data onto an existing \code{"fv"} object.
}
\usage{
\method{cbind}{fv}(...)
bind.fv(x, y, labl = NULL, desc = NULL, preferred = NULL)
}
\arguments{
  \item{\dots}{
    Any number of arguments, which are objects of class \code{"fv"}.
  }
  \item{x}{
    An object of class \code{"fv"}.
  }
  \item{y}{
    Either a data frame or an object of class \code{"fv"}.
  }
  \item{labl}{
    Plot labels (see \code{\link{fv}}) for columns of \code{y}.
    A character vector. 
  }
  \item{desc}{
    Descriptions (see \code{\link{fv}})
    for columns of \code{y}. A character vector.
  }
  \item{preferred}{
    Character string specifying the column which is to be the
    new recommended value of the function.
  }
}
\details{
  This documentation is provided
  for experienced programmers who want to modify the internal
  behaviour of \pkg{spatstat}.

  The function \code{cbind.fv} is a method for the generic
  \R function \code{\link{cbind}}. It combines any number of
  objects of class \code{"fv"} into a single object of
  class \code{"fv"}. The objects must be compatible, in the sense
  that they have identical values of the function argument.
  
  The function \code{bind.fv} is a lower level
  utility which glues additional columns onto an
  existing object \code{x} of class \code{"fv"}.
  It has two modes of use:
  \itemize{
    \item 
    If the additional dataset \code{y} is an object of class \code{"fv"}, then
    \code{x} and \code{y} must be compatible as described above.
    Then the columns of \code{y} that contain function values
    will be appended to the object \code{x}.
    \item 
    Alternatively if \code{y} is a data frame, then \code{y} must have the
    same number of rows as \code{x}. All columns of \code{y} will be
    appended to \code{x}.
  }
  The arguments \code{labl} and \code{desc} provide
  plot labels and description strings (as described in \code{\link{fv}})
  for the \emph{new} columns. If \code{y} is an object of class
  \code{"fv"} then \code{labl} and \code{desc} are optional, and
  default to the relevant entries in the object \code{y}. 
  If \code{y} is a data frame then
  \code{labl} and \code{desc} must be provided.
}
\value{
  An object of class \code{"fv"}.
}
\author{\adrian
  
  
  and \rolf
  
}
\examples{
   data(cells)
   K1 <- Kest(cells, correction="border")
   K2 <- Kest(cells, correction="iso")
   # remove column 'theo' to avoid duplication
   K2 <- K2[, names(K2) != "theo"]

   cbind(K1, K2)

   bind.fv(K1, K2, preferred="iso")

   # constrain border estimate to be monotonically increasing
   bm <- cumsum(c(0, pmax(0, diff(K1$border))))
   bind.fv(K1, data.frame(bmono=bm),
               "\%s[bmo](r)",
               "monotone border-corrected estimate of \%s",
               "bmono") 
}
\seealso{
  \code{\link{fv}}, 
  \code{\link{with.fv}}.

  \emph{Undocumented} functions for modifying an \code{"fv"} object
  include \code{fvnames}, \code{fvnames<-},
  \code{tweak.fv.entry} and \code{rebadge.fv}.
}
\keyword{spatial}
\keyword{attribute}

