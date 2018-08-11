\name{edit.hyperframe}
\alias{edit.hyperframe}
\title{
  Invoke Text Editor on Hyperframe
}
\description{
  Invokes a text editor allowing the user to inspect and change
  entries in a hyperframe.
}
\usage{
\method{edit}{hyperframe}(name, \dots)
}
\arguments{
  \item{name}{
    A hyperframe (object of class \code{"hyperframe"}).
  }
  \item{\dots}{
    Other arguments passed to \code{\link[utils]{edit.data.frame}}.
  }
}
\details{
  The function \code{\link[utils]{edit}} is generic. This function
  is the methods for objects of class \code{"hyperframe"}.

  The hyperframe \code{name} is converted to a data frame or array,
  and the text editor is invoked. The user can change entries in the
  columns of data, and create new columns of data.

  Only the columns of atomic data
  (numbers, characters, factor values etc) can be edited.
  
  Note that the original object \code{name} is not changed;
  the function returns the edited dataset.
}
\value{
  Another hyperframe.
}
\author{
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link[utils]{edit.data.frame}},
  \code{\link{edit.ppp}}
}
\examples{
  if(interactive()) Z <- edit(flu)
}
\keyword{spatial}
\keyword{manip}

