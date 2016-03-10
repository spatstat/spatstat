\name{as.hyperframe} %DontDeclareMethods
\Rdversion{1.1}
\alias{as.hyperframe}
\alias{as.hyperframe.default}
\alias{as.hyperframe.data.frame}
\alias{as.hyperframe.hyperframe}
\alias{as.hyperframe.listof}
\alias{as.hyperframe.anylist}
\title{
  Convert Data to Hyperframe
}
\description{
  Converts data from any suitable format into a hyperframe.
}
\usage{
as.hyperframe(x, \dots)

\method{as.hyperframe}{default}(x, \dots)

\method{as.hyperframe}{data.frame}(x, \dots, stringsAsFactors=FALSE)

\method{as.hyperframe}{hyperframe}(x, \dots)

\method{as.hyperframe}{listof}(x, \dots)

\method{as.hyperframe}{anylist}(x, \dots)
}
\arguments{
  \item{x}{
    Data in some other format.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{hyperframe}}.
  }
  \item{stringsAsFactors}{
    Logical. If \code{TRUE}, any column of the data frame \code{x} that
    contains character strings will be converted to a \code{factor}.
    If \code{FALSE}, no such conversion will occur.
  }
}
\details{
  A hyperframe is like a data frame, except that its entries
  can be objects of any kind.

  The generic function \code{as.hyperframe} converts any suitable kind
  of data into a hyperframe.

  There are methods for the classes
  \code{data.frame}, \code{listof}, \code{anylist} and a default method,
  all of which convert data that is like a hyperframe into
  a hyperframe object. (The method for the class \code{listof}
  and \code{anylist} converts a list of objects, of
  arbitrary type, into a hyperframe with one column.)
  These methods do not discard any information.

  There are also methods for other classes
  (see \code{\link{as.hyperframe.ppx}}) which extract
  the coordinates from a spatial dataset. These methods
  do discard some information.
}
\section{Conversion of Strings to Factors}{
  Note that \code{as.hyperframe.default}
  will convert a character vector to a factor.
  It behaves like \code{\link{as.data.frame}}.

  However \code{as.hyperframe.data.frame} does not convert strings to
  factors; it respects the structure of the data frame \code{x}.

  The behaviour can be changed using the argument \code{stringsAsFactors}.
}
\value{
  An object of class \code{"hyperframe"} created by
  \code{\link{hyperframe}}.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{hyperframe}},
  \code{\link{as.hyperframe.ppx}}
}
\examples{
   df <- data.frame(x=runif(4),y=letters[1:4])
   as.hyperframe(df)

   sims <- list()
   for(i in 1:3) sims[[i]] <- rpoispp(42)
   as.hyperframe(as.listof(sims))
   as.hyperframe(as.solist(sims))
}
\keyword{spatial}
\keyword{manip}
