\name{quasirandom}
\alias{quasirandom} %DoNotExport
\alias{vdCorput}
\alias{Halton}
\alias{Hammersley}
\title{
  Quasirandom Patterns
}
\description{
  Generates quasirandom sequences of numbers and
  quasirandom spatial patterns of points in any dimension.
}
\usage{
vdCorput(n, base)

Halton(n, bases = c(2, 3), raw = FALSE, simplify = TRUE)

Hammersley(n, bases = 2, raw = FALSE, simplify = TRUE)
}
\arguments{
  \item{n}{
    Number of points to generate.
  }
  \item{base}{
    A prime number giving the base of the sequence.
  }
  \item{bases}{
    Vector of prime numbers giving the bases of the sequences
    for each coordinate axis.
  }
  \item{raw}{
    Logical value indicating whether to return the coordinates
    as a matrix (\code{raw=TRUE}) or as a spatial point
    pattern (\code{raw=FALSE}, the default).
  }
  \item{simplify}{
    Argument passed to \code{\link{ppx}} indicating whether
    point patterns of dimension 2 or 3 should be
    returned as objects of class \code{"ppp"} or \code{"pp3"}
    respectively (\code{simplify=TRUE}, the default)
    or as objects of class \code{"ppx"} (\code{simplify=FALSE}).
  }
}
\details{
  The function \code{vdCorput} generates the quasirandom sequence
  of Van der Corput (1935) of length \code{n} with the given
  \code{base}. These are numbers between 0 and 1 which are in
  some sense uniformly distributed over the interval.

  The function \code{Halton} generates the Halton quasirandom sequence
  of points in \code{d}-dimensional space, where
  \code{d = length(bases)}. The values of the \eqn{i}-th coordinate
  of the points are generated using the van der Corput sequence with
  base equal to \code{bases[i]}.

  The function \code{Hammersley} generates the Hammersley set
  of points in \code{d+1}-dimensional space, where 
  \code{d = length(bases)}. The first \code{d} coordinates
  of the points are generated using the van der Corput sequence with
  base equal to \code{bases[i]}. The \code{d+1}-th coordinate
  is the sequence \code{1/n, 2/n, ..., 1}.

  If \code{raw=FALSE} (the default) then the Halton and Hammersley
  sets are interpreted as spatial point patterns of the
  appropriate dimension. They are returned as objects of
  class \code{"ppx"} (multidimensional point patterns)
  unless \code{simplify=TRUE} and \code{d=2} or \code{d=3}
  when they are returned as objects of class \code{"ppp"}
  or \code{"pp3"}.
  If \code{raw=TRUE}, the coordinates are returned as a matrix
  with \code{n} rows and \code{D} columns where \code{D} is the spatial
  dimension.
}
\value{
  For \code{vdCorput}, a numeric vector.

  For \code{Halton} and \code{Hammersley},
  an object of class \code{"ppp"}, \code{"pp3"} or
  \code{"ppx"}; or if \code{raw=TRUE}, a numeric matrix.
}
\references{
  Van der Corput, J. G. (1935) Verteilungsfunktionen.
  \emph{Proc. Ned. Akad. v. Wetensch.} \bold{38}: 813--821.

  Kuipers, L. and Niederreiter, H. (2005)
  \emph{Uniform distribution of sequences},
  Dover Publications.
}
\seealso{
  \code{\link{rQuasi}}
}
\examples{
   vdCorput(10, 2)

   plot(Halton(256, c(2,3)))

   plot(Hammersley(256, 3))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\keyword{spatial}
\keyword{datagen}
