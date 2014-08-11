\name{rmhstart}
\alias{rmhstart}
\alias{rmhstart.default}
\title{Determine Initial State for Metropolis-Hastings Simulation.}
\description{
  Builds a description of the initial state
  for the Metropolis-Hastings
  algorithm. 
}
\usage{
   rmhstart(start, \dots)
   \method{rmhstart}{default}(start=NULL, \dots, n.start=NULL, x.start=NULL)
}
\arguments{
  \item{start}{An existing description of the initial state in some
    format. Incompatible with the arguments listed below.
  }
  \item{\dots}{There should be no other arguments.}
  \item{n.start}{
    Number of initial points (to be randomly generated).
    Incompatible with \code{x.start}.
  }
  \item{x.start}{
    Initial point pattern configuration.
    Incompatible with \code{n.start}.
  }
}
\value{
  An object of class \code{"rmhstart"}, which is essentially
  a list of parameters describing the initial point pattern
  and (optionally) the initial state of the random number generator.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the initial state.
}
\details{
  Simulated realisations of many point process models
  can be generated using the Metropolis-Hastings algorithm
  implemented in \code{\link{rmh}}.
  
  This function \code{rmhstart}
  creates a full description of the initial state of the
  Metropolis-Hastings algorithm,
  \emph{including possibly the initial state of the random number generator},
  for use in a subsequent call to \code{\link{rmh}}. It also
  checks that the initial state is valid.

  The initial state should be specified \bold{either} by the
  first argument \code{start} \bold{or} by the other arguments
  \code{n.start}, \code{x.start} etc.
  
  If \code{start} is a list, then it should have components named 
  \code{n.start} or \code{x.start},
  with the same interpretation as described below.

  The arguments are:
    \describe{
      \item{n.start}{
	The number of \dQuote{initial} points to be randomly
	(uniformly) generated in the simulation window \code{w}.
	Incompatible with \code{x.start}.

	For a multitype point process, \code{n.start} may be a vector
	(of length equal to the number of types) giving the number
	of points of each type to be generated.  

	If expansion of the simulation window is selected (see the argument
	\code{expand} to \code{\link{rmhcontrol}}),
	then the actual number of starting points in the simulation
	will be \code{n.start} multiplied by the expansion factor
	(ratio of the areas of the expanded window and original window).

	For faster convergence of the Metropolis-Hastings algorithm,
	the value of \code{n.start} should be roughly equal to
	(an educated guess at) the expected number of points
	for the point process inside the window.
      }
      \item{x.start}{
	Initial point pattern configuration. Incompatible with
	\code{n.start}.

	\code{x.start} may be a point pattern (an object
	of class \code{ppp}), or an object which can be coerced
	to this class by \code{\link{as.ppp}}, or a dataset containing
	vectors \code{x} and \code{y}.  

	If \code{x.start} is specified, then expansion of the
	simulation window (the argument \code{expand}
	of \code{\link{rmhcontrol}}) is not permitted.
      }
    }
    The parameters \code{n.start} and \code{x.start} are
    \emph{incompatible}.
}  
\seealso{
  \code{\link{rmh}},
  \code{\link{rmhcontrol}},
  \code{\link{rmhmodel}}
}
\examples{
   # 30 random points
   a <- rmhstart(n.start=30)

   # a particular point pattern
   data(cells)
   b <- rmhstart(x.start=cells)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

