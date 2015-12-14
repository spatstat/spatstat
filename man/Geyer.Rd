\name{Geyer}
\alias{Geyer}
\title{Geyer's Saturation Point Process Model}
\description{
Creates an instance of Geyer's saturation point process model
which can then be fitted to point pattern data.
}
\usage{
  Geyer(r,sat)
}
\arguments{
  \item{r}{Interaction radius. A positive real number.}
  \item{sat}{Saturation threshold. A positive real number.}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of Geyer's saturation point process
  with interaction radius \eqn{r} and saturation threshold \code{sat}.
}
\details{
  Geyer (1999) introduced the \dQuote{saturation process},
  a modification of the Strauss process (see \code{\link{Strauss}})
  in which the total contribution
  to the potential  from each point (from its pairwise interaction with all 
  other points) is trimmed to a maximum value \eqn{s}. 
  This model is implemented in the function \code{\link{Geyer}()}.

  The saturation point process with interaction radius \eqn{r},
  saturation threshold \eqn{s}, and 
  parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
  is the point process
  in which each point
  \eqn{x_i}{x[i]} in the pattern \eqn{X}
  contributes a factor
  \deqn{
    \beta \gamma^{\min(s, t(x_i, X))}
  }{
    beta gamma^min(s, t(x[i],X))
  }
  to the probability density of the point pattern,
  where \eqn{t(x_i, X)}{t(x[i],X)} denotes the
  number of \sQuote{close neighbours} of \eqn{x_i}{x[i]} in the pattern
  \eqn{X}. A close neighbour of \eqn{x_i}{x[i]} is a point
  \eqn{x_j}{x[j]} with \eqn{j \neq i}{j != i}
  such that the distance between 
  \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]} is less than or equal to \eqn{r}.

  If the saturation threshold \eqn{s} is set to infinity,
  this model reduces to the Strauss process (see \code{\link{Strauss}})
  with interaction parameter \eqn{\gamma^2}{gamma^2}.
  If \eqn{s = 0}, the model reduces to the Poisson point process.
  If \eqn{s} is a finite positive number, then the interaction parameter
  \eqn{\gamma}{gamma} may take any positive value (unlike the case
  of the Strauss process), with
  values \eqn{\gamma < 1}{gamma < 1}
  describing an \sQuote{ordered} or \sQuote{inhibitive} pattern,
  and 
  values \eqn{\gamma > 1}{gamma > 1}
  describing a \sQuote{clustered} or \sQuote{attractive} pattern.
 
  The nonstationary saturation process is similar except that 
  the value \eqn{\beta}{beta} 
  is replaced by a function \eqn{\beta(x_i)}{beta(x[i])}
  of location.
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the saturation process interaction is
  yielded by \code{Geyer(r, sat)} where the
  arguments \code{r} and \code{sat} specify
  the Strauss interaction radius \eqn{r} and the saturation threshold
  \eqn{s}, respectively. See the examples below.
 
  Note the only arguments are the interaction radius \code{r}
  and the saturation threshold \code{sat}.
  When \code{r} and \code{sat} are fixed,
  the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(beta)}
  and \eqn{\log(\gamma)}{log(gamma)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{Geyer()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{Strauss}},
  \code{\link{SatPiece}}
}
\references{
Geyer, C.J. (1999)
Likelihood Inference for Spatial Point Processes.
Chapter 3 in 
O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. Van Lieshout (eds)
\emph{Stochastic Geometry: Likelihood and Computation},
Chapman and Hall / CRC, 
Monographs on Statistics and Applied Probability, number 80.
Pages 79--140.
}

\examples{
   data(cells) 
   ppm(cells, ~1, Geyer(r=0.07, sat=2))
   # fit the stationary saturation process to `cells'
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
