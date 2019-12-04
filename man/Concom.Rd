\name{Concom}
\alias{Concom}
\title{The Connected Component Process Model}
\description{
  Creates an instance of the Connected Component point process model
  which can then be fitted to point pattern data.
}
\usage{
  Concom(r)
}
\arguments{
  \item{r}{Threshold distance}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the connected component process with disc radius \eqn{r}.
}
\details{
  This function defines the interpoint interaction structure of a point
  process called the connected component process.
  It can be used to fit this model to point pattern data.

  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the connected component interaction is
  yielded by the function \code{Concom()}. See the examples below.

  In \bold{standard form}, the connected component process
  (Baddeley and \Moller, 1989) with disc radius \eqn{r},
  intensity parameter \eqn{\kappa}{\kappa} and interaction parameter
  \eqn{\gamma}{\gamma} is a point process with probability density
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \kappa^{n(x)} \gamma^{-C(x)}
  }{
    f(x[1],\ldots,x[n]) =
    \alpha . \kappa^n(x) . \gamma^(-C(x))
  }
  for a point pattern \eqn{x}, where 
  \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, and \eqn{C(x)} is defined below.
  Here \eqn{\alpha}{\alpha} is a normalising constant.

  To define the term \code{C(x)}, suppose that we construct a planar
  graph by drawing an edge between 
  each pair of points \eqn{x_i,x_j}{x[i],x[j]} which are less than
  \eqn{r} units apart. Two points belong to the same connected component
  of this graph if they are joined by a path in the graph.
  Then \eqn{C(x)} is the number of connected components of the graph.

  The interaction parameter \eqn{\gamma}{\gamma} can be any positive number.
  If \eqn{\gamma = 1}{\gamma = 1} then the model reduces to a Poisson
  process with intensity \eqn{\kappa}{\kappa}.
  If \eqn{\gamma < 1}{\gamma < 1} then the process is regular,
  while if \eqn{\gamma > 1}{\gamma > 1} the process is clustered.
  Thus, a connected-component interaction process can be used to model either
  clustered or regular point patterns.
  
  In \pkg{spatstat}, the model is parametrised in a different form,
  which is easier to interpret.
  In \bold{canonical form}, the probability density is rewritten as
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \gamma^{-U(x)}
  }{
    f(x_1,\ldots,x_n) =
    \alpha . \beta^n(x) \gamma^(-U(x))
  }
  where \eqn{\beta}{\beta} is the new intensity parameter and
  \eqn{U(x) = C(x) - n(x)} is the interaction potential. 
  In this formulation, each isolated point of the pattern contributes a
  factor \eqn{\beta}{\beta} to the probability density (so the
  first order trend is \eqn{\beta}{\beta}). The quantity 
  \eqn{U(x)} is a true interaction potential, in the sense that
  \eqn{U(x) = 0} if the point pattern \eqn{x} does not contain any
  points that lie close together.

  When a new point \eqn{u} is added to an existing point pattern
  \eqn{x}, the rescaled potential \eqn{-U(x)} increases by
  zero or a positive integer.
  The increase is zero if \eqn{u} is not close to any point of \eqn{x}.
  The increase is a positive integer \eqn{k} if there are
  \eqn{k} different connected components of \eqn{x} that lie close to \eqn{u}.
  Addition of the point
  \eqn{u} contributes a factor \eqn{\beta \eta^\delta}{\beta * \eta^\delta}
  to the probability density, where \eqn{\delta}{\delta} is the
  increase in potential.

  If desired, the original parameter \eqn{\kappa}{\kappa} can be recovered from
  the canonical parameter by \eqn{\kappa = \beta\gamma}{\kappa = \beta * \gamma}.

  The \emph{nonstationary} connected component process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{\beta(x[i])}
  of location, rather than a constant beta. 
 
  Note the only argument of \code{Concom()} is the threshold distance \code{r}.
  When \code{r} is fixed, the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(\beta)}
  and \eqn{\log(\gamma)}{log(\gamma)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{Concom()}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\section{Edge correction}{
  The interaction distance of this process is infinite.
  There are no well-established procedures for edge correction
  for fitting such models, and accordingly the model-fitting function
  \code{\link{ppm}} will give an error message saying that the user must
  specify an edge correction. A reasonable solution is
  to use the border correction at the same distance \code{r}, as shown in the
  Examples.
}
\examples{
   # prints a sensible description of itself
   Concom(r=0.1)

   # Fit the stationary connected component process to redwood data
   ppm(redwood, ~1, Concom(r=0.07), rbord=0.07)

   # Fit the stationary connected component process to `cells' data
   ppm(cells, ~1, Concom(r=0.06), rbord=0.06)
   # eta=0 indicates hard core process.

   # Fit a nonstationary connected component model
   # with log-cubic polynomial trend
   \dontrun{
   ppm(swedishpines, ~polynom(x/10,y/10,3), Concom(r=7), rbord=7)
   }
}
\references{
  Baddeley, A.J. and \Moller, J. (1989)
  Nearest-neighbour Markov point processes and random sets.
  \emph{International Statistical Review} \bold{57}, 89--121.
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}
