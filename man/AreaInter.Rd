\name{AreaInter}
\alias{AreaInter}
\title{The Area Interaction Point Process Model}
\description{
  Creates an instance of the Area Interaction point process model
  (Widom-Rowlinson penetrable spheres model) 
  which can then be fitted to point pattern data.
}
\usage{
  AreaInter(r)
}
\arguments{
  \item{r}{The radius of the discs in the area interaction process}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the area-interaction process with disc radius \eqn{r}.
}
\details{
  This function defines the interpoint interaction structure of a point
  process called the Widom-Rowlinson penetrable sphere model
  or area-interaction process. It can be used to fit this model to
  point pattern data.

  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the area interaction structure is
  yielded by the function \code{AreaInter()}. See the examples below.

  In \bold{standard form}, the area-interaction process
  (Widom and Rowlinson, 1970;
  Baddeley and Van Lieshout, 1995) with disc radius \eqn{r},
  intensity parameter \eqn{\kappa}{kappa} and interaction parameter
  \eqn{\gamma}{gamma} is a point process with probability density
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \kappa^{n(x)} \gamma^{-A(x)}
  }{
    f(x[1],\ldots,x[n]) =
    alpha . kappa^n(x) . gamma^(-A(x))
  }
  for a point pattern \eqn{x}, where 
  \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, and \eqn{A(x)} is the area of the region formed by
  the union of discs of radius \eqn{r} centred at the points
  \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]}.
  Here \eqn{\alpha}{alpha} is a normalising constant.

  The interaction parameter \eqn{\gamma}{gamma} can be any positive number.
  If \eqn{\gamma = 1}{gamma = 1} then the model reduces to a Poisson
  process with intensity \eqn{\kappa}{kappa}.
  If \eqn{\gamma < 1}{gamma < 1} then the process is regular,
  while if \eqn{\gamma > 1}{gamma > 1} the process is clustered.
  Thus, an area interaction process can be used to model either
  clustered or regular point patterns. Two points interact if the
  distance between them is less than \eqn{2r}{2 * r}.

  The standard form of the model, shown above, is a little
  complicated to interpret in practical applications.
  For example, each isolated point of the pattern \eqn{x} contributes a factor
  \eqn{\kappa \gamma^{-\pi r^2}}{kappa * gamma^(-pi * r^2)}
  to the probability density. 
  
  In \pkg{spatstat}, the model is parametrised in a different form,
  which is easier to interpret.
  In \bold{canonical scale-free form}, the probability density is rewritten as
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \eta^{-C(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) eta^(-C(x))
  }
  where \eqn{\beta}{beta} is the new intensity parameter,
  \eqn{\eta}{eta} is the new interaction parameter, and
  \eqn{C(x) = B(x) - n(x)} is the interaction potential. Here
  \deqn{
    B(x) = \frac{A(x)}{\pi r^2}
  }{
    B(x) = A(x)/(pi * r^2)
  }
  is the normalised area (so that the discs have unit area).
  In this formulation, each isolated point of the pattern contributes a
  factor \eqn{\beta}{beta} to the probability density (so the
  first order trend is \eqn{\beta}{beta}). The quantity 
  \eqn{C(x)} is a true interaction potential, in the sense that
  \eqn{C(x) = 0} if the point pattern \eqn{x} does not contain any
  points that lie close together (closer than \eqn{2r}{2*r} units
  apart).

  When a new point \eqn{u} is added to an existing point pattern
  \eqn{x}, the rescaled potential \eqn{-C(x)} increases by
  a value between 0 and 1. 
  The increase is zero if \eqn{u} is not close to any point of \eqn{x}.
  The increase is 1 if the disc of radius \eqn{r} centred at \eqn{u}
  is completely contained in the union of discs of radius \eqn{r}
  centred at the data points \eqn{x_i}{x[i]}. Thus, the increase in
  potential is a measure of how close the new point \eqn{u} is to the
  existing pattern \eqn{x}. Addition of the point
  \eqn{u} contributes a factor \eqn{\beta \eta^\delta}{beta * eta^delta}
  to the probability density, where \eqn{\delta}{delta} is the
  increase in potential. 

  The old parameters \eqn{\kappa,\gamma}{kappa,gamma} of the
  standard form are related to
  the new parameters \eqn{\beta,\eta}{beta,eta} of the canonical
  scale-free form, by
  \deqn{
    \beta = \kappa \gamma^{-\pi r^2} = \kappa /\eta
  }{
    beta = kappa * gamma^(-pi * r^2)= kappa / eta
  }
  and
  \deqn{
    \eta = \gamma^{\pi r^2}
  }{
    eta = gamma^(pi * r^2)
  }
  provided \eqn{\gamma}{gamma} and \eqn{\kappa}{kappa} are positive and finite.

  In the canonical scale-free form, the parameter \eqn{\eta}{eta}
  can take any nonnegative value. The value \eqn{\eta = 1}{eta = 1}
  again corresponds to a Poisson process, with intensity \eqn{\beta}{beta}.
  If \eqn{\eta < 1}{eta < 1} then the process is regular,
  while if \eqn{\eta > 1}{eta > 1} the process is clustered.
  The value \eqn{\eta = 0}{eta = 0} corresponds to a hard core process
  with hard core radius \eqn{r} (interaction distance \eqn{2r}).
  
  The \emph{nonstationary} area interaction process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location, rather than a constant beta. 
 
  Note the only argument of \code{AreaInter()} is the disc radius \code{r}.
  When \code{r} is fixed, the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(beta)}
  and \eqn{\log(\eta)}{log(eta)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{AreaInter()}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\section{Warnings}{
  The interaction distance of this process is equal to
  \code{2 * r}. Two discs of radius \code{r} overlap if their centres
  are closer than \code{2 * r} units apart.

  The estimate of the interaction parameter \eqn{\eta}{eta}
  is unreliable if the interaction radius \code{r} is too small
  or too large. In these situations the model is approximately Poisson
  so that \eqn{\eta}{eta} is unidentifiable.
  As a rule of thumb, one can inspect the empty space function
  of the data, computed by \code{\link{Fest}}. The value \eqn{F(r)}
  of the empty space function at the interaction radius \code{r} should
  be between 0.2 and 0.8. 
}
\examples{
   \testonly{op <- spatstat.options(ngrid.disc=8)}

   # prints a sensible description of itself
   AreaInter(r=0.1)

   # Note the reach is twice the radius
   reach(AreaInter(r=1))

   # Fit the stationary area interaction process to Swedish Pines data
   data(swedishpines)
   ppm(swedishpines, ~1, AreaInter(r=7))

   # Fit the stationary area interaction process to `cells'
   data(cells) 
   ppm(cells, ~1, AreaInter(r=0.06))
   # eta=0 indicates hard core process.

   # Fit a nonstationary area interaction with log-cubic polynomial trend
   \dontrun{
   ppm(swedishpines, ~polynom(x/10,y/10,3), AreaInter(r=7))
   }
   \testonly{spatstat.options(op)}
}
\references{
  Baddeley, A.J. and Van Lieshout, M.N.M. (1995).
  Area-interaction point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{47} (1995) 601--619.

  Widom, B. and Rowlinson, J.S. (1970).
  New model for the study of liquid-vapor phase transitions.
  \emph{The Journal of Chemical Physics}
  \bold{52} (1970) 1670--1684.
}
\author{
  \adrian
  and
  \rolf
}
\keyword{spatial}
\keyword{models}
