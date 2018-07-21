\name{rmhmodel.default}
\alias{rmhmodel.default}
\title{Build Point Process Model for Metropolis-Hastings Simulation.}
\description{
  Builds a description of a point process model
  for use in simulating the model by the Metropolis-Hastings
  algorithm. 
}
\usage{
  \method{rmhmodel}{default}(..., 
         cif=NULL, par=NULL, w=NULL, trend=NULL, types=NULL)
}
\arguments{
  \item{\dots}{Ignored.}
  \item{cif}{Character string specifying the choice of model}
  \item{par}{Parameters of the model}
  \item{w}{Spatial window in which to simulate}
  \item{trend}{Specification of the trend in the model}
  \item{types}{A vector of factor levels defining the possible
    marks, for a multitype process.
  }
}
\value{
  An object of class \code{"rmhmodel"}, which is essentially
  a list of parameter values for the model.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the model chosen.
}
\details{
  The generic function \code{\link{rmhmodel}} takes a
  description of a point process model in some format, and
  converts it into an object of class \code{"rmhmodel"}
  so that simulations of the model can be generated using
  the Metropolis-Hastings algorithm \code{\link{rmh}}. 
  
  This function \code{rmhmodel.default} is the default method.
  It builds a description of the point process model
  from the simple arguments listed.

  The argument \code{cif} is a character string specifying the choice of
  interpoint interaction for the point process. The current options are
  \describe{
    \item{\code{'areaint'}}{Area-interaction process.}
    \item{\code{'badgey'}}{Baddeley-Geyer (hybrid Geyer) process.}
    \item{\code{'dgs'}}{Diggle, Gates and Stibbard (1987) process}
    \item{\code{'diggra'}}{Diggle and Gratton (1984) process}
    \item{\code{'fiksel'}}{Fiksel double exponential process (Fiksel, 1984).}
    \item{\code{'geyer'}}{Saturation process (Geyer, 1999).}
    \item{\code{'hardcore'}}{Hard core process}
    \item{\code{'lennard'}}{Lennard-Jones process}
    \item{\code{'lookup'}}{General isotropic pairwise interaction process,
      with the interaction function specified via a ``lookup table''.}
    \item{\code{'multihard'}}{Multitype hardcore process}
    \item{\code{'penttinen'}}{The Penttinen process}
    \item{\code{'strauss'}}{The Strauss process}
    \item{\code{'straush'}}{The Strauss process with hard core}
    \item{\code{'sftcr'}}{The Softcore process}
    \item{\code{'straussm'}}{ The multitype Strauss process}
    \item{\code{'straushm'}}{Multitype Strauss process with hard core}
    \item{\code{'triplets'}}{Triplets process (Geyer, 1999).}
  }
  It is also possible to specify a \emph{hybrid} of these interactions
  in the sense of Baddeley et al (2013).
  In this case, \code{cif} is a character vector containing names from
  the list above. For example, \code{cif=c('strauss', 'geyer')} would
  specify a hybrid of the Strauss and Geyer models.
  
  The argument \code{par} supplies parameter values appropriate to
  the conditional intensity function being invoked.
  For the interactions listed above, these parameters are:
  \describe{
    \item{areaint:}{
      (Area-interaction process.) A \bold{named} list with components
      \code{beta,eta,r} which are respectively the ``base''
      intensity, the scaled interaction parameter and the
      interaction radius.  
    }
    \item{badgey:}{
      (Baddeley-Geyer process.)
      A \bold{named} list with components
      \code{beta} (the ``base'' intensity), \code{gamma} (a vector
      of non-negative interaction parameters), \code{r} (a vector
      of interaction radii, of the same length as \code{gamma},
      in \emph{increasing} order), and \code{sat} (the saturation
      parameter(s); this may be a scalar, or a vector of the same
      length as \code{gamma} and \code{r}; all values should be at
      least 1).  Note that because of the presence of ``saturation''
      the \code{gamma} values are permitted to be larger than 1.
    }
    \item{dgs:}{
      (Diggle, Gates, and Stibbard process.
      See Diggle, Gates, and Stibbard (1987))
      A \bold{named} list with components
      \code{beta} and \code{rho}.  This process has pairwise interaction
      function equal to
      \deqn{
	e(t) = \sin^2\left(\frac{\pi t}{2\rho}\right)
      }{
	e(t) = sin^2((pi * t)/(2 * rho))
      }
      for \eqn{t < \rho}{t < rho}, and equal to 1
      for \eqn{t \ge \rho}{t >= rho}.
    }
    \item{diggra:}{
      (Diggle-Gratton process. See Diggle and Gratton (1984)
      and Diggle, Gates and Stibbard (1987).)
      A \bold{named} list with components \code{beta},
      \code{kappa}, \code{delta} and \code{rho}.  This process has
      pairwise interaction function \eqn{e(t)} equal to 0
      for \eqn{t < \delta}{t < delta}, equal to
      \deqn{
	\left(\frac{t-\delta}{\rho-\delta}\right)^\kappa
      }{
	((t-delta)/(rho-delta))^kappa
      }
      for \eqn{\delta \le t < \rho}{delta <= t < rho},
      and equal to 1 for \eqn{t \ge  \rho}{t >= rho}.
      Note that here we use the symbol
      \eqn{\kappa}{kappa} where Diggle, Gates, and Stibbard use
      \eqn{\beta}{beta} since we reserve the symbol \eqn{\beta}{beta}
      for an intensity parameter.
    }
    \item{fiksel:}{
      (Fiksel double exponential process, see Fiksel (1984))
      A \bold{named} list with components \code{beta},
      \code{r}, \code{hc}, \code{kappa} and \code{a}.  This process has
      pairwise interaction function \eqn{e(t)} equal to 0
      for \eqn{t < hc}, equal to
      \deqn{
	\exp(a \exp(- \kappa t))
      }{
	exp(a * exp( - kappa * t))
      }
      for \eqn{hc \le t < r}{hc <= t < r},
      and equal to 1 for \eqn{t \ge  r}{t >= r}.
    }
    \item{geyer:}{
      (Geyer's saturation process. See Geyer (1999).)
      A \bold{named} list
      with components \code{beta}, \code{gamma}, \code{r}, and \code{sat}.
      The components \code{beta}, \code{gamma}, \code{r} are as for
      the Strauss model, and \code{sat} is the ``saturation''
      parameter.  The model is Geyer's ``saturation'' point process
      model, a modification of the Strauss process in which
      we effectively impose an upper limit (\code{sat}) on the number of
      neighbours which will be counted as close to a given point.
      
      Explicitly, a saturation point process with interaction
      radius \eqn{r}, saturation threshold \eqn{s}, and
      parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
      is the point process in which each point \eqn{x_i}{x[i]}
      in the pattern \eqn{X} contributes a factor
      \deqn{\beta \gamma^{\min(s, t(x_i,X))}}{beta gamma^min(s,t(x[i],X))}
      to the probability density of the point pattern,
      where \eqn{t(x_i,X)}{t(x[i],X)} denotes the number of
      ``\eqn{r}-close neighbours'' of \eqn{x_i}{x[i]} in the
      pattern \eqn{X}.

      If the saturation threshold \eqn{s} is infinite,
      the Geyer process reduces to a Strauss process
      with interaction parameter \eqn{\gamma^2}{gamma^2}
      rather than \eqn{\gamma}{gamma}.
    }
    \item{hardcore:}{
      (Hard core process.) A \bold{named} list
      with components \code{beta} and \code{hc}
      where \code{beta} is the base intensity and \code{hc} is the
      hard core distance.
      This process has pairwise interaction function \eqn{e(t)}
      equal to 1 if \eqn{t > hc} and 0 if \eqn{t <= hc}.
    }
    \item{lennard:}{
      (Lennard-Jones process.) A \bold{named} list
      with components \code{sigma} and \code{epsilon},
      where \code{sigma} is the characteristic diameter
      and \code{epsilon} is the well depth.
      See \code{\link{LennardJones}} for explanation.
    }
    \item{multihard:}{
      (Multitype hard core process.) A \bold{named} list
      with components \code{beta} and \code{hradii},
      where \code{beta} is a vector of base intensities for each type
      of point, and \code{hradii} is a matrix of hard core radii
      between each pair of types. 
    }
    \item{penttinen:}{
      (Penttinen process.) A \bold{named} list with components
      \code{beta,gamma,r} which are respectively the ``base''
      intensity, the pairwise interaction parameter, and the disc radius.
      Note that \code{gamma} must be less than or equal to 1.
      See \code{\link{Penttinen}} for explanation.
      (Note that there is also an algorithm for perfect simulation
      of the Penttinen process, \code{\link{rPenttinen}})
    }
    \item{strauss:}{
      (Strauss process.) A \bold{named} list with components
      \code{beta,gamma,r} which are respectively the ``base''
      intensity, the pairwise interaction parameter and the
      interaction radius.  Note that \code{gamma} must be less than
      or equal to 1.
      (Note that there is also an algorithm for perfect simulation
      of the Strauss process, \code{\link{rStrauss}})
    }
    \item{straush:}{
      (Strauss process with hardcore.) A \bold{named} list with
      entries \code{beta,gamma,r,hc} where \code{beta}, \code{gamma},
      and \code{r} are as for the Strauss process, and \code{hc} is
      the hardcore radius.  Of course \code{hc} must be less than
      \code{r}.
    }
    \item{sftcr:}{
      (Softcore process.) A \bold{named} list with components
      \code{beta,sigma,kappa}.  Again \code{beta} is a ``base''
      intensity. The pairwise interaction between two points
      \eqn{u \neq v}{u != v} is
      \deqn{
	\exp \left \{
	- \left ( \frac{\sigma}{||u-v||} \right )^{2/\kappa}
	\right \}
      }{-(sigma/||u-v||)^(2/kappa)}
      Note that it is necessary that \eqn{0 < \kappa < 1}{0 < kappa <1}.
    }
    \item{straussm:}{
      (Multitype Strauss process.) A \bold{named} list with components
      \itemize{
	\item
	\code{beta}: 
	A vector of ``base'' intensities, one for each possible type.
	\item
	\code{gamma}:
	A \bold{symmetric} matrix of interaction parameters,
	with \eqn{\gamma_{ij}}{gamma_ij} pertaining to the interaction between
	type \eqn{i} and type \eqn{j}.
	\item
	\code{radii}:
	A \bold{symmetric} matrix of interaction radii, with
	entries \eqn{r_{ij}}{r_ij} pertaining to the interaction between type
	\eqn{i} and type \eqn{j}.
      }
    }
    \item{straushm:}{
      (Multitype Strauss process with hardcore.)
      A \bold{named} list with components \code{beta} and \code{gamma}
      as for \code{straussm} and
      \bold{two} ``radii'' components:
      \itemize{
        \item \code{iradii}: the interaction radii
        \item \code{hradii}: the hardcore radii
      }
      which are both symmetric matrices of nonnegative numbers.
      The entries of \code{hradii} must be less than the
      corresponding entries
      of \code{iradii}.
    }
    \item{triplets:}{
      (Triplets process.) A \bold{named} list with components
      \code{beta,gamma,r} which are respectively the ``base''
      intensity, the triplet interaction parameter and the
      interaction radius.  Note that \code{gamma} must be less than
      or equal to 1.
    }
    \item{lookup:}{
      (Arbitrary pairwise interaction process with isotropic interaction.)
      A \bold{named} list with components
      \code{beta}, \code{r}, and \code{h}, or just with components
      \code{beta} and \code{h}.

      This model is the pairwise interaction process
      with an isotropic interaction given by any chosen function \eqn{H}.
      Each pair of points \eqn{x_i, x_j}{x[i], x[j]} in the
      point pattern contributes
      a factor \eqn{H(d(x_i, x_j))}{H(d(x[i],x[j]))}
      to the probability density, where \eqn{d} denotes distance
      and \eqn{H} is the pair interaction function.

      The component \code{beta} is a
      (positive) scalar which determines the ``base'' intensity
      of the process.

      In this implementation, \eqn{H} must be a step function.
      It is specified by the user in one of two ways.
      \itemize{
	\item
	\bold{as a vector of values:}
	If \code{r} is present, then \code{r} is assumed to 
	give the locations of jumps in the function \eqn{H},
	while the vector \code{h} gives the corresponding
	values of the function.

	Specifically, the interaction function
	\eqn{H(t)} takes the value \code{h[1]}
	for distances \eqn{t} in the interval 
	\code{[0, r[1])}; takes the value \code{h[i]}
	for distances \eqn{t} in the interval 
	\code{[r[i-1], r[i])} where
	\eqn{i = 2,\ldots, n}{i = 2, ..., n};
	and takes the value 1 for \eqn{t \ge r[n]}{t >= r[n]}.
	Here \eqn{n} denotes the length of \code{r}.
	    
	The components \code{r} and \code{h}
	must be numeric vectors of equal length.
	The \code{r} values must be strictly positive, and 
	sorted in increasing order.

	The entries of \code{h} must be non-negative. 
	If any entry of \code{h} is greater than 1,
	then the entry \code{h[1]} must be 0 (otherwise the specified
	process is non-existent).

	Greatest efficiency is achieved if the values of
	\code{r} are equally spaced.
	    
	[\bold{Note:} The usage of \code{r} and \code{h}
	has \emph{changed} from the previous usage in \pkg{spatstat}
	versions 1.4-7 to 1.5-1, in which ascending order was not required,
	and in which the first entry of \code{r} had to be 0.]

	\item
	\bold{as a stepfun object:}
	If \code{r} is absent, then \code{h} must be
	an object of class \code{"stepfun"} specifying
	a step function. Such objects are created by
	\code{\link{stepfun}}. 

	The stepfun object \code{h} must be right-continuous
	(which is the default using \code{\link{stepfun}}.)

	The values of the step function must all be nonnegative.
	The values must all be less than 1
	unless the function is identically zero on some initial
	interval \eqn{[0,r)}. The rightmost value (the value of
	\code{h(t)} for large \code{t}) must be equal to 1.

	Greatest efficiency is achieved if the jumps (the
	``knots'' of the step function) are equally spaced.
      }
    }
  }
  For a hybrid model, the argument \code{par} should be a list,
  of the same length as \code{cif}, such that \code{par[[i]]}
  is a list of the parameters required for the interaction
  \code{cif[i]}. See the Examples.
  
  The optional argument \code{trend} determines the spatial trend in the model,
  if it has one. It should be a function or image
  (or a list of such, if the model is multitype)
  to provide the value of the trend at an arbitrary point.
  \describe{
    \item{trend given as a function:}{A trend
      function may be a function of any number of arguments,
      but the first two must be the \eqn{x,y} coordinates of
      a point.  Auxiliary arguments may be passed
      to the \code{trend} function at the time of simulation,
      via the \code{\dots} argument to \code{\link{rmh}}.
      
      The function \bold{must} be \bold{vectorized}.
      That is, it must be capable of accepting vector valued
      \code{x} and \code{y} arguments.  Put another way,
      it must be capable of calculating the trend value at a
      number of points, simultaneously, and should return the
      \bold{vector} of corresponding trend values.
    }
    \item{trend given as an image:}{
      An image (see \code{\link{im.object}})
      provides the trend values at a grid of
      points in the observation window and determines the trend
      value at other points as the value at the nearest grid point.
    }
  }
  Note that the trend or trends must be \bold{non-negative};
  no checking is done for this.
  
  The optional argument \code{w} specifies the window
  in which the pattern is to be generated.  If specified, it must be in
  a form which can be coerced to an object of class \code{owin}
  by \code{\link{as.owin}}.

  The optional argument \code{types} specifies the possible
  types in a multitype point process. If the model being simulated
  is multitype, and \code{types} is not specified, then this vector
  defaults to \code{1:ntypes} where \code{ntypes} is the number of
  types.
}
\references{
  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  \emph{Journal of Statistical Software} \bold{55}:11, 1--43.
  \url{http://www.jstatsoft.org/v55/i11/}

   Diggle, P. J. (2003) \emph{Statistical Analysis of Spatial Point
   Patterns} (2nd ed.) Arnold, London.

   Diggle, P.J. and Gratton, R.J. (1984)
   Monte Carlo methods of inference for implicit statistical models.
   \emph{Journal of the Royal Statistical Society, series B}
   \bold{46}, 193 -- 212.

   Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.

   Fiksel, T. (1984)
   Estimation of parameterized pair potentials
   of marked and non-marked Gibbsian point processes.
   \emph{Electronische Informationsverabeitung und Kybernetika}
   \bold{20}, 270--278.

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\section{Warnings in Respect of ``lookup''}{

For the \code{lookup} cif, 
the entries of the \code{r} component of \code{par}
must be \emph{strictly positive} and sorted into ascending order.

Note that if you specify the \code{lookup} pairwise interaction
function via \code{\link{stepfun}()} the arguments \code{x}
and \code{y} which are passed to \code{stepfun()} are slightly
different from \code{r} and \code{h}:  \code{length(y)} is equal
to \code{1+length(x)}; the final entry of \code{y} must be equal
to 1 --- i.e. this value is explicitly supplied by the user rather
than getting tacked on internally.

The step function returned by \code{stepfun()} must be right
continuous (this is the default behaviour of \code{stepfun()})
otherwise an error is given.
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{ppm}},
  \rmhInteractionsList.
}
\examples{
   # Strauss process:
   mod01 <- rmhmodel(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   mod01
   # The above could also be simulated using 'rStrauss'

   # Strauss with hardcore:
   mod04 <- rmhmodel(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=owin(c(0,10),c(0,5)))

   # Hard core:
   mod05 <- rmhmodel(cif="hardcore",par=list(beta=2,hc=0.3),
              w=square(5))

   # Soft core:
   w    <- square(10)
   mod07 <- rmhmodel(cif="sftcr",
                     par=list(beta=0.8,sigma=0.1,kappa=0.5),
                     w=w)
   
   # Area-interaction process:
   mod42 <- rmhmodel(cif="areaint",par=list(beta=2,eta=1.6,r=0.7),
                 w=c(0,10,0,10))

   # Baddeley-Geyer process:
   mod99 <- rmhmodel(cif="badgey",par=list(beta=0.3,
                     gamma=c(0.2,1.8,2.4),r=c(0.035,0.07,0.14),sat=5),
                     w=unit.square())

   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- rmhmodel(cif="straussm",
                     par=list(beta=beta,gamma=gmma,radii=r),
                     w=square(250))
   # specify types
   mod09 <- rmhmodel(cif="straussm",
                     par=list(beta=beta,gamma=gmma,radii=r),
                     w=square(250),
                     types=c("A", "B"))

   # Multitype Hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod08hard <- rmhmodel(cif="multihard",
                     par=list(beta=beta,hradii=rhc),
                     w=square(250),
                     types=c("A", "B"))

   
   # Multitype Strauss hardcore with trends for each type:
   beta  <- c(0.27,0.08)
   ri    <- matrix(c(45,45,45,45),2,2)
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   tr4   <- function(x,y){x <- x/250; y <- y/250;
                         exp(-0.6*x+0.5*y)}
                        # log linear trend
   mod10 <- rmhmodel(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=ri,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4))

   # Triplets process:
   mod11 <- rmhmodel(cif="triplets",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
      r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
      h <- 20*(r-0.05)
      h[r<0.05] <- 0
      h[r>0.10] <- 1
      mod17 <- rmhmodel(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))

  # hybrid model
  modhy <- rmhmodel(cif=c('strauss', 'geyer'),
                    par=list(list(beta=100,gamma=0.5,r=0.05),
                             list(beta=1, gamma=0.7,r=0.1, sat=2)),
                    w=square(1))
  modhy
}
\author{
  \adrian
  and
  \rolf
}
\keyword{spatial}
\keyword{datagen}

