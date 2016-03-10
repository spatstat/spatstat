\name{rmhmodel.list}
\alias{rmhmodel.list}
\title{Define Point Process Model for Metropolis-Hastings Simulation.}
\description{
  Given a list of parameters,
  builds a description of a point process model
  for use in simulating the model by the Metropolis-Hastings
  algorithm. 
}
\usage{
   \method{rmhmodel}{list}(model, ...)
}
\arguments{
  \item{model}{A list of parameters. See Details.}
  \item{\dots}{
    Optional list of additional named parameters.
  }
}
\value{
  An object of class \code{"rmhmodel"}, which is essentially
  a validated list of parameter values for the model.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the model chosen.
}
\details{
  The generic function \code{\link{rmhmodel}} takes a
  description of a point process model in some format, and
  converts it into an object of class \code{"rmhmodel"}
  so that simulations of the model can be generated using
  the Metropolis-Hastings algorithm \code{\link{rmh}}. 
  
  This function \code{rmhmodel.list} is the method
  for lists. The argument \code{model} should be a named list of parameters
  of the form
  
  \code{list(cif, par, w, trend, types)}

  where \code{cif} and \code{par} are required and the others are
  optional. For details about these components,
  see \code{\link{rmhmodel.default}}.

  The subsequent arguments \code{\dots} (if any) may also
  have these names, and they will take precedence over
  elements of the list \code{model}. 
}  
\references{
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

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\seealso{
  \code{\link{rmhmodel}},
  \code{\link{rmhmodel.default}},
  \code{\link{rmhmodel.ppm}},
  \code{\link{rmh}},
  \code{\link{rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{ppm}},
  \code{\link{Strauss}},
  \code{\link{Softcore}},
  \code{\link{StraussHard}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{DiggleGratton}},
  \code{\link{PairPiece}}
}
\examples{
   # Strauss process:
   mod01 <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   mod01 <- rmhmodel(mod01)

   # Strauss with hardcore:
   mod04 <- list(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=owin(c(0,10),c(0,5)))
   mod04 <- rmhmodel(mod04)

   # Soft core:
   w    <- square(10)
   mod07 <- list(cif="sftcr",
                     par=list(beta=0.8,sigma=0.1,kappa=0.5),
                     w=w)
   mod07 <- rmhmodel(mod07)
   
   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- list(cif="straussm",
                     par=list(beta=beta,gamma=gmma,radii=r),
                     w=square(250))
   mod08 <- rmhmodel(mod08)

   # specify types
   mod09 <- rmhmodel(list(cif="straussm",
                     par=list(beta=beta,gamma=gmma,radii=r),
                     w=square(250),
                     types=c("A", "B")))

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
   mod10 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=ri,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4))
   mod10 <- rmhmodel(mod10)

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
   r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
   h <- 20*(r-0.05)
   h[r<0.05] <- 0
   h[r>0.10] <- 1
   mod17 <- list(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))
   mod17 <- rmhmodel(mod17)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}

