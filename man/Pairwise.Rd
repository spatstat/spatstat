\name{Pairwise}
\alias{Pairwise}
\title{Generic Pairwise Interaction model}
\description{
Creates an instance of a pairwise interaction point process model
which can then be fitted to point pattern data.
}
\usage{
  Pairwise(pot, name, par, parnames, printfun)
}
\arguments{
  \item{pot}{An R language function giving the user-supplied
    pairwise interaction potential.}
  \item{name}{Character string.}
  \item{par}{List of numerical values for irregular parameters}
  \item{parnames}{Vector of names of irregular parameters}
  \item{printfun}{Do not specify this argument: for internal use only.}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  This code constructs a member of the
  pairwise interaction family \code{\link{pairwise.family}}
  with arbitrary pairwise interaction potential given by
  the user.

  Each pair of points in the point pattern contributes a factor
  \eqn{h(d)} to the probability density, where \eqn{d} is the distance
  between the two points. The factor term \eqn{h(d)} is
  \deqn{h(d) = \exp(-\theta \mbox{pot}(d))}{h(d) = exp(-theta * pot(d))}
  provided \eqn{\mbox{pot}(d)}{pot(d)} is finite, 
  where \eqn{\theta}{theta} is the coefficient vector in the model.  

  The function \code{pot} must take as its first argument
  a matrix of interpoint distances, and evaluate the
  potential for each of these distances. The result must be
  either a matrix with the same dimensions as its input,
  or an array with its first two dimensions the same as its input
  (the latter case corresponds to a vector-valued potential).

  If irregular parameters are present, then the second argument
  to \code{pot} should be a vector of the same type as \code{par}
  giving those parameter values.

  The values returned by \code{pot} may be finite numeric values,
  or \code{-Inf} indicating a hard core (that is, the corresponding
  interpoint distance is forbidden). We define
  \eqn{h(d) = 0} if \eqn{\mbox{pot}(d) = -\infty}{pot(d) = -Inf}.
  Thus, a potential value of minus infinity is \emph{always} interpreted
  as corresponding to \eqn{h(d) = 0}, regardless of the sign
  and magnitude of \eqn{\theta}{theta}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\examples{
   #This is the same as StraussHard(r=0.7,h=0.05)
   strpot <- function(d,par) {
         r <- par$r
         h <- par$h
         value <- (d <= r)
         value[d < h] <- -Inf
         value
   }
   mySH <- Pairwise(strpot, "StraussHard process", list(r=0.7,h=0.05),
           c("interaction distance r", "hard core distance h"))
   data(cells)
   ppm(cells, ~ 1, mySH, correction="isotropic")

   # Fiksel (1984) double exponential interaction
   # see Stoyan, Kendall, Mecke 1987 p 161

   fikspot <- function(d, par) {
      r <- par$r
      h <- par$h
      zeta <- par$zeta
      value <- exp(-zeta * d)
      value[d < h] <- -Inf
      value[d > r] <- 0
      value
   }
   Fiksel <- Pairwise(fikspot, "Fiksel double exponential process",
                      list(r=3.5, h=1, zeta=1),
                      c("interaction distance r",
                        "hard core distance h",
                        "exponential coefficient zeta"))
   data(spruces)
   fit <- ppm(unmark(spruces), ~1, Fiksel, rbord=3.5)
   fit
   plot(fitin(fit), xlim=c(0,4))
   coef(fit)
   # corresponding values obtained by Fiksel (1984) were -1.9 and -6.0
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{models}

