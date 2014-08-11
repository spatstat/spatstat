\name{markcorrint}
\alias{markcorrint}
\title{Mark Correlation Integral}
\description{
  Estimates the mark correlation integral of a marked point pattern.
}
\usage{
  markcorrint(X, f = NULL, r = NULL,
              correction = c("isotropic", "Ripley", "translate"), ...,
              f1 = NULL, normalise = TRUE, returnL = FALSE, fargs = NULL)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. 
  }
  \item{f}{Optional. Test function \eqn{f} used in the definition of the
    mark correlation function. An \R function with at least two
    arguments. There is a sensible default.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the mark correlation function 
    \eqn{k_f(r)}{k[f](r)} should be evaluated.
    There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{\dots}{
    Ignored.
  }
  \item{f1}{
    An alternative to \code{f}. If this argument is given,
    then \eqn{f} is assumed to take the form
    \eqn{f(u,v)=f_1(u)f_1(v)}{f(u,v)=f1(u) * f1(v)}.
  }
  \item{normalise}{
    If \code{normalise=FALSE},
    compute only the numerator of the expression for the
    mark correlation.
  }
  \item{returnL}{
    Compute the analogue of the K-function if \code{returnL=FALSE}
    or the analogue of the L-function if \code{returnL=TRUE}.
  }
  \item{fargs}{
    Optional. A list of extra arguments to be passed to the function
    \code{f} or \code{f1}.
  }
}
\details{
  Given a marked point pattern \code{X},
  this command estimates the weighted indefinite integral
  \deqn{
    K_f(r) = 2 \pi \int_0^r s k_f(s) ds
  }{
    K[f](r) = 2 * pi * integral[0,r] (s * k[f](s)) ds
  }
  of the mark correlation function \eqn{k_f(r)}{k[f](r)}.
  See \code{\link{markcorr}} for a definition of the
  mark correlation function.

  The use of the weighted indefinite integral
  was advocated by Penttinen et al (1992).
  The relationship between
  \eqn{K_f}{K[f]} and \eqn{k_f}{k[f]} is analogous to the
  relationship between the classical K-function \eqn{K(r)} and the
  pair correlation function \eqn{g(r)}.

  If \code{returnL=FALSE} then the function \eqn{K_f(r)}{K[f](r)} is returned;
  otherwise the function
  \deqn{
    L_f(r) = \sqrt{K_f(r)/pi}
  }{
    L[f](r) = sqrt(K[f](r)/pi)
  }
  is returned.
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the mark correlation integral \eqn{K_f(r)}{K[f](r)}
    has been  estimated
  }
  \item{theo}{the theoretical value of \eqn{K_f(r)}{K[f](r)}
    when the marks attached to different points are independent,
    namely \eqn{\pi r^2}{pi * r^2}
  }
  together with a column or columns named 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the mark correlation integral \eqn{K_f(r)}{K[f](r)}
  obtained by the edge corrections named (if \code{returnL=FALSE}).
}
\references{
  Penttinen, A., Stoyan, D. and Henttonen, H. M. (1992)
  Marked point processes in forest statistics.
  \emph{Forest Science} \bold{38} (1992) 806-824.

  Illian, J., Penttinen, A., Stoyan, H. and Stoyan, D. (2008)
  \emph{Statistical analysis and modelling of spatial point patterns}.
  Chichester: John Wiley.
}
\seealso{
  \code{\link{markcorr}} to estimate the mark correlation function.
}
\examples{
    # CONTINUOUS-VALUED MARKS:
    # (1) Spruces
    # marks represent tree diameter
    data(spruces)
    # mark correlation function
    ms <- markcorrint(spruces)
    plot(ms)

    # (2) simulated data with independent marks
    X <- rpoispp(100)
    X <- X \%mark\% runif(X$n)
    Xc <- markcorrint(X)
    plot(Xc)
    
    # MULTITYPE DATA:
    # Hughes' amacrine data
    # Cells marked as 'on'/'off'
    data(amacrine)
    M <- markcorrint(amacrine, function(m1,m2) {m1==m2},
                  correction="translate")
    plot(M)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
