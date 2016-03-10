\name{Kmark}
\alias{Kmark}
\alias{markcorrint}
\title{Mark-Weighted K Function}
\description{
  Estimates the mark-weighted \eqn{K} function
  of a marked point pattern.
}
\usage{
  Kmark(X, f = NULL, r = NULL,
        correction = c("isotropic", "Ripley", "translate"), ...,
        f1 = NULL, normalise = TRUE, returnL = FALSE, fargs = NULL)

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
  The functions \code{Kmark} and \code{markcorrint} are identical.
  (Eventually \code{markcorrint} will be deprecated.)
  
  The \emph{mark-weighted \eqn{K} function} \eqn{K_f(r)}{K[f](r)}
  of a marked point process (Penttinen et al, 1992)
  is a generalisation of Ripley's \eqn{K} function, in which the contribution
  from each pair of points is weighted by a function of their marks.
  If the marks of the two points are \eqn{m_1, m_2}{m1, m2} then
  the weight is proportional to \eqn{f(m_1, m_2)}{f(m1, m2)} where
  \eqn{f} is a specified \emph{test function}.

  The mark-weighted \eqn{K} function is defined so that
  \deqn{
    \lambda K_f(r) = \frac{C_f(r)}{E[ f(M_1, M_2) ]}
  }{
    lambda * K_f(r) = C[f](r)/E[f(M1, M2)]
  }
  where 
  \deqn{
    C_f(r) = 
    E \left[
    \sum_{x \in X}
    f(m(u), m(x))
    1{0 < ||u - x|| \le r}
    \;  \big| \;
    u \in X
    \right]
  }{
    C[f](r) = E[ sum[x in X] f(m(u), m(x)) 1(0 < d(u,x) <= r) | u in X]
  }
  for any spatial location \eqn{u} taken to be a typical point of
  the point process \eqn{X}. Here \eqn{||u-x||}{d(u,x)} is the
  euclidean distance between \eqn{u} and \eqn{x}, so that the sum
  is taken over all random points \eqn{x} that lie within a distance
  \eqn{r} of the point \eqn{u}. The function \eqn{C_f(r)}{C[f](r)} is
  the \emph{unnormalised} mark-weighted \eqn{K} function.
  To obtain \eqn{K_f(r)}{K[f](r)} we standardise \eqn{C_f(r)}{C[f](r)}
  by dividing by \eqn{E[f(M_1,M_2)]}{E[f(M1,M2)]}, the expected value of
  \eqn{f(M_1,M_2)}{f(M1,M2)} when \eqn{M_1}{M1} and \eqn{M_2}{M2} are
  independent random marks with the same distribution as the marks in
  the point process. 

  Under the hypothesis of random labelling, the
  mark-weighted \eqn{K} function 
  is equal to Ripley's \eqn{K} function,
  \eqn{K_f(r) = K(r)}{K[f](r) = K(r)}.

  The mark-weighted \eqn{K} function is sometimes called the 
  \emph{mark correlation integral} because it is related to the
  mark correlation function \eqn{k_f(r)}{k[f](r)}
  and the pair correlation function \eqn{g(r)} by
  \deqn{
    K_f(r) = 2 \pi \int_0^r s k_f(s) \, g(s) \, {\rm d}s
  }{
    K[f](r) = 2 * pi * integral[0,r] (s * k[f](s) * g(s) ) ds
  }
  See \code{\link{markcorr}} for a definition of the
  mark correlation function.

  Given a marked point pattern \code{X},
  this command computes edge-corrected estimates
  of the mark-weighted \eqn{K} function.
  If \code{returnL=FALSE} then the estimated
  function \eqn{K_f(r)}{K[f](r)} is returned;
  otherwise the function
  \deqn{
    L_f(r) = \sqrt{K_f(r)/\pi}
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
  estimates of the mark-weighted \eqn{K} function \eqn{K_f(r)}{K[f](r)}
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
    # mark correlation function
    ms <- Kmark(spruces)
    plot(ms)

    # (2) simulated data with independent marks
    X <- rpoispp(100)
    X <- X \%mark\% runif(npoints(X))
    Xc <- Kmark(X)
    plot(Xc)
    
    # MULTITYPE DATA:
    # Hughes' amacrine data
    # Cells marked as 'on'/'off'
    M <- Kmark(amacrine, function(m1,m2) {m1==m2},
                         correction="translate")
    plot(M)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
