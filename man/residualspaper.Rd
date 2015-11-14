\name{residualspaper}
\alias{residualspaper}
\docType{data}
\title{
  Data and Code From JRSS Discussion Paper on Residuals
}
\description{
  This dataset contains the point patterns
  used as examples in the paper of Baddeley et al (2005).
  [Figure 2 is already available in \pkg{spatstat}
  as the \code{\link{copper}} dataset.] 

  R code is also provided to reproduce all
  the Figures displayed in Baddeley et al (2005).
  The component \code{plotfig}
  is a function, which can be called
  with a numeric or character argument specifying the Figure or Figures
  that should be plotted. See the Examples.
} 
\format{
  \code{residualspaper} is a list with the following components:
  \describe{
    \item{Fig1}{
      The locations of Japanese pine seedlings and saplings
      from Figure 1 of the paper.
      A point pattern (object of class \code{"ppp"}).
    }
    \item{Fig3}{
      The Chorley-Ribble data from Figure 3 of the paper.
      A list with three components, \code{lung}, \code{larynx}
      and \code{incin}. Each is a matrix with 2 columns
      giving the coordinates of the lung cancer cases,
      larynx cancer cases, and the incinerator, respectively.
      Coordinates are Eastings and Northings in km.
    }
    \item{Fig4a}{
      The synthetic dataset in Figure 4 (a) of the paper.
    }
    \item{Fig4b}{
      The synthetic dataset in Figure 4 (b) of the paper.
    }
    \item{Fig4c}{
      The synthetic dataset in Figure 4 (c) of the paper.
    }
    \item{Fig11}{
      The covariate displayed in Figure 11. A pixel image (object of
      class \code{"im"}) whose pixel values are distances to the
      nearest line segment in the \code{copper} data.
    }
    \item{plotfig}{A function which will compute and plot
      any of the Figures from the paper. The argument of
      \code{plotfig} is either a numeric vector or a character vector,
      specifying the Figure or Figures to be plotted. See the Examples.
    }
  }
}
\usage{data(residualspaper)}
\examples{

\dontrun{
  data(residualspaper)
  
  X <- residualspaper$Fig4a
  summary(X)
  plot(X)

  # reproduce all Figures
  residualspaper$plotfig()

  # reproduce Figures 1 to 10
  residualspaper$plotfig(1:10)

  # reproduce Figure 7 (a)
  residualspaper$plotfig("7a")
}
}
\source{
  Figure 1: Prof M. Numata. Data kindly supplied by Professor Y. Ogata
  with kind permission of Prof M. Tanemura.
  
  Figure 3:  Professor P.J. Diggle (rescaled by \adrian)

  Figure 4 (a,b,c): \adrian
}
\references{
  Baddeley, A., Turner, R., \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\keyword{datasets}
\keyword{spatial}
\keyword{models}
