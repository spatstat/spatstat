\name{copper}
\alias{copper}
\docType{data}
\title{
  Berman-Huntington points and lines data
}
\description{
  These data come from an intensive geological survey of
  a 70 x 158 km region in central Queensland, Australia.
  They consist of 67 points representing copper ore deposits,
  and 146 line segments representing geological `lineaments'.
  Lineaments are linear features, visible on a satellite image,
  that are believed to consist largely of geological faults (Berman, 1986,
  p. 55).
  It would be of great interest to predict the occurrence of copper deposits
  from the lineament pattern, since the latter can easily be observed on
  satellite images. 

  These data were introduced and analysed by Berman (1986).
  They have also been studied by Berman and Diggle (1989),
  Berman and Turner (1992),
  Baddeley and Turner (2000, 2005), Foxall and Baddeley (2002)
  and Baddeley et al (2005).
  
  Many analyses have been performed on the southern half of the data only.
  This subset is also provided.
} 
\format{
  \code{copper} is a list with the following entries:
  \describe{
    \item{Points}{a point pattern (object of class \code{"ppp"})
      representing the full point pattern of copper deposits.
      See \code{\link{ppp.object}} for details of the format.
    }
    \item{Lines}{a line segment pattern (object of class \code{"psp"})
      representing the lineaments in the full dataset.
      See \code{\link{psp.object}} for details of the format.
    }
    \item{SouthWindow}{the window delineating the southern half of
      the study region. An object of class \code{"owin"}.
    }
    \item{SouthPoints}{the point pattern of copper deposits in the
      southern half of the study region. An object of class
      \code{"ppp"}.
    }
    \item{SouthLines}{the line segment pattern of the lineaments in the
      southern half of the study region. An object of class \code{"psp"}.
    }
  }
}
\usage{data(copper)}
\examples{

  data(copper)

  # Plot full dataset

  plot(copper$Points)
  plot(copper$Lines, add=TRUE)

  # Plot southern half of data
  plot(copper$SouthPoints)
  plot(copper$SouthLines, add=TRUE)

  \dontrun{
    Z <- distmap(copper$SouthLines)
    plot(Z)
    X <- copper$SouthPoints
    ppm(X, ~D, covariates=list(D=Z))
  }
}
\source{
  Dr Jonathan Huntington, CSIRO Earth Science and Resource Engineering,
  Sydney, Australia.
  Coordinates kindly provided by Dr. Mark Berman
  and Dr. Andy Green, CSIRO, Sydney, Australia.
}
\references{
Baddeley, A. and Turner, R. (2000)
Practical maximum pseudolikelihood for spatial point patterns.
\emph{Australian and New Zealand Journal of Statistics}
\bold{42}, 283--322.
 
Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society, Series B}
\bold{67}, 617--666.

Baddeley, A. and Turner, R. (2005)
Modelling spatial point patterns in R.
In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
editors, \emph{Case Studies in Spatial Point Pattern Modelling},
Lecture Notes in Statistics number 185. Pages 23--74.
Springer-Verlag, New York, 2006. 
ISBN: 0-387-28311-0.  

Berman, M. (1986).
Testing for spatial association between a point process and another
  stochastic process.
\emph{Applied Statistics} \bold{35}, 54--62.

Berman, M. and Diggle, P.J. (1989)
Estimating Weighted Integrals of the Second-order Intensity of a
Spatial Point Process.
\emph{Journal of the Royal Statistical Society, series B}
\bold{51}, 81--92.

Berman, M. and Turner, T.R. (1992)
Approximating point process likelihoods with GLIM.
\emph{Applied Statistics} \bold{41}, 31--38.
 
Foxall, R. and Baddeley, A. (2002)
Nonparametric measures of association between a
spatial point process and a random set, with
geological applications. \emph{Applied Statistics} \bold{51}, 165--182.

}
\keyword{datasets}
\keyword{spatial}
