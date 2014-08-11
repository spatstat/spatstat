\name{spruces}
\alias{spruces}
\docType{data}
\title{
  Spruces Point Pattern
}
\description{
  The data give the locations of Norwegian spruce trees 
  in a natural forest stand in Saxonia, Germany.
  Each tree is marked with its diameter at breast height.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations
  in a 56 x 38 metre sampling region. Each tree is marked
  with its diameter at breast height. All values are given in metres.

  See \code{\link{ppp.object}} for details of the format of a
  point pattern object. The marks are numeric.

  These data have been analysed by Fiksel (1984, 1988),
  Stoyan et al (1987), Penttinen et al (1992) and
  Goulard et al (1996).
}
\usage{data(spruces)}
\source{Stoyan et al (1987). Original source unknown.}
\examples{
     data(spruces)
     plot(spruces)
     # To reproduce Goulard et al. Figure 3
     plot(spruces, maxsize=5*max(spruces$marks))
     plot(unmark(spruces), add=TRUE)
}
\references{
  Fiksel, T. (1984)
  Estimation of parameterized pair potentials of marked and
  nonmarked Gibbsian point processes.
  \emph{Elektron. Informationsverarb. u. Kybernet.}
  \bold{20}, 270--278.
  
  Fiksel, T. (1988)
  Estimation of interaction potentials of Gibbsian point processes.
  \emph{Statistics}
  \bold{19}, 77-86
  
  Goulard, M., S\"arkk\"a, A. and Grabarnik, P. (1996)
  Parameter estimation for marked Gibbs point
  processes through the maximum pseudolikelihood method.
  \emph{Scandinavian Journal of Statistics}
  \bold{23}, 365--379.

  Penttinen, A., Stoyan, D. and Henttonen, H. (1992)
  Marked point processes in forest statistics.
  \emph{Forest Science} \bold{38}, 806--824.

  Stoyan, D., Kendall, W.S. and Mecke, J. (1987)
  \emph{Stochastic Geometry and its Applications}.
  Wiley.
}
\keyword{datasets}
\keyword{spatial}
