\name{miplot}
\alias{miplot}
\title{Morisita Index Plot}
\description{
  Displays the Morisita Index Plot of a spatial point pattern.
}
\usage{
miplot(X, ...)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}) or something
    acceptable to \code{\link{as.ppp}}.
  }
  \item{\dots}{Optional arguments to control the appearance of the plot.}
}
\details{
  Morisita (1959) defined an index of spatial aggregation for a spatial
  point pattern based on quadrat counts. The spatial domain of the point
  pattern is first divided into \eqn{Q} subsets (quadrats) of equal size and
  shape. The numbers of points falling in each quadrat are counted.
  Then the Morisita Index is computed as
  \deqn{
    \mbox{MI} = Q \frac{\sum_{i=1}^Q n_i (n_i - 1)}{N(N-1)}
  }{
    MI = Q * sum(n[i] (n[i]-1))/(N(N-1))
  }
  where \eqn{n_i}{n[i]} is the number of points falling in the \eqn{i}-th
  quadrat, and \eqn{N} is the total number of points.
  If the pattern is completely random, \code{MI} should be approximately
  equal to 1. Values of \code{MI} greater than 1 suggest clustering.

  The \emph{Morisita Index plot} is a plot of the Morisita Index
  \code{MI} against the linear dimension of the quadrats. 
  The point pattern dataset is divided into \eqn{2 \times 2}{2 * 2}
  quadrats, then \eqn{3 \times 3}{3 * 3} quadrats, etc, and the
  Morisita Index is computed each time. This plot is an attempt to
  discern different scales of dependence in the point pattern data.
}
\value{
  None.
}
\references{
  M. Morisita (1959) Measuring of the dispersion of individuals and
  analysis of the distributional patterns.
  Memoir of the Faculty of Science, Kyushu University, Series E: Biology. 
  \bold{2}: 215--235. 
}
\seealso{
  \code{\link{quadratcount}}
}
\examples{
 data(longleaf)
 miplot(longleaf)
 opa <- par(mfrow=c(2,3))
 data(cells)
 data(japanesepines)
 data(redwood)
 plot(cells)
 plot(japanesepines)
 plot(redwood)
 miplot(cells)
 miplot(japanesepines)
 miplot(redwood)
 par(opa)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
