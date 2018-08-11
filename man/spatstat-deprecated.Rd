\name{spatstat-deprecated}
\alias{as.psp.owin}
\alias{which.max.im}
\alias{delaunay.distance}
\alias{delaunay.network}
\alias{dirichlet.edges}
\alias{dirichlet.network}
\alias{dirichlet.vertices}
\alias{dirichlet.weights}
\alias{circumradius}
\alias{circumradius.owin}
\alias{circumradius.ppp}
\alias{circumradius.linnet}
\title{Deprecated spatstat functions}
\description{
  Deprecated spatstat functions.
}
\usage{
\method{as.psp}{owin}(x, \dots, window=NULL, 
       check=spatstat.options("checksegments"), fatal=TRUE)
which.max.im(x)
delaunay.distance(\dots)
delaunay.network(\dots)
dirichlet.edges(\dots)
dirichlet.network(\dots)
dirichlet.vertices(\dots)
dirichlet.weights(\dots)
circumradius(x, \dots)
\method{circumradius}{owin}(x, \dots)
\method{circumradius}{ppp}(x, \dots)
\method{circumradius}{linnet}(x, \dots)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat} package.
  
  \code{as.psp.owin} has been replaced by \code{\link{edges}}.
  
  \code{which.max.im(x)} is replaced by
  \code{\link{im.apply}(x, which.max)}.

  \code{delaunay.distance}   is replaced by \code{delaunayDistance}.
  
  \code{delaunay.network} is replaced by \code{delaunayNetwork}.
  
  \code{dirichlet.edges} is replaced by \code{dirichletEdges}.
  
  \code{dirichlet.network} is replaced by \code{dirichletNetwork}.
  
  \code{dirichlet.vertices} is replaced by \code{dirichletVertices}.
  
  \code{dirichlet.weights} is replaced by \code{dirichletWeights}.

  \code{circumradius} is replaced by the more appropriately named
  \code{boundingradius}.
}
\keyword{internal}
