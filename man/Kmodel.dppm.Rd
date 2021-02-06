\name{Kmodel.dppm}
\alias{Kmodel.detpointprocfamily}
\alias{pcfmodel.detpointprocfamily}
\alias{Kmodel.dppm}
\alias{pcfmodel.dppm}
\title{
  K-function or Pair Correlation Function of a
  Determinantal Point Process Model
}
\description{Returns the theoretical \eqn{K}-function
  or theoretical pair correlation function of a determinantal point
  process model as a function of one argument \eqn{r}.
}
\usage{
   \method{Kmodel}{dppm}(model, \dots)

   \method{pcfmodel}{dppm}(model, \dots)

   \method{Kmodel}{detpointprocfamily}(model, \dots)

   \method{pcfmodel}{detpointprocfamily}(model, \dots)
}
\arguments{
  \item{model}{Model of class \code{"detpointprocfamily"} or \code{"dppm"}.}
  \item{\dots}{Ignored (not quite true -- there is some undocumented internal use)}
}
\value{
  A function in the \R language, with one numeric argument \code{r},
  that can be used to evaluate the theoretical \eqn{K}-function or
  pair correlation function of the model at distances \code{r}.
}
\author{
  \spatstatAuthors.
}
\examples{
  model <- dppMatern(lambda=100, alpha=.01, nu=1, d=2)
  KMatern <- Kmodel(model)
  pcfMatern <- pcfmodel(model)
  plot(KMatern, xlim = c(0,0.05))
  plot(pcfMatern, xlim = c(0,0.05))
}
