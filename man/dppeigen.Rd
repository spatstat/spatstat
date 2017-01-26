\name{dppeigen}
\alias{dppeigen}
\title{Internal function calculating eig and index}
\description{This function is mainly for internal package use and is usually
not called by the user.}
\usage{dppeigen(model, trunc, Wscale, stationary = FALSE)}
\arguments{
  \item{model}{object of class \code{"detpointprocfamily"} }
  \item{trunc}{numeric giving the truncation}
  \item{Wscale}{numeric giving the scale of the window relative to a unit box}
  \item{stationary}{logical indicating whether the stationarity of the model should be used (only works in dimension 2).}
}
\value{A list}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{models}





