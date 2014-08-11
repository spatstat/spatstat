\name{infline}
\alias{infline}
\alias{plot.infline}
\alias{print.infline}
\title{Infinite Straight Lines}
\description{
  Define the coordinates of one or more straight lines in the plane
}
\usage{
infline(a = NULL, b = NULL, h = NULL, v = NULL, p = NULL, theta = NULL)
\method{print}{infline}(x, \dots)
\method{plot}{infline}(x, \dots)
}
\arguments{
  \item{a,b}{Numeric vectors of equal length giving the
    intercepts \eqn{a} and slopes \eqn{b} of the lines.
    Incompatible with \code{h,v,p,theta}
  }
  \item{h}{Numeric vector giving the positions of horizontal lines when
    they cross the \eqn{y} axis.
    Incompatible with \code{a,b,v,p,theta}
  }
  \item{v}{Numeric vector giving the positions of vertical lines when
    they cross the \eqn{x} axis.
    Incompatible with \code{a,b,h,p,theta}
  }
  \item{p,theta}{Numeric vectors of equal length
    giving the polar coordinates of the line.
    Incompatible with \code{a,b,h,v}
  }
  \item{x}{An object of class \code{"infline"}}
  \item{\dots}{
    Extra arguments passed to \code{\link[base]{print}}
    for printing or \code{\link[graphics]{abline}} for plotting
  }
}
\details{
  The class \code{infline} is a convenient way to handle
  infinite straight lines in the plane.

  The position of a line can be specified in several ways:
  \itemize{
    \item
    its intercept \eqn{a} and slope \eqn{b} 
    in the equation \eqn{y = a + b x}{y = a + b * x}
    can be used unless the line is vertical.
    \item
    for vertical lines we can use the
    position \eqn{v} where the line crosses the \eqn{y} axis
    \item
    for horizontal lines we can use the
    position \eqn{h} where the line crosses the \eqn{x} axis
    \item
    the polar coordinates \eqn{p} and \eqn{\theta}{theta}
    can be used for any line. The line equation is
    \deqn{
      y \cos\theta + x \sin\theta = p
    }{
      y * cos(theta) + x * sin(theta) = p
    }
  }
  The command \code{infline} will accept line coordinates in any
  of these formats. The arguments \code{a,b,h,v} have the same interpretation
  as they do in the line-plotting function
  \code{\link[graphics]{abline}}.

  The command \code{infline} converts between different coordinate
  systems (e.g. from \code{a,b} to \code{p,theta}) and returns an
  object of class \code{"infline"}
  that contains a representation of the lines in
  each appropriate coordinate system. This object can be printed
  and plotted.
}
\value{
  The value of \code{infline} is an object of class \code{"infline"}
  which is basically a data frame with columns \code{a,b,h,v,p,theta}.
  Each row of the data frame represents one line. 
  Entries may be \code{NA} if a coordinate is not applicable to
  a particular line. 
}
\examples{
  infline(a=10:13,b=1)
  infline(p=1:3, theta=pi/4)
  plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="", asp=1)
  plot(infline(p=0.4, theta=seq(0,pi,length=20)))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
