\name{Gfox}
\alias{Gfox}
\alias{Jfox}
\title{
  Foxall's Distance Functions
}
\description{
  Given a point pattern \code{X} and a spatial object \code{Y},
  compute estimates of Foxall's  \eqn{G} and \eqn{J} functions.
}
\usage{
Gfox(X, Y, r=NULL, breaks=NULL, correction=c("km", "rs", "han"), W, \dots)
Jfox(X, Y, r=NULL, breaks=NULL, correction=c("km", "rs", "han"), W, \dots,
     warn.trim=TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    from which distances will be measured.
  }
  \item{Y}{
    An object of class \code{"ppp"}, \code{"psp"} or \code{"owin"}
    to which distances will be measured. Alternatively a pixel image
    (class \code{"im"}) with logical values.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which \eqn{Gfox(r)} or \eqn{Jfox(r)}
    should be evaluated. There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    Optional.
    The edge correction(s) to be used to estimate
    \eqn{Gfox(r)} or \eqn{Jfox(r)}.
    A vector of character strings selected from
    \code{"none"}, \code{"rs"}, \code{"km"}, \code{"cs"}
    and \code{"best"}.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{W}{
    Optional. A window (object of class \code{"owin"})
    to be taken as the window of observation.
    The distribution function will be estimated from data inside \code{W}.
    The default is \code{W=Frame(Y)} when \code{Y} is a window,
    and \code{W=Window(Y)} otherwise.
  }
  \item{\dots}{
    Extra arguments affecting the discretisation of distances.
    These arguments are ignored by \code{Gfox}, but
    \code{Jfox} passes them to \code{\link{Hest}} to determine
    the discretisation of the spatial domain.
  }
  \item{warn.trim}{
    Logical value indicating whether a warning should be issued
    by \code{Jfox} when the window of \code{X} had to be trimmed
    in order to be a subset of the frame of \code{Y}.
  }
}
\details{
  Given a point pattern \code{X} and another spatial object \code{Y},
  these functions compute two nonparametric measures of association
  between \code{X} and \code{Y}, introduced by Foxall
  (Foxall and Baddeley, 2002).
  
  Let the random variable \eqn{R} be the distance from a typical point
  of \code{X} to the object \code{Y}.
  Foxall's \eqn{G}-function is the cumulative distribution function
  of \eqn{R}:
  \deqn{G(r) = P(R \le r)}{P(R <= r)}
  
  Let the random variable \eqn{S} be the distance from a \emph{fixed} point
  in space to the object \code{Y}. The cumulative distribution function
  of \eqn{S} is the (unconditional) spherical contact distribution
  function
  \deqn{H(r) = P(S \le r)}{H(r) = P(S <= r)}
  which is computed by \code{\link{Hest}}.

  Foxall's \eqn{J}-function is the ratio
  \deqn{
    J(r) = \frac{1-G(r)}{1-H(r)}
  }{
    J(r) = (1-G(r))/(1-H(r))
  }
  For further interpretation, see Foxall and Baddeley (2002).
  
  Accuracy of \code{Jfox} depends on the pixel resolution,
  which is controlled by the
  arguments \code{eps}, \code{dimyx} and \code{xy} passed to
  \code{\link{as.mask}}. For example, use \code{eps=0.1} to specify
  square pixels of side 0.1 units, and \code{dimyx=256} to specify a
  256 by 256 grid of pixels.
}
\value{
  A function value table (object of class \code{"fv"})
  which can be printed, plotted, or converted to a data frame of values.
}
\references{
  Foxall, R. and Baddeley, A. (2002)
  Nonparametric measures of association between a
  spatial point process and a random set, with
  geological applications. \emph{Applied Statistics} \bold{51}, 165--182.
}
\seealso{
  \code{\link{Gest}}, 
  \code{\link{Hest}}, 
  \code{\link{Jest}}, 
  \code{\link{Fest}}
}
\examples{
  data(copper)
  X <- copper$SouthPoints
  Y <- copper$SouthLines
  G <- Gfox(X,Y)
  J <- Jfox(X,Y, correction="km")
  \testonly{
  J <- Jfox(X,Y, correction="km", eps=1)
  }
  \dontrun{
  J <- Jfox(X,Y, correction="km", eps=0.25)
  }
}
\author{Rob Foxall and
  \adrian
  
  
}
\keyword{spatial}
\keyword{nonparametric}
