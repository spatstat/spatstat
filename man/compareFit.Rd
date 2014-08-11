\name{compareFit}
\alias{compareFit}
\title{
  Residual Diagnostics for Multiple Fitted Models
}
\description{
  Compares several fitted point process models using the
  same residual diagnostic. 
}
\usage{
compareFit(object, Fun, r = NULL, breaks = NULL, ...,
         trend = ~1, interaction = Poisson(), rbord = NULL,
         modelnames = NULL, same = NULL, different = NULL)
}
\arguments{
  \item{object}{
    Object or objects to be analysed.
    Either a fitted point process model (object of class \code{"ppm"}),
    a point pattern (object of class \code{"ppp"}),
    or a list of these objects.
  }
  \item{Fun}{
    Diagnostic function to be computed for each model.
    One of the functions \code{Kcom}, \code{Kres}, \code{Gcom},
    \code{Gres}, \code{psst}, \code{psstA} or \code{psstG}
    or a string containing one of these names.
}
  \item{r}{
    Optional. 
    Vector of values of the argument \eqn{r} at which the
    diagnostic should be computed.
    This argument is usually not specified. There is a sensible default.
  }
  \item{breaks}{
    Optional alternative to \code{r} for advanced use. 
  }
  \item{\dots}{
    Extra arguments passed to \code{Fun}.
  }
  \item{trend,interaction,rbord}{
    Optional. Arguments passed to \code{\link{ppm}}
    to fit a point process model to the data,
    if \code{object} is a point pattern or list of point patterns.
    See \code{\link{ppm}} for details.
    Each of these arguments can be a list, specifying different
    \code{trend}, \code{interaction} and/or \code{rbord}
    values to be used to generate different fitted models.
  }
  \item{modelnames}{
    Character vector. Short descriptive names for the different models.
}
  \item{same,different}{
    Character strings or character vectors passed to
    \code{\link{collapse.fv}} to
    determine the format of the output.
  }
}
\details{
  This is a convenient way to collect diagnostic information
  for several different point process models fitted to the same
  point pattern dataset, or for point process models of the same form fitted to
  several different datasets, etc.

  The first argument, \code{object}, is usually a list of
  fitted point process models
  (objects of class \code{"ppm"}), obtained from the
  model-fitting function \code{\link{ppm}}.

  For convenience, \code{object} can also be a list of point patterns
  (objects of class \code{"ppp"}).
  In that case, point process models will be fitted to
  each of the point pattern datasets,
  by calling \code{\link{ppm}} using the arguments
  \code{trend} (for the first order trend),
  \code{interaction} (for the interpoint interaction)
  and \code{rbord} (for the erosion distance in the border correction
  for the pseudolikelihood). See \code{\link{ppm}} for details
  of these arguments.

  Alternatively \code{object} can be a single point pattern
  (object of class \code{"ppp"}) and one or more of the arguments
  \code{trend}, \code{interaction} or \code{rbord}
  can be a list. In this case, point process models will be fitted to
  the same point pattern dataset, using each of the model specifications
  listed.

  The diagnostic function \code{Fun} will be applied to each of the
  point process models. The results will be collected into a single
  function value table. The \code{modelnames} are used to label the
  results from each fitted model.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Ege Rubak, 
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Jesper Moller.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{Kcom}},
  \code{\link{Kres}},
  \code{\link{Gcom}},
  \code{\link{Gres}},
  \code{\link{psst}},
  \code{\link{psstA}},
  \code{\link{psstG}},
  \code{\link{collapse.fv}}
}
\examples{
   nd <- 40
   \testonly{
        nd <- 10
   }
   ilist <- list(Poisson(), Geyer(7, 2), Strauss(7))
   iname <- c("Poisson", "Geyer", "Strauss")
   \testonly{
      ilist <- ilist[c(1,3)]
      iname <- iname[c(1,3)]
   }
   K <- compareFit(swedishpines, Kcom, interaction=ilist, rbord=9, 
            same="iso", different="icom", modelnames=iname, nd=nd)
   K
}
\keyword{spatial}
\keyword{models}
