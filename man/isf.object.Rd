\name{isf.object}
\alias{isf.object} %DoNotExport
\title{Interaction Structure Family Objects}
\description{
  Objects of class \code{"isf"} are used internally by the
  \pkg{spatstat} package to represent the structure
  of the interpoint interactions in a family of point process models.
}
\details{
  \emph{Advanced Use Only!}
  
  An object of class \code{"isf"} (Interaction Structure Family)
  is used internally by the \pkg{spatstat} package to
  represent the common mathematical and algorithmic structure of
  the interpoint interactions in a family of point process models.

  The existing objects of class \code{"isf"} are:
  \tabular{ll}{
    \code{\link{pairwise.family}} \tab
    pairwise interaction \cr
    \code{\link{triplet.family}} \tab
    triplet interaction \cr
    \code{\link{pairsat.family}} \tab
    saturated pairwise interaction \cr
    \code{\link{hierpair.family}} \tab
    hierarchical pairwise interaction \cr
    \code{\link{inforder.family}} \tab
    infinite order interaction \cr
    \code{\link{hybrid.family}} \tab
    hybrids of several interactions \cr
    \code{\link{ord.family}} \tab Ord interactions \cr
  }

  The information contained in these objects enables the
  \pkg{spatstat} package to select the appropriate algorithm
  for fitting, predicting and simulating each point process model.

  For example, in order to fit a model that involves pairwise
  interactions, the model-fitting function \code{\link{ppm}} would use
  information contained in \code{\link{pairwise.family}} to select the
  appropriate algorithms.
  
  An object of class \code{"isf"} is essentially a list of functions
  for various tasks.
  The internal format is undocumented and may be changed
  without notice.
}
\value{
  An object of class \code{"isf"}, essentially a list of functions
  for various tasks.

  The internal format is undocumented and may be changed
  without notice.
}
\author{\adrian.}
\keyword{spatial}

 
 
