\name{dendrite}
\alias{dendrite}
\docType{data}
\title{
  Dendritic Spines Data
}
\description{
  Dendrites are branching filaments which extend from the
  main body of a neuron (nerve cell) to propagate electrochemical
  signals. Spines are small protrusions on the dendrites.

  This dataset gives the locations of 566 spines
  observed on one branch of the dendritic tree of a rat neuron.
  The spines are classified according to their shape into three types:
  mushroom, stubby or thin.

  The data have been analysed in Jammalamadaka et al (2013) and
  Baddeley et al (2014). Please cite these papers and
  acknowledge the Kosik Lab, UC Santa Barbara, in any use of the data.
}
\usage{data("dendrite")}
\format{
  Object of class \code{"lpp"}. 
  See \code{\link{lpp}}.
}
\source{
  Kosik Lab, UC Santa Barbara (Dr Kenneth Kosik, Dr Sourav Banerjee).
  Formatted for \code{spatstat} by Dr Aruna Jammalamadaka.
}
\references{
Baddeley, A., Jammalamadaka, A. and Nair, G. (2014)
Multitype point process analysis of spines on the
dendrite network of a neuron.
\emph{Applied Statistics (Journal of the Royal Statistical
 Society, Series C)}, In press.
\code{doi: 10.1111/rssc.12054}

Jammalamadaka, A., Banerjee, S., Manjunath, B.S. and Kosik, K. (2013)
Statistical Analysis of Dendritic Spine Distributions in
Rat Hippocampal Cultures.
\emph{BMC Bioinformatics} \bold{14}, 287.
}
\examples{
plot(dendrite,leg.side="bottom", main="", cex=0.75, cols=2:4)
}
\keyword{datasets}
