\name{Kovesi}
\alias{Kovesi}
\docType{data}
\title{
  Colour Sequences with Uniform Perceptual Contrast 
}
\description{
  A collection of 41 different sequences of colours,
  each sequence having a uniform perceptual contrast over its whole
  range. These sequences make very good colour maps which avoid
  introducing artefacts when displaying image data.
}
\usage{data(Kovesi)}
\format{
  A \code{\link{hyperframe}} with the following columns:

  \tabular{ll}{
    \code{linear} \tab Logical: whether the sequence is linear. \cr
    \code{diverging} \tab Logical: whether the sequence is diverging. \cr
    \code{rainbow} \tab Logical: whether the sequence is a rainbow. \cr
    \code{cyclic} \tab Logical: whether the sequence is cyclic. \cr
    \code{isoluminant} \tab Logical: whether the sequence is isoluminant. \cr
    \code{ternary} \tab Logical: whether the sequence is ternary. \cr
    \code{colsig} \tab Character: colour signature (see Details) \cr
    \code{l1}, \code{l2} \tab Numeric: lightness parameters \cr
    \code{chro} \tab Numeric: average chroma (percent) \cr
    \code{n} \tab Numeric: length of colour sequence \cr
    \code{cycsh} \tab Numeric: cyclic shift (percent) \cr
    \code{values} \tab: Character: the colour values.
  }
}
\details{
  Kovesi (2014, 2015) presented a collection of colour sequences
  that have uniform perceptual contrast over their whole range.

  The dataset \code{Kovesi} provides these data. It is a
  \code{hyperframe} with 41 rows, in which each row provides information
  about one colour sequence.

  Additional information in each row specifies whether the
  colour sequence is \sQuote{linear}, \sQuote{diverging},
  \sQuote{rainbow}, \sQuote{cyclic}, \sQuote{isoluminant}
  and/or \sQuote{ternary} as defined by Kovesi (2014, 2015).

  The \sQuote{colour signature} is a string composed of letters
  representing the successive hues, using the following code:
  \tabular{ll}{
    r \tab red \cr
    g \tab green\cr
    b \tab blue \cr
    c \tab cyan \cr
    m \tab magenta \cr
    y \tab yellow \cr
    o \tab orange \cr
    v \tab violet \cr
    k \tab black \cr
    w \tab white \cr
    j \tab grey (j rhymes with grey)
  }
  For example \code{kryw} is the sequence from black to red to yellow to
  white.

  The column \code{values} contains the colour data themselves.
  The \code{i}th colour sequence is \code{Kovesi$values[[i]]},
  a character vector of length 256.
}
\source{
  Dr Peter Kovesi, Centre for Exploration Targeting,
  University of Western Australia.
}
\references{
  Kovesi, P. (2014)
  Website \emph{CET Uniform Perceptual Contrast Colour Maps}
  \url{www.peterkovesi.com/projects/colourmaps/}

  Kovesi, P. (2015)
  Designing colour maps with uniform perceptual contrast.
  Manuscript submitted for publication.
}
\examples{
  Kovesi
  LinearBMW <- Kovesi$values[[28]]
  plot(colourmap(LinearBMW, range=c(0,1)))

  ## The following would be suitable for spatstat.options(image.colfun)
  BMWfun <- function(n) { interp.colours(LinearBMW, n) }
}
\keyword{datasets}
