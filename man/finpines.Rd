\name{finpines}
\alias{finpines}
\alias{finpines.extra}
\docType{data}
\title{
  Pine saplings in Finland.
}
\description{
  The data record the locations of 126 pine saplings
  in a Finnish forest, their heights and their diameters.

  The dataset \code{finpines} is a marked point pattern
  containing the locations of the saplings marked by their heights
  and their diameters.

  Sapling locations are given in metres (to six significant digits);
  heights are in metres (rounded to the nearest 0.1 metre,
  except in one case to the nearest 0.05 metres);
  diameters are in centimetres (rounded to the nearest centimetre).
  
  The data were recorded by Professor Seppo Kellomaki, Faculty of
  Forestry, University of Joensuu, Finland,
  and subsequently massaged by Professor Antti Penttinen,
  Department of Statistics, University of Jyv\"askyl\"a, Finland.
  
  Originally the point locations were observed in polar coordinates with
  rather poor angular precision. Hence the coordinates are imprecise for large
  radius because of rounding errors: indeed the alignments can be observed by
  eye.
  
  The data were manipulated by Prof Penttinen by making small angular
  perturbations at random. After this transformation, the original data
  (in a circular plot) were clipped to a square window, for convenience.

  Professor Penttinen emphasises that the data were intended only
  for initial experimentation. They have some strange features.
  For example, if the height is less than 1.3 metres then the diameter
  can be uncertain. Also there are some very close pairs of points.
  Some pairs of trees (namely (58,59), (78,79), (96,97) and (102,103))
  violate the requirement that the interpoint distance should be
  greater than half the sum of their diameters.

  These data have subsequently been analysed by Van Lieshout (2004).
} 
\format{
  Object of class \code{"ppp"} 
  representing the point pattern of sapling locations marked by
  their heights and diameters.
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(finpines)}
\examples{
    data(finpines)
    plot(unmark(finpines), main="Finnish pines: locations")
    plot(finpines, which.marks="height", main="heights")
    plot(finpines, which.marks="diameter", main="diameters")
    plot(finpines, which.marks="diameter", 
              main="diameters to scale", markscale=1/200)
}
\source{Prof Antti Penttinen}
\references{
  Van Lieshout, M.N.M. (2004)
  A J-function for marked point patterns.
  Research Report PNA-R0404, June 2004.
  Centrum voor Wiskunde en Informatica (CWI), Amsterdam, 2004.
}
\keyword{datasets}
\keyword{spatial}
