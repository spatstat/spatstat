\name{flu}
\alias{flu}
\docType{data}
\title{
  Influenza Virus Proteins
}
\description{
  Replicated spatial point patterns giving the locations of two
  different virus proteins on the membranes of cells infected with
  influenza virus.
}
\usage{data(flu)}
\format{
  A \code{\link{hyperframe}} with 41 rows and four columns:
  \describe{
    \item{pattern}{
      List of spatial point patterns
      (objects of class \code{"ppp"})
      with points of two types, identifying the locations of 
      two different proteins on a membrane sheet.
    }
    \item{virustype}{
      Factor identifying whether the infecting virus was
      the wild type (\code{wt}) or mutant (\code{mut1}).
    }
    \item{stain}{
      Factor identifying whether the membrane sheet was stained
      for the proteins \emph{M2} and \emph{M1}
      (\code{stain="M2-M1"})
      or stained for the proteins \emph{M2} and \emph{HA}
      (\code{stain="M2-HA"}).
    }
    \item{frameid}{
      Integer. Serial number of the microscope frame
      in the original experiment. Frame identifier is not unique
      across different values of \code{virustype} and \code{stain}.
    }
  }
  The row names of the hyperframe can be used as succinct labels
  in plots.
}
\details{
  The data consist of 41 spatial point patterns, each
  giving the locations of two different virus proteins
  on the membranes of cells infected with influenza virus.
  
  Chen et al (2008) conducted the experiment and used spatial analysis
  to establish evidence for an interaction between the influenza virus
  proteins M1 and M2 that is important for the study of viral replication.

  Canine kidney cells were infected with human influenza, Udorn strain,
  either the wild type or a mutant which encodes a defective M2 protein.
  At twelve hours post-infection, membrane sheets were prepared
  and stained for viral proteins, using two antibodies conjugated to
  gold particles of two sizes (6 nanometre and 12 nanometre diameter)
  enabling localisation of two different proteins on each sheet.
  The 6 nm particles were stained for M2 (ion channel protein),
  while the 12 nm particles were stained either for M1 (matrix protein)
  or for HA (hemagglutinin). Membrane sheets were visualised in
  electron microscopy.
  
  Experimental technique and spatial analysis of the membranes
  stained for M2 and M1 is reported in Chen et al (2008).
  Analysis of the membranes stained for
  M2 and HA is reported in Rossman et al (2010).
  The M2-HA data shows a stronger association
  between the two proteins which has also been observed biochemically
  and functionally (Rossman et al, 2010).

  The dataset \code{flu} is a \code{\link{hyperframe}}
  with one row for each membrane sheet. The column named \code{pattern}
  contains the spatial point patterns of gold particle locations,
  with two types of points (either \code{M1} and \code{M2} or
  \code{HA} and \code{M2}). The column named \code{virustype}
  is a factor identifying the virus: either wild type \code{wt}
  or mutant \code{mut1}. The column named \code{stain} is a factor
  identifying whether the membrane was stained for
  M1 and M2 (\code{stain="M2-M1"}) or stained for HA and M2
  (\code{stain="M2-HA"}).
  The row names of the hyperframe are a succinct summary of
  the experimental conditions and can be used as labels
  in plots. See the Examples.
}
\source{
  Data generously provided by Dr G.P. Leser and Dr R.A. Lamb.
  Please cite Chen et al (2008) in any use of these data.
}
\references{
  Chen, B.J., Leser, G.P., Jackson, D. and Lamb, R.A. (2008)
  The influenza virus M2 protein cytoplasmic tail interacts with the M1 protein
  and influences virus assembly at the site of virus budding.
  \emph{Journal of Virology} \bold{82}, 10059--10070.

  Rossman, J.S., Jing, X.H.,  Leser, G.P. and Lamb, R.A. (2010)
  Influenza virus M2 protein mediates ESCRT-independent membrane scission
  \emph{Cell} \bold{142}, 902--913.
}
\examples{
data(flu)
flu
Y <- flu$pattern[10]
Y <- flu[10, 1, drop=TRUE]
wildM1 <- with(flu, virustype == "wt" & stain == "M2-M1")
plot(flu[wildM1, 1, drop=TRUE], 
     main=c("flu data", "wild type virus, M2-M1 stain"),
     pch=c(3,16), cex=0.4, cols=2:3)
}
\keyword{datasets}
