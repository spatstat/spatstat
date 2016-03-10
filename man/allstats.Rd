\name{allstats}
\alias{allstats}
\title{Calculate four standard summary functions of a point pattern.}
\description{
  Calculates the \eqn{F}, \eqn{G}, \eqn{J}, and \eqn{K}
  summary functions for an unmarked point pattern.
  Returns them as a function array (of class \code{"fasp"}, see
  \code{\link{fasp.object}}).
}
\usage{
  allstats(pp, \dots, dataname=NULL, verb=FALSE)
}
\arguments{
  \item{pp}{The observed point pattern, for which summary function
    estimates are required.  An object of class \code{"ppp"}.
    It must not be marked.
  }
  \item{\dots}{
    Optional arguments passed to the summary functions
    \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}}
    and \code{\link{Kest}}.
  }
  \item{dataname}{A character string giving an optional (alternative)
    name for the point pattern.
  }
  \item{verb}{A logical value meaning ``verbose''. If \code{TRUE},
    progress reports are printed during calculation.
  }
}
\details{
  This computes four standard summary statistics for a
  point pattern: the empty space function \eqn{F(r)},
  nearest neighbour distance distribution function \eqn{G(r)},
  van Lieshout-Baddeley function \eqn{J(r)}
  and Ripley's function \eqn{K(r)}.
  The real work is done by 
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}} and
  \code{\link{Kest}} respectively. Consult the help files for these functions
  for further information about the statistical interpretation
  of \eqn{F}, \eqn{G}, \eqn{J} and \eqn{K}.

  If \code{verb} is \code{TRUE}, then ``progress reports''
  (just indications of completion) are printed out when the
  calculations are finished for each of the four function types.

  The overall title of the array of four functions
  (for plotting by \code{\link{plot.fasp}})
  will be formed from the argument \code{dataname}.
  If this is not given, it defaults to the expression
  for \code{pp} given in the call to \code{allstats}.
}
\value{
  A list of length 4 containing the \eqn{F}, \eqn{G}, \eqn{J} and
  \eqn{K} functions respectively. 

  The list can be plotted directly using \code{plot} (which dispatches to
  \code{\link{plot.solist}}). 

  Each list entry retains the format
  of the output of the relevant estimating routine
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}} or
  \code{\link{Kest}}. Thus each entry in the list is
  a function value table (object of class \code{"fv"},
  see \code{\link{fv.object}}).

  The default formulae for plotting these functions are 
  \code{cbind(km,theo) ~ r} for F, G, and J, and
  \code{cbind(trans,theo) ~ r} for K.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
         \code{\link{plot.solist}},
         \code{\link{plot.fv}},
         \code{\link{fv.object}},
         \code{\link{Fest}},
         \code{\link{Gest}},
         \code{\link{Jest}},
         \code{\link{Kest}}
}

\examples{
	data(swedishpines)
        a <- allstats(swedishpines,dataname="Swedish Pines")
        \dontrun{
        plot(a)
        plot(a, subset=list("r<=15","r<=15","r<=15","r<=50"))
        }
}
\keyword{spatial}
\keyword{nonparametric}
