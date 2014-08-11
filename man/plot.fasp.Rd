\name{plot.fasp}
\alias{plot.fasp}
\title{Plot a Function Array}
\description{
   Plots an array of summary functions, usually associated with a
   point pattern, stored in an object of class \code{"fasp"}.
   A method for \code{plot}.
}

\usage{
   \method{plot}{fasp}(x,formule=NULL, \dots,
                        subset=NULL, title=NULL, banner=TRUE,
                        samex=FALSE, samey=FALSE,
                        mar.panel=NULL,
                        outerlabels=TRUE, cex.outerlabels=1.25,
                        legend=FALSE)
}

\arguments{
  \item{x}{An object of class \code{"fasp"} representing a
    function array.
  }
  \item{formule}{
    A formula or list of formulae indicating what
    variables are to be plotted against what variable. Each formula is
    either an R language formula object, or a string that can be parsed
    as a formula. If \code{formule} is a list, its \eqn{k^{th}}{k-th} component
    should be applicable to the \eqn{(i,j)^{th}}{(i,j)-th}
    plot where \code{x$which[i,j]=k}.  If the formula is left
    as \code{NULL}, then \code{plot.fasp} attempts to use the component
    \code{default.formula} of \code{x}.  If that component is NULL
    as well, it gives up.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.fv}} to control 
    the individual plot panels. 
  }
  \item{subset}{
    A logical vector, or a vector of indices, or an
    expression or a character string, or a \bold{list} of such,
    indicating a subset of the data to be included in each plot.
    If \code{subset} is a list, its \eqn{k^{th}}{k-th} component
    should be applicable to the \eqn{(i,j)^{th}}{(i,j)-th} plot
    where \code{x$which[i,j]=k}.
  }
  \item{title}{
    Overall title for the plot.
  }
  \item{banner}{
    Logical. If \code{TRUE}, the overall title is plotted.
    If \code{FALSE}, the overall title is not plotted
    and no space is allocated for it.
  }
  \item{samex,samey}{
    Logical values indicating whether all individual plot panels should have the
    same x axis limits and the same y axis limits, respectively.
    This makes it easier to compare the plots.
  }
  \item{mar.panel}{
    Vector of length 4 giving the value of the
    graphics parameter \code{mar} controlling the size of plot margins
    for each individual plot panel. See \code{\link{par}}.
  }
  \item{outerlabels}{Logical.
    If \code{TRUE}, the row and column names of the array of functions
    are plotted in the margins of the array of plot panels.
    If \code{FALSE}, each individual plot panel is labelled by its
    row and column name.
  }
  \item{cex.outerlabels}{
    Character expansion factor for row and column labels of array.
  }
  \item{legend}{
    Logical flag determining whether to plot a legend in each panel.
  }
}

\details{
  An object of class \code{"fasp"} represents
  an array of summary functions, usually associated with a point
  pattern. See \code{\link{fasp.object}} for details.
  Such an object is created, for example, by \code{\link{alltypes}}.

  The function \code{plot.fasp} is
  a method for \code{plot}.  It calls \code{\link{plot.fv}} to plot the
  individual panels.

  For information about the interpretation of the
  arguments \code{formule} and \code{subset},
  see \code{\link{plot.fv}}.

  Arguments that are often passed through \code{...} include
  \code{col} to control the colours of the different lines in a panel,
  and \code{lty} and \code{lwd} to control the line type and line width
  of the different lines in a panel. The argument \code{shade}
  can also be used to display confidence intervals or significance bands
  as filled grey shading. See \code{\link{plot.fv}}.
  
  The argument \code{title}, if present, will determine the
  overall title of the plot. If it is absent, it defaults to \code{x$title}.
  Titles for the individual plot panels will be taken from
  \code{x$titles}.
}

\value{None.}

\section{Warnings}{
  (Each component of) the \code{subset} argument may be a
  logical vector (of the same length as the vectors of data which
  are extracted from \code{x}), or a vector of indices, or an
  \bold{expression} such as \code{expression(r<=0.2)}, or a text string,
  such as \code{"r<=0.2"}.

  Attempting a syntax such as \code{subset = r<=0.2} (without
  wrapping \code{r<=0.2} either in quote marks or in \code{expression()})
  will cause this function to fall over.

  Variables referred to in any formula must exist in the data frames
  stored in \code{x}.  What the names of these variables are will
  of course depend upon the nature of \code{x}.
}

\seealso{
  \code{\link{alltypes}},
  \code{\link{plot.fv}},
  \code{\link{fasp.object}}
}

\examples{
   \dontrun{
   # Bramble Canes data.
   data(bramblecanes)

   X.G <- alltypes(bramblecanes,"G",dataname="Bramblecanes",verb=TRUE)
   plot(X.G)
   plot(X.G,subset="r<=0.2")
   plot(X.G,formule=asin(sqrt(cbind(km,theo))) ~ asin(sqrt(theo)))
   plot(X.G,fo=cbind(km,theo) - theo~r,subset="r<=0.2")

   # Simulated data.
   pp <- runifpoint(350, owin(c(0,1),c(0,1)))
   pp <- pp \%mark\% factor(c(rep(1,50),rep(2,100),rep(3,200)))
   X.K <- alltypes(pp,"K",verb=TRUE,dataname="Fake Data")
   plot(X.K,fo=cbind(border,theo)~theo,subset="theo<=0.75")
   }
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  }
\keyword{spatial}
\keyword{hplot}
