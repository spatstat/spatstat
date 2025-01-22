\name{latest.changes}
\alias{latest.changes}
\title{
  List Recent Significant Changes to a Function
}
\description{
  List all the changes to a particular function in a package,
  starting from a certain date or version of the package.
  The default is to list changes in the latest
  version of the \pkg{spatstat} family of packages.  
}
\usage{
  latest.changes(x,
          sinceversion = NULL, sincedate = NULL,
          package = spatstat.family(), show = TRUE)
}
\arguments{
  \item{x}{
    Character string giving the name
    of the function of interest, or a search pattern to be matched.
    A character vector is permitted.
  }
  \item{sinceversion}{
    Earliest version of \code{package}
    for which changes should be listed.
    A character string.
    Default is the most recent version.
  }
  \item{sincedate}{
    Earliest release date of \code{package}
    for which changes should be listed.
    A character string or a date-time object.
    Default is the date of the most recent version.
  }
  \item{package}{
    The name of the package (or packages) for which changes are to be listed.
    A character string or character vector.
  }
  \item{show}{
    Logical value indicating whether to display the table of changes
    on the terminal.
  }
}
\details{
  Details of changes are extracted from the NEWS file
  of the specified \code{package} under the heading
  \sQuote{Significant User-Visible Changes}.
  All entries for which the first line matches \code{x} are selected.
  The table of changes
  is displayed (if \code{show=TRUE}) and returned invisibly.
  
  The argument \code{sinceversion} should be a character string
  like \code{"1.2-3"}. The default is the version string of the most
  recent version.

  The argument \code{sincedate} should be a 
  character string like \code{"2015-05-27"}, or a date-time object.
  The default is the date of the most recent version.
  
  If \code{sinceversion="all"} or \code{sincedate="all"}
  then all recorded changes will be listed. 
  
  The special options \code{sinceversion="book"} and \code{sincedate="book"}
  are interpreted to mean \code{sincedate="2015-06-05"},
  which gives all changes reported after publication of
  the book by Baddeley, Rubak and Turner (2015).

  By default, changes in the \emph{extension} packages
  \pkg{spatstat.local}, \pkg{spatstat.Knet}, \pkg{spatstat.gui}
  are \emph{not} reported. To include these changes as well,
  set \code{package=spatstat.family(TRUE,TRUE)}.
}
\section{Pattern matching}{
  The first argument \code{x} is interpreted as a pattern to be matched
  using \code{\link[base]{grepl}}.
  
  For example if \code{x = "ppp"}, then this will match any news items
  that are about \code{ppp}, but will also match any news items
  about \code{plot.ppp}, \code{print.ppp}, \code{ppplist},
  \code{pppmatch}, \code{[.ppp} and so on.

  The pattern \code{x = "^ppp"} will match only those strings
  which begin with \code{"ppp"}, such as \code{ppp}, \code{ppplist},
  \code{pppmatch}. The symbol \code{^} means the start of the string.

  The pattern \code{x = "^ppp$"} will match only the string
  \code{"ppp"}. The symbol \code{$} means the end of the string.

  The pattern \code{x = "\\.ppp$"} will match only those strings which
  end with \code{".ppp"}, such as \code{plot.ppp}, \code{print.ppp},
  \code{[.ppp}. The symbol \code{.} has a special meaning, but here we
  want  the character \code{"."},
  so we use \code{\\} to escape the special
  meaning, so that \code{\\.} means the character \code{"."},
  and \code{"\\.ppp$"} means any string ending in \code{".ppp"}.
}
\value{
  (Invisibly) a data frame, belonging to the class \code{"changetable"},
  which has a \code{print} method.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{bugfixes}},
  \code{\link{latest.news}},
  \code{\link[utils]{news}}.
}
\references{
  \baddrubaturnbook
}
\examples{
   latest.changes("plot.symbolmap")
}
\keyword{documentation}
