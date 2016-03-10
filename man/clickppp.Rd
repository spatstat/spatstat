\name{clickppp}
\alias{clickppp}
\title{Interactively Add Points}
\description{
  Allows the user to create a point pattern by
  point-and-click in the display.
}
\usage{
  clickppp(n=NULL, win=square(1), types=NULL, \dots, add=FALSE,
    main=NULL, hook=NULL)
}
\arguments{
  \item{n}{
    Number of points to be added (if this is predetermined).
  }
  \item{win}{
    Window in which to create the point pattern.
    An object of class \code{"owin"}.
  }
  \item{types}{
    Vector of types, when creating a multitype point pattern.
  }
  \item{\dots}{
    Optional extra arguments to be passed to \code{\link[graphics]{locator}}
    to control the display.
  }
  \item{add}{
    Logical value indicating whether to create a new plot
    (\code{add=FALSE}) or draw over the existing plot (\code{add=TRUE}).
  }
  \item{main}{
    Main heading for plot.
  }
  \item{hook}{For internal use only. Do not use this argument.}
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\details{
  This function allows the user to create a point pattern
  by interactively clicking on the screen display.

  First the window \code{win} is plotted on the current screen device.
  Then the user is prompted to point the mouse at any desired locations
  and click the left mouse button to add each point.
  Interactive input stops after \code{n} clicks (if \code{n} was given)
  or when the middle mouse button is pressed.

  The return value is a point pattern
  containing the locations of all the clicked points
  inside the original window \code{win},
  provided that all of the clicked locations were
  inside this window. Otherwise, the window is expanded to a box
  large enough to contain all the points (as well as containing
  the original window).

  If the argument \code{types} is given, then a multitype point pattern
  will be created. The user is
  prompted to input the locations of points of type \code{type[i]},
  for each successive index \code{i}. (If the argument \code{n} was
  given, there will be \code{n} points of \emph{each} type.)
  The return value is a multitype point pattern.
  
  This function uses the \R{} command \code{\link[graphics]{locator}} to
  input the mouse clicks. It only works on screen devices such as
  \sQuote{X11}, \sQuote{windows} and \sQuote{quartz}. Arguments that can be
  passed to \code{\link[graphics]{locator}} through \code{\dots} include
  \code{pch} (plotting character), \code{cex} (character expansion
  factor) and \code{col} (colour). See \code{\link[graphics]{locator}}
  and \code{\link{par}}.
}
\seealso{
  \code{\link{identify.ppp}},
  \code{\link[graphics]{locator}},
  \code{\link{clickpoly}},
  \code{\link{clickbox}},
  \code{\link{clickdist}}
}
\author{Original by Dominic Schuhmacher.
  Adapted by Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{iplot}
