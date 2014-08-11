\name{simplepanel}
\alias{simplepanel}
\alias{grow.simplepanel}
\title{Simple Point-and-Click Interface Panels}
\description{
  These functions enable the user to create a simple, robust, 
  point-and-click interface to any \R code.
}
\usage{
   simplepanel(title, B, boxes, clicks,
      redraws=NULL, exit = NULL, env)

   grow.simplepanel(P, side = c("right", "left", "top", "bottom"),
      len = NULL, new.clicks, new.redraws=NULL, \dots, aspect)
}
\arguments{
  \item{title}{
    Character string giving the title of the interface panel.
  }
  \item{B}{
    Bounding box of the panel coordinates.
    A rectangular window (object of class \code{"owin"})
  }
  \item{boxes}{
    A list of rectangular windows (objects of class \code{"owin"})
    specifying the placement of the buttons
    and other interactive components of the panel.
  }
  \item{clicks}{
    A list of \R functions, of the same length as \code{boxes},
    specifying the operations to be performed when each button
    is clicked. See Details.
  }
  \item{redraws}{
    Optional list of \R functions, of the same length as \code{boxes},
    specifying how to redraw each button. See Details.
  }
  \item{exit}{
    An \R function specifying actions to be taken when the
    interactive panel terminates.
  }
  \item{env}{
    An \code{environment} that will be passed as an argument
    to all the functions in \code{clicks}, \code{redraws} and
    \code{exit}.
  }
  \item{P}{
    An existing interaction panel (object of class \code{"simplepanel"}).
  }
  \item{side}{
    Character string identifying which side of the panel \code{P}
    should be grown to accommodate the new buttons.
  }
  \item{len}{
    Optional. Thickness of the new panel area that should be grown
    to accommodate the new buttons. A single number in the same units
    as the coordinate system of \code{P}.
  }
  \item{new.clicks}{
    List of \R functions defining the operations to be performed
    when each of the new buttons is clicked.
  }
  \item{new.redraws}{
    Optional. List of \R functions, of the same length as
    \code{new.clicks}, defining how to redraw each of the new buttons.
  }
  \item{\dots}{
    Arguments passed to \code{\link{layout.boxes}} to determine the
    layout of the new buttons.
  }
  \item{aspect}{
    Optional. Aspect ratio (height/width) of the new buttons.
  }
}
\details{
  These functions enable the user to create a simple, robust,
  point-and-click interface to any \R code.

  The functions \code{simplepanel} and \code{grow.simplepanel}
  create an object of class \code{"simplepanel"}. Such an object defines
  the graphics to be displayed and the actions to be performed
  when the user interacts with the panel.

  The panel is activated by calling \code{\link{run.simplepanel}}.

  The function \code{simplepanel} creates a panel object
  from basic data.
  The function \code{grow.simplepanel} modifies an existing
  panel object \code{P} by growing an additional row or column
  of buttons. 

  For \code{simplepanel},
  \itemize{
    \item 
    The spatial layout of the panel is determined by the rectangles
    \code{B} and \code{boxes}.
    \item 
    The argument \code{clicks} must be a list of functions
    specifying the action to be taken when each button is clicked.
    The list entries should have names.
    Each function should be of the form \code{function(env, xy)} where
    \code{env} is an \code{environment} that may contain shared data,
    and \code{xy} gives the coordinates of the mouse click, in the format
    \code{list(x, y)}.
    The function returns \code{TRUE} if the
    panel should continue running, and \code{FALSE} if the panel
    should terminate.
    \item 
    The argument \code{redraws}, if given, must be a list of functions
    specifying the action to be taken when each button is to be redrawn.
    Each function should be of the form \code{function(button, name, env)} where
    \code{button} is a rectangle specifying the location of the button
    in the current coordinate system; \code{name} is a character string
    giving the name of the button; and \code{env} is the
    \code{environment} that may contain shared data.
    The function returns \code{TRUE} if the
    panel should continue running, and \code{FALSE} if the panel
    should terminate. 
    If \code{redraws} is not given, the default action is to draw a pink
    rectangle showing the button position,
    draw the name of the button in the middle of this rectangle,
    and return \code{TRUE}. 
    \item
    The argument \code{exit}, if given, must be a function
    specifying the action to be taken when the panel terminates.
    (Termination occurs when one of the \code{clicks} functions
    returns \code{FALSE}).
    The \code{exit} function should be of the form \code{function(env)} where
    \code{env} is the \code{environment} that may contain shared data.
    Its return value will be used as the return value
    of \code{\link{run.simplepanel}}.
    \item
    The argument \code{env} should be an \R environment.
    The panel buttons will have access to this environment,
    and will be able to read and write data in it. This mechanism is used
    to exchange data between the panel and other \R code.
  }
  For \code{grow.simplepanel},
  \itemize{
    \item the spatial layout of the new boxes
    is determined by the arguments \code{side}, \code{len},
    \code{aspect} and by the additional \code{\dots} arguments passed to
    \code{\link{layout.boxes}}.
    \item the argument \code{new.clicks} 
    should have the same format as \code{clicks}.
    It implicitly specifies the number of new buttons to be added,
    and the actions to be performed when they are clicked.
    \item the optional argument \code{new.redraws}, if given,
    should have the same format as \code{redraws}.
    It specifies the actions to be performed when the
    new buttons are clicked.
  }
}
\value{
  An object of class \code{"simplepanel"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{run.simplepanel}},
  \code{\link{layout.boxes}}
}
\examples{
  # make boxes (alternatively use layout.boxes())
  Bminus <- square(1)
  Bvalue <- shift(Bminus, c(1.2, 0))
  Bplus <- shift(Bvalue, c(1.2, 0))
  Bdone <- shift(Bplus, c(1.2, 0))
  myboxes <- list(Bminus, Bvalue, Bplus, Bdone)
  myB <- do.call(boundingbox,myboxes)

  # make environment containing an integer count
  myenv <- new.env()
  assign("answer", 0, envir=myenv)

  # what to do when finished: return the count.
  myexit <- function(e) { return(get("answer", envir=e)) }

  # button clicks
  # decrement the count
  Cminus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans - 1, envir=e)
    return(TRUE)
  }
  # display the count (clicking does nothing)
  Cvalue <- function(...) { TRUE }
  # increment the count
  Cplus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans + 1, envir=e)
    return(TRUE)
  }
  # quit button
  Cdone <- function(e, xy) { return(FALSE) }

  myclicks <- list("-"=Cminus,
                   value=Cvalue,
                   "+"=Cplus,
                   done=Cdone)

  # redraw the button that displays the current value of the count
  Rvalue <- function(button, nam, e) {
     plot(button, add=TRUE)
     ans <- get("answer", envir=e)
     text(centroid.owin(button), labels=ans)
     return(TRUE)
  }

  # make the panel
  P <- simplepanel("Counter",
                   B=myB, boxes=myboxes,
                   clicks=myclicks,
                   redraws = list(NULL, Rvalue, NULL, NULL),
                   exit=myexit, env=myenv)
  P
}
\keyword{iplot}
\keyword{utilities}
