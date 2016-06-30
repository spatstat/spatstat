\name{run.simplepanel}
\alias{clear.simplepanel}
\alias{redraw.simplepanel}
\alias{run.simplepanel}
\title{
  Run Point-and-Click Interface
}
\description{
  Execute various operations in a simple
  point-and-click user interface.
}
\usage{
run.simplepanel(P, popup=TRUE, verbose = FALSE)
clear.simplepanel(P)
redraw.simplepanel(P, verbose = FALSE)
}
\arguments{
  \item{P}{
    An interaction panel (object of class \code{"simplepanel"},
    created by \code{\link{simplepanel}} or \code{\link{grow.simplepanel}}).
  }
  \item{popup}{
    Logical. If \code{popup=TRUE} (the default),
    the panel will be displayed in a
    new popup window. If \code{popup=FALSE},
    the panel will be displayed on the current graphics
    window if it already exists, and on a new window otherwise.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, debugging information will be printed.
  }
}
\details{
  These commands enable the user to run a simple, robust,
  point-and-click interface to any \R code.
  The interface is implemented
  using only the basic graphics package in \R.

  The argument \code{P} is an object of class \code{"simplepanel"},
  created by \code{\link{simplepanel}} or \code{\link{grow.simplepanel}},  
  which specifies the graphics to be displayed and the actions to be performed
  when the user interacts with the panel.

  The command \code{run.simplepanel(P)} activates the panel:
  the display is initialised and the graphics system waits for the
  user to click the panel. 
  While the panel is active, the user can only interact with the panel;
  the \R command line interface and the \R GUI cannot be used.
  When the panel terminates (typically because the user
  clicked a button labelled Exit), control returns to the \R command line
  interface and the \R GUI. 

  The command \code{clear.simplepanel(P)} clears all the display
  elements in the panel, resulting in a blank display except for the
  title of the panel.

  The command \code{redraw.simplepanel(P)} redraws all the buttons
  of the panel, according to the \code{redraw} functions contained
  in the panel.

  If \code{popup=TRUE} (the default), \code{run.simplepanel} begins by
  calling \code{\link[grDevices]{dev.new}} so that a new popup window
  is created; this window is closed using
  \code{\link[grDevices]{dev.off}} when \code{run.simplepanel} terminates.
  If \code{popup=FALSE}, the panel will be displayed on the current graphics
  window if it already exists, and on a new window otherwise;
  this window is not closed when \code{run.simplepanel} terminates.

  For more sophisticated control of the graphics focus (for example, to
  use the panel to control the display on another window),
  initialise the graphics devices yourself using
  \code{\link[grDevices]{dev.new}} or similar commands; save these devices
  in the shared environment \code{env} of the panel \code{P};
  and write the click/redraw functions of \code{P} in such a way that
  they access these devices using \code{\link[grDevices]{dev.set}}.
  Then use \code{run.simplepanel} with \code{popup=FALSE}.
}
\value{
  The return value of \code{run.simplepanel(P)} is the value returned
  by the \code{exit} function of \code{P}. See \code{\link{simplepanel}}.

  The functions \code{clear.simplepanel} and \code{redraw.simplepanel}
  return \code{NULL}.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{simplepanel}}
}
\examples{
  if(interactive()) {
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

  run.simplepanel(P)
  }
}
\keyword{iplot}
\keyword{utilities}
