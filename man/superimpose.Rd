\name{superimpose}  %DontDeclareMethods
\alias{superimpose}
\alias{superimpose.ppp}
\alias{superimpose.splitppp}
\alias{superimpose.ppplist}
\alias{superimpose.psp}
\alias{superimpose.default}
\title{Superimpose Several Geometric Patterns}
\description{
  Superimpose any number of point patterns or line segment patterns.
}
\usage{
  superimpose(\dots)

  \method{superimpose}{ppp}(\dots, W=NULL, check=TRUE)

  \method{superimpose}{psp}(\dots, W=NULL, check=TRUE)

  \method{superimpose}{splitppp}(\dots, W=NULL, check=TRUE)

  \method{superimpose}{ppplist}(\dots, W=NULL, check=TRUE)

  \method{superimpose}{default}(\dots)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which represents either a point
    pattern or a line segment pattern or a list of point patterns.
  }
  \item{W}{
    Optional. Data determining the window for the resulting pattern.
    Either a window (object of class \code{"owin"}, or something
    acceptable to \code{\link{as.owin}}), or a function
    which returns a window, or one of the strings
    \code{"convex"}, \code{"rectangle"}, \code{"bbox"} or \code{"none"}.
  }
  \item{check}{
    Logical value (passed to \code{\link{ppp}} or \code{\link{psp}}
    as appropriate) determining whether to check the geometrical
    validity of the resulting pattern.
  }
}
\value{
  For \code{superimpose.ppp}, a point pattern (object of class \code{"ppp"}).
  For \code{superimpose.default}, either a point pattern
  (object of class \code{"ppp"}) or a \code{list(x,y)}.
  For \code{superimpose.psp}, a line segment pattern (object of class
  \code{"psp"}).
}
\details{
  This function is used to superimpose several geometric patterns
  of the same kind, producing a single pattern of the same kind.
  
  The function \code{superimpose} is generic, with methods
  for the class \code{ppp} of point patterns, the class \code{psp}
  of line segment patterns, and a default method. There is also
  a method for \code{lpp}, described separately in
  \code{\link{superimpose.lpp}}.

  The dispatch to a method is initially determined
  by the class of the \emph{first} argument in \code{\dots}.
  \itemize{
    \item \code{default}:  If the first argument is 
    \emph{not} an object of class \code{ppp} or \code{psp}, then the
    default method \code{superimpose.default} is executed.
    This checks the class of all arguments, and dispatches to the
    appropriate method. Arguments of class \code{ppplist} can be handled.
    \item \code{ppp}:
    If the first \code{\dots} argument is an object of
    class \code{ppp} then the method \code{superimpose.ppp}
    is executed. All arguments in \code{\dots}
    must be either \code{ppp} objects or lists
    with components \code{x} and \code{y}. The result will
    be an object of class \code{ppp}.
    \item psp:
    If the first \code{\dots} argument is an object of
    class \code{psp} then the \code{psp} method is dispatched and all
    \code{\dots} arguments must be \code{psp} objects.
    The result is a \code{psp} object.
  }

  The patterns are \emph{not} required to have the same window
  of observation. 

  The window for the superimposed pattern is controlled
  by the argument \code{W}.
  \itemize{
    \item
    If \code{W} is a window (object of class \code{"W"}
    or something acceptable to \code{\link{as.owin}})
    then this determines the window for the superimposed pattern.
    \item
    If \code{W} is \code{NULL}, or the character string \code{"none"},
    then windows are extracted from the geometric patterns,
    as follows. 
    For \code{superimpose.psp}, all arguments \code{\dots}
    are line segment patterns (objects of class \code{"psp"});
    their observation windows are extracted; the union of these
    windows is computed; and this union is taken to be the
    window for the superimposed pattern.
    For \code{superimpose.ppp} and \code{superimpose.default},
    the arguments \code{\dots} are inspected, and 
    any arguments which are point patterns (objects of class
    \code{"ppp"}) are selected; their observation windows
    are extracted, and the union of these windows
    is taken to be the window for the superimposed point pattern.
    For \code{superimpose.default} if none of the arguments
    is of class \code{"ppp"} then no window is computed
    and the result of \code{superimpose} is a \code{list(x,y)}.
    \item
    If \code{W} is one of the strings
    \code{"convex"}, \code{"rectangle"} or \code{"bbox"}
    then a window for the superimposed pattern
    is computed from the coordinates of the
    points or the line segments as follows.
    \describe{
      \item{\code{"bbox"}:}{the bounding box of the
      points or line segments (see \code{\link{bounding.box.xy}});}
      \item{\code{"convex"}:}{the Ripley-Rasson estimator
      of a convex window (see \code{\link{ripras}});}
      \item{\code{"rectangle"}:}{the Ripley-Rasson estimator
      of a rectangular window
      (using \code{\link{ripras}} with argument
      \code{shape="rectangle"}).}
    }
    \item
    If \code{W} is a function,
    then this function is used to compute
    a window for the superimposed pattern
    from the coordinates of the
    points or the line segments.
    The function should accept input of the form \code{list(x,y)}
    and is expected to return an object of class \code{"owin"}.
    Examples of such functions are
    \code{\link{ripras}} and \code{\link{bounding.box.xy}}.
  }
  The arguments \code{\dots} may be \emph{marked} patterns.
  The marks of each component pattern must have the same format.
  Numeric and character marks may be ``mixed''.  If there is such
  mixing then the numeric marks are coerced to character in the
  combining process. If the mark structures are all data frames,
  then these data frames must have the same number of columns and
  identical column names.

  If the arguments \code{\dots} are given in the form \code{name=value},
  then the \code{name}s will be used as an extra column of marks
  attached to the elements of the corresponding patterns.
}

\seealso{
  \code{\link{superimpose.lpp}},
  \code{\link{concatxy}}, \code{\link{quadscheme}}.
}
\examples{
  # superimposing point patterns
  p1  <- runifrect(30)
  p2  <- runifrect(42)
  s1  <- superimpose(p1,p2) # Unmarked pattern.
  p3  <- list(x=rnorm(20),y=rnorm(20))
  s2  <- superimpose(p3,p2,p1) # Default method gets called.
  s2a <- superimpose(p1,p2,p3) # Same as s2 except for order of points.
  s3  <- superimpose(clyde=p1,irving=p2) # Marked pattern; marks a factor
                                         # with levels "clyde" and "irving";
                                         # warning given.
  marks(p1) <- factor(sample(LETTERS[1:3],30,TRUE))
  marks(p2) <- factor(sample(LETTERS[1:3],42,TRUE))
  s5  <- superimpose(clyde=p1,irving=p2) # Marked pattern with extra column
  marks(p2) <- data.frame(a=marks(p2),b=runif(42))
  s6  <- try(superimpose(p1,p2)) # Gives an error.
  marks(p1) <- data.frame(a=marks(p1),b=1:30)
  s7  <- superimpose(p1,p2) # O.K.

  # how to make a 2-type point pattern with types "a" and "b"
  u <- superimpose(a = rpoispp(10), b = rpoispp(20))

  # how to make a 2-type point pattern with types 1 and 2
  u <- superimpose("1" = rpoispp(10), "2" = rpoispp(20))
 
  # superimposing line segment patterns
  X <- rpoisline(10)
  Y <- as.psp(matrix(runif(40), 10, 4), window=owin())
  Z <- superimpose(X, Y)

  # being unreasonable
  \dontrun{
   if(FALSE) {
    crud <- try(superimpose(p1,p2,X,Y)) # Gives an error, of course!
   }
  }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
