\name{rmpoint}
\alias{rmpoint}
\title{Generate N Random Multitype Points}
\description{
  Generate a random multitype point pattern
  with a fixed number of points, or a fixed number of points of each type.
}
\usage{
 rmpoint(n, f=1, fmax=NULL, win=unit.square(),
         types, ptypes,
         \dots, giveup=1000, verbose=FALSE,
         nsim=1, drop=TRUE)
}
\arguments{
  \item{n}{
    Number of marked points to generate.
    Either a single number specifying the total number of points,
    or a vector specifying the number of points of each type.
  }
  \item{f}{
    The probability density of the multitype points,
    usually un-normalised.
    Either a constant, a vector,
    a function \code{f(x,y,m, ...)}, a pixel image,
    a list of functions \code{f(x,y,...)}
    or a list of pixel images.
  }
  \item{fmax}{
    An upper bound on the values of \code{f}.
    If missing, this number will be estimated.
  }
  \item{win}{
    Window in which to simulate the pattern.
    Ignored if \code{f} is a pixel image or list of pixel images.
  }
  \item{types}{
    All the possible types for the multitype pattern. 
  }
  \item{ptypes}{
    Optional vector of probabilities for each type.
  }
  \item{\dots}{
    Arguments passed to \code{f} if it is a function.
  }
  \item{giveup}{
    Number of attempts in the rejection method after which
    the algorithm should stop trying to generate new points.
  }
  \item{verbose}{
    Flag indicating whether to report details of performance
    of the simulation algorithm.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1},
  or a list of point patterns if \code{nsim > 1}.
}
\details{
  This function generates random multitype point patterns
  consisting of a fixed number of points.
  
  Three different models are available:
  \describe{
    \item{I. Random location and type:}{
      If \code{n} is a single number and the argument \code{ptypes} is missing,
      then \code{n} independent, identically distributed
      random multitype points are generated.
      Their locations
      \code{(x[i],y[i])} and types \code{m[i]} have
      joint probability density proportional to \eqn{f(x,y,m)}.
    }
    \item{II. Random type, and random location given type:}{
      If \code{n} is a single number and \code{ptypes} is given,
      then \code{n} independent, identically distributed
      random multitype points are generated.
      Their types \code{m[i]} have probability distribution
      \code{ptypes}. Given the types, the locations \code{(x[i],y[i])}
      have conditional probability density proportional to
      \eqn{f(x,y,m)}.
    }
    \item{III. Fixed types, and random location given type:}{
      If \code{n} is a vector, then we generate \code{n[i]}
      independent, identically distributed random points of type
      \code{types[i]}. For points of type \eqn{m} the conditional probability
      density of location \eqn{(x,y)} is proportional to
      \eqn{f(x,y,m)}.
    }
  }
  Note that the density \code{f} is normalised in different ways
  in Model I and Models II and III. In Model I the normalised
  joint density is \eqn{g(x,y,m)=f(x,y,m)/Z} where
  \deqn{
    Z = \sum_m \int\int \lambda(x,y,m) {\rm d}x \, {\rm d}y
  }{
    Z = sum_[m] integral lambda(x,y,m) dx dy
  }
  while in Models II and III the normalised conditional density
  is \eqn{g(x,y\mid m) = f(x,y,m)/Z_m}{g(x,y|m) = f(x,y,m)/Z[m]}
  where
  \deqn{
    Z_m = \int\int \lambda(x,y,m)  {\rm d}x \, {\rm d}y.
  }{
    Z[m] = integral lambda(x,y,m) dx dy.
  }
  In Model I, the marginal distribution of types
  is \eqn{p_m = Z_m/Z}{p[m] = Z[m]/Z}.
  
  The unnormalised density \code{f} may be specified
  in any of the following ways.
  \describe{
    \item{single number:}{
      If \code{f} is a single number, the conditional density of
      location given type is uniform. That is, the points of each type
      are uniformly distributed.
      In Model I, the marginal distribution of types is also uniform
      (all possible types have equal probability).
    }
    \item{vector:}{
      If \code{f} is a numeric vector, the conditional density of
      location given type is uniform. That is, the points of each type
      are uniformly distributed.
      In Model I, the marginal distribution of types is
      proportional to the vector \code{f}. In Model II, the marginal
      distribution of types is \code{ptypes}, that is, the values in
      \code{f} are ignored.
      The argument \code{types} defaults to \code{names(f)},
      or if that is null, \code{1:length(f)}.
    }
    \item{function:}{
      If \code{f} is a function, it will be called in the form
      \code{f(x,y,m,\dots)} at spatial location \code{(x,y)}
      for points of type \code{m}.
      In Model I, the joint probability density of location and type is
      proportional to \code{f(x,y,m,\dots)}.
      In Models II and III, the conditional probability density of
      location \code{(x,y)} given type \code{m} is
      proportional to \code{f(x,y,m,\dots)}.
      The function \code{f} must work correctly with vectors \code{x},
      \code{y} and \code{m}, returning a vector of function values.
      (Note that \code{m} will be a factor
      with levels \code{types}.)
      The value \code{fmax} must be given and must be an upper bound on the 
      values of \code{f(x,y,m,\dots)} for all locations \code{(x, y)}
      inside the window \code{win} and all types \code{m}.
      The argument \code{types} must be given.
    }
    \item{list of functions:}{
      If \code{f} is a list of functions, then the functions will be
      called in the form \code{f[[i]](x,y,\dots)} at spatial
      location \code{(x,y)} for points of type \code{types[i]}.
      In Model I, the joint probability density of location and type is
      proportional to \code{f[[m]](x,y,\dots)}.
      In Models II and III, the conditional probability density of
      location \code{(x,y)} given type \code{m} is
      proportional to \code{f[[m]](x,y,\dots)}.
      The function \code{f[[i]]} must work correctly with vectors
      \code{x} and \code{y}, returning a vector of function values.
      The value \code{fmax} must be given and must be an upper bound on the 
      values of \code{f[[i]](x,y,\dots)} for all locations \code{(x, y)}
      inside the window \code{win}. 
      The argument \code{types} defaults to \code{names(f)},
      or if that is null, \code{1:length(f)}.
    }
    \item{pixel image:}{
      If \code{f} is a pixel image object of class \code{"im"}
      (see \code{\link{im.object}}), the unnormalised density at a location
      \code{(x,y)} for points of any type is equal to the pixel value
      of \code{f} for the pixel nearest to \code{(x,y)}.
      In Model I, the marginal distribution of types is uniform.
      The argument \code{win} is ignored;
      the window of the pixel image is used instead.
      The argument \code{types} must be given.
    }
    \item{list of pixel images:}{
      If \code{f} is a list of pixel images,
      then the image \code{f[[i]]} determines the density values
      of points of type \code{types[i]}. 
      The argument \code{win} is ignored;
      the window of the pixel image is used instead.
      The argument \code{types} defaults to \code{names(f)},
      or if that is null, \code{1:length(f)}.
    }
  }
  The implementation uses the rejection method.
  For Model I, \code{\link{rmpoispp}} is called repeatedly
  until \code{n} points have been generated.
  It gives up after \code{giveup} calls
  if there are still fewer than \code{n} points.
  For Model II, the types are first generated according to
  \code{ptypes}, then
  the locations of the points of each type
  are generated using \code{\link{rpoint}}.
  For Model III,   the locations of the points of each type
  are generated using \code{\link{rpoint}}.
}
\seealso{
\code{\link{ppp.object}},
\code{\link{owin.object}}
}
\examples{

abc <- c("a","b","c")

##### Model I

rmpoint(25, types=abc)
rmpoint(25, 1, types=abc)
# 25 points, equal probability for each type, uniformly distributed locations

rmpoint(25, function(x,y,m) {rep(1, length(x))}, types=abc)
# same as above
rmpoint(25, list(function(x,y){rep(1, length(x))},
                 function(x,y){rep(1, length(x))},
                 function(x,y){rep(1, length(x))}),
             types=abc)
# same as above

rmpoint(25, function(x,y,m) { x }, types=abc)
# 25 points, equal probability for each type,
# locations nonuniform with density proportional to x

rmpoint(25, function(x,y,m) { ifelse(m == "a", 1, x) }, types=abc)
rmpoint(25, list(function(x,y) { rep(1, length(x)) },
                function(x,y) { x },
                function(x,y) { x }),
                types=abc)
# 25 points, UNEQUAL probabilities for each type,
# type "a" points uniformly distributed,
# type "b" and "c" points nonuniformly distributed.

##### Model II

rmpoint(25, 1, types=abc, ptypes=rep(1,3)/3)
rmpoint(25, 1, types=abc, ptypes=rep(1,3))
# 25 points, equal probability for each type,
# uniformly distributed locations

rmpoint(25, function(x,y,m) {rep(1, length(x))}, types=abc, ptypes=rep(1,3))
# same as above
rmpoint(25, list(function(x,y){rep(1, length(x))},
                 function(x,y){rep(1, length(x))},
                 function(x,y){rep(1, length(x))}),
             types=abc, ptypes=rep(1,3))
# same as above

rmpoint(25, function(x,y,m) { x }, types=abc, ptypes=rep(1,3))
# 25 points, equal probability for each type,
# locations nonuniform with density proportional to x

rmpoint(25, function(x,y,m) { ifelse(m == "a", 1, x) }, types=abc, ptypes=rep(1,3))
# 25 points, EQUAL probabilities for each type,
# type "a" points uniformly distributed,
# type "b" and "c" points nonuniformly distributed.

###### Model III

rmpoint(c(12, 8, 4), 1, types=abc)
# 12 points of type "a",
# 8 points of type "b",
# 4 points of type "c",
# each uniformly distributed

rmpoint(c(12, 8, 4), function(x,y,m) { ifelse(m=="a", 1, x)}, types=abc)
rmpoint(c(12, 8, 4), list(function(x,y) { rep(1, length(x)) },
                               function(x,y) { x },
                               function(x,y) { x }),
             types=abc)

# 12 points of type "a", uniformly distributed
# 8 points of type "b", nonuniform
# 4 points of type "c", nonuniform


#########

## Randomising an existing point pattern:
# same numbers of points of each type, uniform random locations (Model III)
rmpoint(table(marks(demopat)), 1, win=Window(demopat))

# same total number of points, distribution of types estimated from X,
# uniform random locations (Model II)
rmpoint(npoints(demopat), 1, types=levels(marks(demopat)), win=Window(demopat),
       ptypes=table(marks(demopat)))

}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}

