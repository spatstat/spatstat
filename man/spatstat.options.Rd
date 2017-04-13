\name{spatstat.options}
\alias{spatstat.options}
\alias{reset.spatstat.options}
\title{Internal Options in Spatstat Package}
\description{
  Allows the user to examine and reset the values
  of global parameters which control actions in the
  \pkg{spatstat} package.
}
\usage{
  spatstat.options(...)
  reset.spatstat.options()
}
\arguments{
  \item{\dots}{
    Either empty,
    or a succession of parameter names in quotes,
    or a succession of \code{name=value} pairs.
    See below for the parameter names.
  }
}
\value{
  Either a list of parameters and their values,
  or a single value. See Details.
}
\details{
  The function \code{spatstat.options}
  allows the user to examine and reset the values
  of global parameters which control actions in the
  \pkg{spatstat} package.
  It is analogous to the system function \code{\link[base]{options}}.

  The function \code{reset.spatstat.options} resets all the
  global parameters in \pkg{spatstat} to their original,
  default values.

  The global parameters of interest to the user are:
  \describe{
    \item{checkpolygons}{
      Logical flag indicating whether the functions
      \code{\link[spatstat]{owin}} and \code{\link[spatstat]{as.owin}} should
      apply very strict checks on the validity of polygon data. 
      These strict checks are no longer necessary, and the default is
      \code{checkpolygons=FALSE}.
      See also \code{fixpolygons} below.
    }
    \item{checksegments}{
      Logical flag indicating whether the functions
      \code{\link[spatstat]{psp}} and \code{\link[spatstat]{as.psp}} should check
      the validity of line segment data (in particular, checking that
      the endpoints of the line segments are inside the specified
      window). It is advisable to leave this flag set to \code{TRUE}. 
    }
    \item{eroded.intensity}{
      Logical flag affecting the behaviour of the
      score and pseudo-score residual functions
      \code{\link[spatstat]{Gcom}}, \code{\link[spatstat]{Gres}}
      \code{\link[spatstat]{Kcom}}, \code{\link[spatstat]{Kres}}, 
      \code{\link[spatstat]{psstA}}, 
      \code{\link[spatstat]{psstG}}, 
      \code{\link[spatstat]{psst}}.
      The flag indicates whether to compute intensity estimates
      on an eroded window (\code{eroded.intensity=TRUE})
      or on the original data window
      (\code{eroded.intensity=FALSE}, the default). 
    }
    \item{expand}{
      The default expansion factor (area inflation factor)
      for expansion of the simulation window in
      \code{\link[spatstat]{rmh}} (see \code{\link[spatstat]{rmhcontrol}}).
      Initialised to \code{2}.
    }
    \item{expand.polynom}{
      Logical. Whether expressions involving \code{\link[spatstat]{polynom}}
      in a model formula should be expanded, so that
      \code{polynom(x,2)} is replaced by \code{x + I(x^2)} and so on.
      Initialised to \code{TRUE}.
    }
    \item{fastpois}{
      Logical. Whether to use a fast algorithm
      (introduced in \pkg{spatstat 1.42-3}) for
      simulating the Poisson point process in \code{\link[spatstat]{rpoispp}}
      when the argument \code{lambda} is a pixel image.
      Initialised to \code{TRUE}.
      Should be set to \code{FALSE} if needed to guarantee repeatability
      of results computed using earlier versions of \pkg{spatstat}.
    }
    \item{fastthin}{
      Logical. Whether to use a fast C language algorithm
      (introduced in \pkg{spatstat 1.42-3}) for
      random thinning in \code{\link[spatstat]{rthin}} when the argument
      \code{P} is a single number.
      Initialised to \code{TRUE}.
      Should be set to \code{FALSE} if needed to guarantee repeatability
      of results computed using earlier versions of \pkg{spatstat}.
    }
    \item{fastK.lgcp}{
      Logical. Whether to use fast or slow algorithm to compute the
      (theoretical) \eqn{K}-function of a log-Gaussian Cox process
      for use in \code{\link[spatstat]{lgcp.estK}} or \code{\link[spatstat]{Kmodel}}.
      The slow algorithm uses accurate numerical integration; the
      fast algorithm uses Simpson's Rule for numerical integration,
      and is about two orders of magnitude faster.
      Initialised to \code{FALSE}. 
    }
    \item{fixpolygons}{
      Logical flag indicating whether the functions
      \code{\link[spatstat]{owin}} and \code{\link[spatstat]{as.owin}} should
      repair errors in polygon data. For example,
      self-intersecting polygons and overlapping polygons will
      be repaired. The default is \code{fixpolygons=TRUE}.
    }
    \item{fftw}{
      Logical value indicating whether the two-dimensional
      Fast Fourier Transform
      should be computed using the package \pkg{fftwtools},
      instead of the \code{fft} function in the \pkg{stats} package.
      This affects the speed of \code{\link[spatstat]{density.ppp}},
      \code{\link[spatstat]{density.psp}}, \code{\link[spatstat]{blur}}
      \code{\link[spatstat]{setcov}} and \code{\link[spatstat]{Smooth.ppp}}.
    }
    \item{gpclib}{
      Defunct. This parameter was used to permit or forbid the use of the
      package \pkg{gpclib}, because of its restricted software licence.
      This package is no longer needed.
    }
    \item{huge.npoints}{
      The maximum value of \code{n} for which \code{runif(n)} will
      not generate an error (possible errors include failure to allocate
      sufficient memory, and integer overflow of \code{n}).
      An attempt to generate more than this number of random points
      triggers a warning from \code{\link[spatstat]{runifpoint}} and other
      functions. Defaults to \code{1e6}.
    }
    \item{image.colfun}{
      Function determining the default colour map for
      \code{\link[spatstat]{plot.im}}. When called with one integer argument
      \code{n}, this function should return a character vector of length
      \code{n} specifying \code{n} different colours.
    }
    \item{Kcom.remove.zeroes}{
      Logical value, determining whether the algorithm in
      \code{\link[spatstat]{Kcom}} and \code{\link[spatstat]{Kres}}
      removes or retains the contributions to the function
      from pairs of points that are identical. If these are
      retained then the function has a jump at \eqn{r=0}.
      Initialised to \code{TRUE}.
    }
    \item{maxedgewt}{
      Edge correction weights will be trimmed
      so as not to exceed this value.
      This applies to the weights computed by
      \code{\link[spatstat]{edge.Trans}} or \code{\link[spatstat]{edge.Ripley}}
      and used in \code{\link[spatstat]{Kest}} and its relatives.
    }
    \item{maxmatrix}{
      The maximum permitted size (rows times columns)
      of matrices generated by \pkg{spatstat}'s internal code.
      Used by \code{\link[spatstat]{ppm}} and \code{\link[spatstat]{predict.ppm}}
      (for example) to decide when to
      split a large calculation into blocks.
      Defaults to \code{2^24=16777216}. 
    }
    \item{monochrome}{
      Logical flag indicating whether graphics should be
      plotted in grey scale (\code{monochrome=TRUE}) or
      in colour (\code{monochrome=FALSE}, the default).
    }
    \item{n.bandwidth}{
      Integer. Number of trial values of smoothing bandwidth to use for
      cross-validation in \code{\link[spatstat]{bw.relrisk}} and similar
      functions.
    }
    \item{ndummy.min}{
      The minimum number of dummy points in a quadrature scheme
      created by \code{\link[spatstat]{default.dummy}}.
      Either an integer or a pair of integers
      giving the minimum number of dummy points in the \code{x} and \code{y}
      directions respectively. 
    }
    \item{ngrid.disc}{
      Number of points in the square grid used to compute
      a discrete approximation to the areas of discs
      in \code{\link[spatstat]{areaLoss}} and \code{\link[spatstat]{areaGain}}
      when exact calculation is not available.
      A single integer.
    }
    \item{npixel}{
      Default number of pixels in a binary mask or pixel image.
      Either an integer, or a pair of integers,
      giving the number of pixels in the \code{x} and \code{y}
      directions respectively.
    }
    \item{nvoxel}{
      Default number of voxels in a 3D image,
      typically for calculating the distance transform in \code{\link[spatstat]{F3est}}.
      Initialised to 4 megavoxels: \code{nvoxel = 2^22 = 4194304}.
    }
    \item{par.binary}{
      List of arguments to be passed to the function
      \code{\link[graphics]{image}}
      when displaying a binary image mask (in \code{\link[spatstat]{plot.owin}}
      or \code{\link[spatstat]{plot.ppp}}).
      Typically used to reset the colours of foreground and background.
    }
    \item{par.contour}{
      List of arguments controlling contour plots of pixel images
      by \code{\link[spatstat]{contour.im}}.
    }
    \item{par.fv}{
      List of arguments controlling the plotting of functions 
      by \code{\link[spatstat]{plot.fv}} and its relatives.
    }
    \item{par.persp}{
      List of arguments to be passed to the function
      \code{\link[graphics]{persp}}
      when displaying a real-valued image, such as the fitted surfaces
      in \code{\link[spatstat]{plot.ppm}}.
    }
    \item{par.points}{
      List of arguments controlling the plotting of point patterns
      by \code{\link[spatstat]{plot.ppp}}.
    }
    \item{par.pp3}{
      List of arguments controlling the plotting of three-dimensional
      point patterns by \code{\link[spatstat]{plot.pp3}}.
    }
    \item{print.ppm.SE}{
      Default rule used by \code{\link[spatstat]{print.ppm}}
      to decide whether to calculate and print standard errors
      of the estimated coefficients of the model.
      One of the strings \code{"always"}, \code{"never"} or
      \code{"poisson"} (the latter indicating that standard errors
      will be calculated only for Poisson models). The default is
      \code{"poisson"} because the calculation for non-Poisson models
      can take a long time.
    }
    \item{progress}{
      Character string determining the style of progress reports
      printed by \code{\link[spatstat]{progressreport}}. Either
      \code{"tty"}, \code{"tk"} or \code{"txtbar"}.
      For explanation of these options, see \code{\link[spatstat]{progressreport}}.
    }
    \item{project.fast}{
      Logical. If \code{TRUE}, the algorithm of
      \code{\link[spatstat]{project.ppm}}
      will be accelerated using a shorcut.
      Initialised to \code{FALSE}.
    }
    \item{psstA.ngrid}{
      Single integer, 
      controlling the accuracy of the discrete approximation
      of areas computed in the function \code{\link[spatstat]{psstA}}.
      The area of a disc is approximated by counting points on
      an \eqn{n \times n}{n * n} grid.
      Initialised to 32.
    }
    \item{psstA.nr}{
      Single integer,
      determining the number of distances \eqn{r}
      at which the function \code{\link[spatstat]{psstA}} will be evaluated
      (in the default case where argument \code{r} is absent).
      Initialised to 30.
    }
    \item{psstG.remove.zeroes}{
      Logical value, determining whether the algorithm in
      \code{\link[spatstat]{psstG}} 
      removes or retains the contributions to the function
      from pairs of points that are identical. If these are
      retained then the function has a jump at \eqn{r=0}.
      Initialised to \code{TRUE}.
    }
    \item{rmh.p, rmh.q, rmh.nrep}{
      New default values for the parameters \code{p}, \code{q}
      and \code{nrep} in the Metropolis-Hastings simulation
      algorithm. These override the defaults
      in \code{\link[spatstat]{rmhcontrol.default}}.
    }
    \item{scalable}{
      Logical flag indicating whether the new code in \code{rmh.default}
      which makes the results scalable (invariant to change of units)
      should be used.  In order to recover former behaviour (so that
      previous results can be reproduced) set this option equal to
      \code{FALSE}.  See the \dQuote{Warning} section in the help for
      \code{\link[spatstat]{rmh}()} for more detail.
    }
    \item{terse}{
      Integer between 0 and 4.
      The level of terseness (brevity) in printed output
      from many functions in \pkg{spatstat}.
      Higher values mean shorter output.
      A rough guide is the following:
      \tabular{ll}{
	0 \tab Full output\cr
        1 \tab Avoid wasteful output \cr
        2 \tab Remove space between paragraphs\cr
        3 \tab Suppress extras such as standard errors \cr
        4 \tab Compress text, suppress internal warnings
      }
      The value of \code{terse} is initialised to 0.
    }
    \item{transparent}{
      Logical value indicating whether default colour maps
      are allowed to include semi-transparent colours, where possible.
      Default is \code{TRUE}.
      Currently this only affects \code{\link[spatstat]{plot.ppp}}. 
    }
    \item{units.paren}{
      The kind of parenthesis which encloses the text that
      explains a \code{unitname}. This text is seen in
      the text output of functions like \code{\link[spatstat]{print.ppp}} and
      in the graphics generated by \code{\link[spatstat]{plot.fv}}.
      The value should be one of the character strings \code{'('}, \code{'['},
      \code{'{'} or \code{''}. The default is \code{'('}.
    }
  }
  If no arguments are given, 
  the current values of all parameters are returned,
  in a list.

  If one parameter name is given, the current value of this
  parameter is returned (\bold{not} in a list, just the value).
  
  If several parameter names are given,
  the current values of these parameters are returned, in a list.
  
  If \code{name=value} pairs are given, the named parameters
  are reset to the given values, and the \bold{previous} values of
  these parameters are returned, in a list.
}
\section{Internal parameters}{
  The following parameters may also be specified to
  \code{spatstat.options} but are intended for
  software development or testing purposes.
  \describe{
    \item{closepairs.newcode}{
      Logical. Whether to use new version of the
      code for \code{\link[spatstat]{closepairs}}.
      Initialised to \code{TRUE}.
    }
    \item{crossing.psp.useCall}{
      Logical. Whether to use new version of the
      code for \code{\link[spatstat]{crossing.psp}}.
      Initialised to \code{TRUE}.
    }
    \item{crosspairs.newcode}{
      Logical. Whether to use new version of the
      code for \code{\link[spatstat]{crosspairs}}.
      Initialised to \code{TRUE}.
    }
    \item{densityC}{
      Logical. 
      Indicates whether to use accelerated C code
      (\code{densityC=TRUE}) or interpreted R code (\code{densityC=FALSE})
      to evaluate \code{density.ppp(X, at="points")}.
      Initialised to \code{TRUE}.
    }
    \item{exactdt.checks.data}{
      Logical.
      Do not change this value, unless you are
      \adrian. 
    }
    \item{fasteval}{
      One of the strings \code{'off'}, \code{'on'} or
      \code{'test'} determining whether to use accelerated C code
      to evaluate the conditional intensity of a Gibbs model.
      Initialised to \code{'on'}.
    }
    \item{old.morpho.psp}{
      Logical. Whether to use old R code for
      morphological operations.
      Initialise to \code{FALSE}.
    }
    \item{selfcrossing.psp.useCall}{
      Logical. Whether to use new version of the
      code for \code{\link[spatstat]{selfcrossing.psp}}.
      Initialised to \code{TRUE}.
    }
    \item{use.Krect}{
      Logical. Whether to use new code for the
      K-function in a rectangular window.
      Initialised to \code{TRUE}.
    }
  }
}
\seealso{
  \code{\link[base]{options}}
}
\examples{
  # save current values
  oldopt <- spatstat.options()

  spatstat.options("npixel")
  spatstat.options(npixel=150)
  spatstat.options(npixel=c(100,200))

  spatstat.options(par.binary=list(col=grey(c(0.5,1))))

  spatstat.options(par.persp=list(theta=-30,phi=40,d=4))
  # see help(persp.default) for other options

  # revert
  spatstat.options(oldopt)

}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
