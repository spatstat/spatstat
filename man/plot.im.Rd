\name{plot.im}
\alias{plot.im}
\alias{image.im}
\title{Plot a Pixel Image}
\description{
  Plot a pixel image.
}
\usage{
  \method{plot}{im}(x, \dots,
                   main,
                   add=FALSE, clipwin=NULL,
                   col=NULL, valuesAreColours=NULL, log=FALSE,
                   ribbon=show.all, show.all=!add,
                   ribside=c("right", "left", "bottom", "top"),
                   ribsep=0.15, ribwid=0.05, ribn=1024,
                   ribscale=1, ribargs=list(), colargs=list(),
                   useRaster=NULL,
                   do.plot=TRUE) 

  \method{image}{im}(x, \dots,
                   main,
                   add=FALSE, clipwin=NULL,
                   col=NULL, valuesAreColours=NULL, log=FALSE,
                   ribbon=show.all, show.all=!add,
                   ribside=c("right", "left", "bottom", "top"),
                   ribsep=0.15, ribwid=0.05, ribn=1024,
                   ribscale=1, ribargs=list(), colargs=list(),
                   useRaster=NULL,
                   do.plot=TRUE) 
}
\arguments{
  \item{x}{
    The pixel image to be plotted.
    An object of class \code{"im"} (see \code{\link{im.object}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link[graphics]{image.default}}
    to control the plot. See Details.
  }
  \item{main}{Main title for the plot.}
  \item{add}{
    Logical value indicating whether to superimpose the image on the
    existing plot (\code{add=TRUE})
    or to initialise a new plot (\code{add=FALSE}, the default).
  }
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the image will be displayed.
  }
  \item{col}{
    Colours for displaying the pixel values.
    Either a character vector of colour values,
    an object of class \code{\link{colourmap}},
    or a \code{function} as described under Details.
  }
  \item{valuesAreColours}{
    Logical value. If \code{TRUE}, the pixel values of \code{x}
    are to be interpreted as colour values.
  }
  \item{log}{
    Logical value. If \code{TRUE}, the colour map will be
    evenly-spaced on a logarithmic scale.
  }
  \item{ribbon}{
    Logical flag indicating whether to display a ribbon
    showing the colour map. Default is \code{TRUE}
    for new plots and \code{FALSE} for added plots.
  }
  \item{show.all}{
    Logical value indicating whether to display all plot elements
    including the main title and colour ribbon. Default is \code{TRUE}
    for new plots and \code{FALSE} for added plots.
  }
  \item{ribside}{
    Character string indicating where to display the ribbon
    relative to the main image.
  }
  \item{ribsep}{
    Factor controlling the space between the ribbon and the image.
  }
  \item{ribwid}{
    Factor controlling the width of the ribbon.
  }
  \item{ribn}{
    Number of different values to display in the ribbon.
  }
  \item{ribscale}{
    Rescaling factor for tick marks. The values on the numerical scale
    printed beside the ribbon will be multiplied by this rescaling factor.
  }
  \item{ribargs}{
    List of additional arguments passed to
    \code{\link[graphics]{image.default}} and
    \code{\link[graphics]{axis}} 
    to control the display of the ribbon and its scale axis. These may override
    the \code{\dots} arguments. 
  }
  \item{colargs}{
    List of additional arguments passed to
    \code{col} if it is a function.
  }
  \item{useRaster}{
    Logical value, passed to \code{\link[graphics]{image.default}}.
    Images are plotted using a bitmap raster if \code{useRaster=TRUE}
    or by drawing polygons if \code{useRaster=FALSE}.
    Bitmap raster display tends to produce better results,
    but is not supported on all graphics devices.
    The default is to use bitmap raster display if it is supported.
  }
  \item{do.plot}{
    Logical value indicating whether to actually plot the image
    and colour ribbon.
    Setting \code{do.plot=FALSE} will simply return the
    colour map and the bounding box that were chosen for the plot.
  }
}
\value{
  The colour map used. An object of class \code{"colourmap"}.

  Also has an attribute \code{"bbox"} giving a bounding box
  for the colour image (including the ribbon if present).
}
\details{
  This is the \code{plot} method for the class \code{"im"}.
  [It is also the \code{image} method for \code{"im"}.]

  The pixel image \code{x} is displayed on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  If \code{ribbon=TRUE}, a legend will be plotted.
  The legend consists of a colour ribbon and an axis with tick-marks,
  showing the correspondence between the pixel values and the colour map.

  By default, the ribbon is placed at the right of the main image.
  This can be changed using the argument \code{ribside}.
  
  Arguments \code{ribsep, ribwid, ribn} control the appearance of the
  ribbon.
  The width of the ribbon is \code{ribwid} times the size of the pixel
  image, where `size' means the larger of the width and the height.
  The distance separating the ribbon and the image is \code{ribsep} times
  the size of the pixel image. The ribbon contains \code{ribn}
  different numerical values, evenly spaced between the minimum and
  maximum pixel values in the image \code{x}, rendered according to
  the chosen colour map.

  Arguments \code{ribscale, ribargs} control the annotation of the
  colour ribbon. To plot the colour ribbon without the axis and
  tick-marks, use \code{ribargs=list(axes=FALSE)}.

  Normally the pixel values are displayed using the colours given in the
  argument \code{col}. This may be either
  \itemize{
    \item an explicit colour map (an object of class
    \code{"colourmap"}, created by the command \code{\link{colourmap}}).
    This is the best way to ensure
    that when we plot different images, the colour maps are consistent.
    \item a character vector or integer vector
    that specifies a set of colours.
    The colour mapping will be stretched to match the range of
    pixel values in the image \code{x}. The mapping of pixel values
    to colours is determined as follows.
    \describe{
      \item{logical-valued images:}{the values \code{FALSE} and
	\code{TRUE} are mapped to the colours \code{col[1]} and
	\code{col[2]} respectively.
	The vector \code{col} should have length 2.
      }
      \item{factor-valued images:}{the factor levels \code{levels(x)}
	are mapped to the entries of \code{col} in order. The vector
	\code{col} should have the same length as \code{levels(x)}.
      }
      \item{numeric-valued images:}{
	By default, the range of pixel values in \code{x}
	is divided into \code{n = length(col)}
	equal subintervals, which are mapped to the colours in \code{col}.
	(If \code{col} was not specified, it defaults to a vector of 255
	colours.)
      
	Alternatively if the argument \code{zlim} is given, it should be
	a vector of length 2 specifying an interval of real numbers.
	This interval will be used instead of the range of pixel
	values. The interval from \code{zlim[1]} to \code{zlim[2]} will be
	mapped to the colours in \code{col}. This facility enables the user to
	plot several images using a consistent colour map.

	Alternatively if the argument \code{breaks} is given,
	then this specifies
	the endpoints of the subintervals that are mapped to each colour.
	This is incompatible with \code{zlim}.
      
	The arguments \code{col} and \code{zlim} or \code{breaks}
	are then passed to the function \code{\link{image.default}}.
	For examples of the use of these arguments,
	see \code{\link{image.default}}.
      }
    }
    \item {
      a \code{function} in the \R language
      with an argument named \code{range}
      or \code{inputs}.
    
      If \code{col} is a function with an argument named \code{range},
      and if the pixel values of \code{x} are numeric values, 
      then the colour values will be determined by evaluating
      \code{col(range=range(x))}. The result of this evaluation
      should be a character vector containing colour values, or
      a \code{"colourmap"} object. Examples of such functions
      are \code{\link{beachcolours}} and \code{\link{beachcolourmap}}.
  
      If \code{col} is a function with an argument named \code{inputs},
      and if the pixel values of \code{x} are discrete values (integer,
      logical, factor or character), 
      then the colour values will be determined by evaluating
      \code{col(inputs=p)} where \code{p} is the set of possible pixel
      values. The result should be a character vector 
      containing colour values, or a \code{"colourmap"} object.
    }
    \item{
      a \code{function} in the \R language with first argument
      named \code{n}.
      The colour values will be determined by evaluating
      \code{col(n)} where \code{n} is the
      number of distinct pixel values, up to a maximum of 128.
      The result of this evaluation
      should be a character vector containing color values.
      Examples of such functions are
      \code{\link[grDevices]{heat.colors}},
      \code{\link[grDevices]{terrain.colors}},
      \code{\link[grDevices]{topo.colors}} and
      \code{\link[grDevices]{cm.colors}}. 
   }
  }

  If \code{spatstat.options("monochrome")} has been set to \code{TRUE}
  then \bold{all colours will be converted to grey scale values}.
  
  Other graphical parameters controlling the display of both the pixel image
  and the ribbon can be passed through the \code{...} arguments
  to the function \code{\link[graphics]{image.default}}.
  A parameter is handled only if it is one of the following:
  \itemize{
    \item a formal argument of \code{\link[graphics]{image.default}}
    that is operative when \code{add=TRUE}.
    \item one of the
    parameters \code{"main", "asp", "sub", "axes", "ann",
      "cex", "font", "cex.axis", "cex.lab", "cex.main", "cex.sub",
      "col.axis", "col.lab", "col.main", "col.sub",
      "font.axis", "font.lab", "font.main", "font.sub"}
    described in \code{\link[graphics]{par}}.
    \item the argument \code{box}, a logical value specifying whether
    a box should be drawn.
  }
  Images are plotted using a bitmap raster if \code{useRaster=TRUE}
  or by drawing polygons if \code{useRaster=FALSE}.
  Bitmap raster display (performed by \code{\link[graphics]{rasterImage}})
  tends to produce better results, but is not supported on all graphics devices.
  The default is to use bitmap raster display if it is
  supported according to \code{\link[grDevices]{dev.capabilities}}.

  Alternatively, the pixel values could be directly 
  interpretable as colour values in \R. That is, the pixel values
  could be character strings that represent colours, or
  values of a factor whose levels are character strings representing
  colours.
  \itemize{
    \item If \code{valuesAreColours=TRUE}, then the pixel values will
    be interpreted as colour values and displayed using these colours.
    \item If \code{valuesAreColours=FALSE}, then the pixel values will
    \emph{not} be interpreted as colour values, even if they could be.
    \item If \code{valuesAreColours=NULL}, the algorithm will guess
    what it should do. If the argument \code{col} is given,
    the pixel values will \emph{not} be interpreted as colour values. Otherwise,
    if all the pixel values are strings that represent colours, then
    they will be interpreted and displayed as colours.
  }
  If pixel values are interpreted as colours, 
  the arguments \code{col} and \code{ribbon} will be ignored,
  and a ribbon will not be plotted.

}
\section{Complex-valued images}{
  If the pixel values in \code{x} are complex numbers,
  they will be converted into four images containing the real and
  imaginary parts and the modulus and argument, 
  and plotted side-by-side using \code{\link{plot.imlist}}.
}
\section{Monochrome colours}{
  If \code{spatstat.options("monochrome")} has been set to \code{TRUE},
  then \bold{the image will be plotted in greyscale}.
  The colours are converted to grey scale values using
  \code{\link{to.grey}}.
  The choice of colour map still has an effect, since it determines
  the final grey scale values.

  Monochrome display can also be achieved by
  setting the graphics device parameter \code{colormodel="grey"}
  when starting a new graphics device, or in a call to
  \code{\link{ps.options}} or \code{\link{pdf.options}}.
}
\section{Image Rendering Errors and Problems}{
  The help for \code{\link[graphics]{image.default}}
  and \code{\link[graphics]{rasterImage}} explains that
  errors may occur, or images may be rendered incorrectly, on some
  devices, depending on the availability of colours and other
  device-specific constraints.

  If the image is not displayed at all,
  try setting \code{useRaster=FALSE} in the call to \code{plot.im}.
  If the ribbon colours are not displayed, set
  \code{ribargs=list(useRaster=FALSE)}.
  
  Errors may occur on some graphics devices if the image is very
  large. If this happens, try setting \code{useRaster=FALSE} in the
  call to \code{plot.im}.

  The error message
  \code{useRaster=TRUE can only be used with a regular grid}
  means that the \eqn{x} and \eqn{y} coordinates of the pixels in the
  image are not perfectly equally spaced, due to numerical rounding.
  This occurs with some images created by earlier versions of \pkg{spatstat}.
  To repair the coordinates in an image \code{X}, type
  \code{X <- as.im(X)}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{colourmap}},
  \code{\link{contour.im}},
  \code{\link{persp.im}},
  \code{\link[graphics]{image.default}},
  \code{\link{spatstat.options}}
}
\examples{
   # an image
   Z <- setcov(owin())
   plot(Z)
   plot(Z, ribside="bottom")
   # stretchable colour map
   plot(Z, col=terrain.colors(128), axes=FALSE)
   # fixed colour map
   tc <- colourmap(rainbow(128), breaks=seq(-1,2,length=129))
   plot(Z, col=tc)
   # colour map function, with argument 'range'
   plot(Z, col=beachcolours, colargs=list(sealevel=0.5))
   # tweaking the plot
   plot(Z, main="La vie en bleu", col.main="blue", cex.main=1.5,
        box=FALSE,
        ribargs=list(col.axis="blue", col.ticks="blue", cex.axis=0.75))
   # log scale
   V <- eval.im(exp(exp(Z+2))/1e4)
   plot(V, log=TRUE, main="Log scale")
   # it's complex
   Y <- exp(Z + V * 1i)
   plot(Y)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
 
\keyword{hplot}
