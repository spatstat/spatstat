\name{spatstat-package}
\alias{spatstat-package} 
\alias{spatstat} 
\docType{package}
\title{The Spatstat Package}
\description{
  This is a summary of the features of 
  \pkg{spatstat}, a package in \R 
  for the statistical analysis of spatial point patterns.
}
\details{
  \pkg{spatstat} is a package for the statistical analysis
  of spatial data. Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.
  The points may carry auxiliary data (`marks'),
  and the spatial region in which the points were recorded 
  may have arbitrary shape. 

  The package is designed to support a complete statistical analysis
  of spatial data. It supports
  \itemize{
    \item creation, manipulation and plotting of point patterns;
    \item exploratory data analysis;
    \item spatial random sampling;
    \item simulation of point process models;
    \item parametric model-fitting;
    \item non-parametric smoothing and regression;
    \item formal inference (hypothesis tests, confidence intervals);
    \item model diagnostics.
  }
  
  Apart from two-dimensional point patterns and point processes,
  \pkg{spatstat} also supports point patterns in three dimensions, 
  point patterns in multidimensional space-time,
  point patterns on a linear network,
  patterns of line segments in two dimensions, 
  and spatial tessellations and random sets in two dimensions.

  The package can fit several types of point process models
  to a point pattern dataset:
  \itemize{
    \item Poisson point process models
    (by Berman-Turner approximate maximum likelihood
    or by spatial logistic regression)
    \item Gibbs/Markov point process models
    (by Baddeley-Turner approximate maximum pseudolikelihood,
    Coeurjolly-Rubak logistic likelihood, or
    Huang-Ogata approximate maximum likelihood)
    \item Cox/cluster point process models
    (by Waagepetersen's two-step fitting procedure
    and minimum contrast, composite likelihood, or Palm likelihood)
    \item determinantal point process models
    (by Waagepetersen's two-step fitting procedure
    and minimum contrast, composite likelihood, or Palm likelihood)
  }
  The models may include spatial trend,
  dependence on covariates, and complicated interpoint interactions.
  Models are specified by 
  a \code{formula} in the \R  language, and are fitted using
  a function analogous to \code{\link{lm}} and \code{\link{glm}}.
  Fitted models can be printed, plotted, predicted, simulated and so on.
}
\section{Getting Started}{
  For a quick introduction to \pkg{spatstat}, read
  the package vignette \emph{Getting started with spatstat}
  installed with \pkg{spatstat}. To read that document, you can either
  \itemize{
    \item visit \url{cran.r-project.org/web/packages/spatstat}
    and click on \code{Getting Started with Spatstat}
    \item
    start \R, type \code{library(spatstat)} and \code{vignette('getstart')}
    \item
    start \R, type \code{help.start()} to open the help
    browser, and navigate to \code{Packages > spatstat > Vignettes}.
  }
  Once you have installed \pkg{spatstat}, start \R and type
  \code{library(spatstat)}. Then type \code{beginner}
  for a beginner's introduction, or 
  \code{demo(spatstat)} for a demonstration of the package's capabilities.

  For a complete course on \pkg{spatstat},
  and on statistical analysis of spatial point patterns,
  read the book by Baddeley, Rubak and Turner (2015).
  Other recommended books on spatial point process methods
  are Diggle (2014), Gelfand et al (2010) and Illian et al (2008).

  The \pkg{spatstat} package includes over 50 datasets,
  which can be useful when learning the package.
  Type \code{demo(data)} to see plots of all datasets
  available in the package.
  Type \code{vignette('datasets')} for detailed background information
  on these datasets, and plots of each dataset.
      
  For information on converting your data into \pkg{spatstat} format,
  read Chapter 3 of Baddeley, Rubak and Turner (2015).
  This chapter is available free online, as one of the sample chapters
  at the book companion website, \url{spatstat.github.io/book}.

  For information about handling data in \bold{shapefiles},
  see Chapter 3, or the Vignette
  \emph{Handling shapefiles in the spatstat package},
  installed with \pkg{spatstat}, accessible as
  \code{vignette('shapefiles')}.
}
\section{Updates}{
  New versions of \pkg{spatstat} are released every 8 weeks.
  Users are advised to update their installation of \pkg{spatstat}
  regularly.
  
  Type \code{latest.news} to read the news documentation about
  changes to the current installed version of \pkg{spatstat}.

  See the Vignette \emph{Summary of recent updates},
  installed with \pkg{spatstat}, which describes the main changes
  to \pkg{spatstat} since the book (Baddeley, Rubak and Turner, 2015)
  was published. It is accessible as \code{vignette('updates')}.
  
  Type \code{news(package="spatstat")} to read news documentation about
  all previous versions of the package.
}
\section{FUNCTIONS AND DATASETS}{
  Following is a summary of the main functions and datasets
  in the \pkg{spatstat} package.
  Alternatively an alphabetical list of all functions and
  datasets is available by typing \code{library(help=spatstat)}.

  For further information on any of these,
  type \code{help(name)} or \code{?name}
  where \code{name} is the name of the function
  or dataset.
}
\section{CONTENTS:}{
  \tabular{ll}{
    I. \tab Creating and manipulating data \cr
    II. \tab Exploratory Data Analysis \cr
    III. \tab Model fitting (Cox and cluster models) \cr
    IV. \tab Model fitting (Poisson and Gibbs models) \cr
    V. \tab Model fitting (determinantal point processes)\cr
    VI. \tab Model fitting (spatial logistic regression)\cr
    VII. \tab Simulation \cr
    VIII. \tab Tests and diagnostics\cr
    IX. \tab Documentation
  }
}

\section{I. CREATING AND MANIPULATING DATA}{

  \bold{Types of spatial data:}

  The main types of spatial data supported by \pkg{spatstat} are:

  \tabular{ll}{
    \code{\link{ppp}} \tab point pattern \cr
    \code{\link{owin}} \tab window (spatial region) \cr
    \code{\link{im}} \tab pixel image \cr
    \code{\link{psp}} \tab line segment pattern \cr
    \code{\link{tess}} \tab tessellation \cr
    \code{\link{pp3}} \tab three-dimensional point pattern \cr
    \code{\link{ppx}} \tab point pattern in any number of dimensions \cr
    \code{\link{lpp}} \tab point pattern on a linear network
  }

  \bold{To create a point pattern:}
  
  \tabular{ll}{
    \code{\link{ppp}} \tab
    create a point pattern from \eqn{(x,y)} and window information
    \cr
    \tab
    \code{ppp(x, y, xlim, ylim)} for rectangular window\cr
    \tab
    \code{ppp(x, y, poly)} for polygonal window \cr
    \tab
    \code{ppp(x, y, mask)} for binary image window \cr
    \code{\link{as.ppp}} \tab
    convert other types of data to a \code{ppp} object \cr
    \code{\link{clickppp}} \tab
    interactively add points to a plot \cr
    \code{\link{marks<-}}, \code{\%mark\%}  \tab
    attach/reassign marks to a point pattern
  }
      
  \bold{To simulate a random point pattern:}
  
  \tabular{ll}{
	\code{\link{runifpoint}} \tab
	generate \eqn{n} independent uniform random points \cr
	\code{\link{rpoint}} \tab
	generate \eqn{n} independent random points \cr
	\code{\link{rmpoint}} \tab
	generate \eqn{n} independent multitype random points \cr
	\code{\link{rpoispp}} \tab
	simulate the (in)homogeneous Poisson point process \cr
	\code{\link{rmpoispp}} \tab
	simulate the (in)homogeneous multitype Poisson point process \cr
	\code{\link{runifdisc}} \tab
	generate \eqn{n} independent uniform random points in disc\cr
	\code{\link{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link{rsyst}} \tab
	systematic random sample of points \cr
	\code{\link{rjitter}} \tab
	apply random displacements to points in a pattern\cr
	\code{\link{rMaternI}}  \tab
	simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Model I inhibition process\cr
	\code{\link{rMaternII}} \tab
	simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Model II inhibition process\cr
	\code{\link{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link{rStrauss}} \tab
	simulate Strauss process (perfect simulation)\cr
	\code{\link{rHardcore}} \tab
	simulate Hard Core process (perfect simulation)\cr
	\code{\link{rStraussHard}} \tab
	simulate Strauss-hard core process (perfect simulation)\cr
	\code{\link{rDiggleGratton}} \tab
	simulate Diggle-Gratton process (perfect simulation)\cr
	\code{\link{rDGS}} \tab
	simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
	\code{\link{rPenttinen}} \tab
	simulate Penttinen process (perfect simulation)\cr
	\code{\link{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link{rPoissonCluster}} \tab
	simulate a general Poisson cluster process\cr
	\code{\link{rMatClust}} \tab
	simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Cluster process\cr
	\code{\link{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link{rGaussPoisson}}  \tab
	simulate the Gauss-Poisson cluster process\cr
	\code{\link{rCauchy}} \tab
	simulate Neyman-Scott Cauchy cluster process \cr
	\code{\link{rVarGamma}} \tab
	simulate Neyman-Scott Variance Gamma cluster process \cr
	\code{\link{rthin}} \tab  random thinning  \cr
	\code{\link{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link{simulate.ppm}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link{runifpointOnLines}} \tab
	generate \eqn{n} random points along specified line segments \cr
	\code{\link{rpoisppOnLines}} \tab
	generate Poisson random points along specified line segments 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link{rshift}} \tab random shifting of points \cr
        \code{\link{rjitter}} \tab
	apply random displacements to points in a pattern\cr
	\code{\link{rthin}} \tab  random thinning \cr
	\code{\link{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
        \code{\link{quadratresample}} \tab block resampling 
      }

      \bold{Standard point pattern datasets:}

      Datasets in \pkg{spatstat} are lazy-loaded, so you can simply
      type the name of the dataset to use it; there is no need
      to type \code{\link{data}(amacrine)} etc.

      Type \code{demo(data)} to see a display of all the datasets
      installed with the package.

      Type \code{vignette('datasets')} for a document giving an overview
      of all datasets, including background information, and plots.
      
      \tabular{ll}{
	\code{\link{amacrine}} \tab Austin Hughes' rabbit amacrine cells \cr
	\code{\link{anemones}} \tab Upton-Fingleton sea anemones data\cr
	\code{\link{ants}} \tab Harkness-Isham ant nests data\cr
	\code{\link{bdspots}} \tab Breakdown spots in microelectrodes \cr
	\code{\link{bei}} \tab Tropical rainforest trees \cr
	\code{\link{betacells}} \tab Waessle et al. cat retinal ganglia data \cr
	\code{\link{bramblecanes}} \tab Bramble Canes data \cr
	\code{\link{bronzefilter}} \tab Bronze Filter Section data \cr
	\code{\link{cells}} \tab Crick-Ripley biological cells data \cr
	\code{\link{chicago}} \tab Chicago street crimes \cr
	\code{\link{chorley}} \tab Chorley-Ribble cancer data \cr
	\code{\link{clmfires}} \tab Castilla-La Mancha forest fires \cr
	\code{\link{copper}} \tab Berman-Huntington copper deposits data \cr
	\code{\link{dendrite}} \tab Dendritic spines \cr
	\code{\link{demohyper}} \tab Synthetic point patterns\cr
	\code{\link{demopat}} \tab Synthetic point pattern \cr
	\code{\link{finpines}} \tab Finnish Pines data \cr
	\code{\link{flu}} \tab Influenza virus proteins \cr
	\code{\link{gordon}} \tab People in Gordon Square, London \cr
	\code{\link{gorillas}} \tab Gorilla nest sites \cr
	\code{\link{hamster}} \tab Aherne's hamster tumour data \cr
	\code{\link{humberside}} \tab North Humberside childhood leukaemia data \cr
	\code{\link{hyytiala}} \tab {Mixed forest in
              \ifelse{latex}{\out{Hyyti{\"a}l{\"a}}}{Hyytiala}, Finland}\cr
	\code{\link{japanesepines}} \tab Japanese Pines data \cr
	\code{\link{lansing}} \tab Lansing Woods data \cr
	\code{\link{longleaf}} \tab Longleaf Pines data \cr
	\code{\link{mucosa}} \tab Cells in gastric mucosa \cr
	\code{\link{murchison}} \tab Murchison gold deposits \cr
	\code{\link{nbfires}} \tab New Brunswick fires data \cr
	\code{\link{nztrees}} \tab Mark-Esler-Ripley trees data \cr
	\code{\link{osteo}} \tab Osteocyte lacunae (3D, replicated) \cr
	\code{\link{paracou}} \tab Kimboto trees in Paracou, French Guiana \cr
	\code{\link{ponderosa}} \tab Getis-Franklin ponderosa pine trees data \cr
	\code{\link{pyramidal}} \tab Pyramidal neurons from 31 brains \cr
	\code{\link{redwood}} \tab Strauss-Ripley redwood saplings data \cr
	\code{\link{redwoodfull}} \tab Strauss redwood saplings data (full set) \cr
	\code{\link{residualspaper}} \tab Data from Baddeley et al (2005) \cr
	\code{\link{shapley}} \tab Galaxies in an astronomical survey \cr
	\code{\link{simdat}} \tab Simulated point pattern (inhomogeneous, with interaction) \cr
	\code{\link{spiders}} \tab Spider webs on mortar lines of brick wall \cr
	\code{\link{sporophores}} \tab Mycorrhizal fungi around a tree \cr
	\code{\link{spruces}} \tab Spruce trees in Saxonia \cr
	\code{\link{swedishpines}} \tab Strand-Ripley Swedish pines data \cr
	\code{\link{urkiola}} \tab Urkiola Woods data \cr
	\code{\link{waka}} \tab Trees in Waka national park \cr
	\code{\link{waterstriders}} \tab Insects on water surface 
      }

      \bold{To manipulate a point pattern:}

      \tabular{ll}{
	\code{\link{plot.ppp}} \tab
	plot a point pattern (e.g. \code{plot(X)}) \cr
	\code{\link{iplot}} \tab
	plot a point pattern interactively \cr
	\code{\link{edit.ppp}} \tab interactive text editor \cr
	\code{\link{[.ppp}} \tab
	extract or replace a subset of a point pattern \cr
        \tab \code{pp[subset]} or \code{pp[subwindow]} \cr
	\code{\link{subset.ppp}} \tab
	extract subset of point pattern satisfying a condition \cr
	\code{\link{superimpose}} \tab
	combine several point patterns  \cr
	\code{\link{by.ppp}} \tab
	apply a function to sub-patterns of a point pattern \cr
	\code{\link{cut.ppp}} \tab
	classify the points in a point pattern \cr
	\code{\link{split.ppp}} \tab
	divide pattern into sub-patterns \cr
	\code{\link{unmark}} \tab
	remove marks  \cr
	\code{\link{npoints}} \tab
	count the number of points  \cr
	\code{\link{coords}} \tab
	extract coordinates, change coordinates  \cr
	\code{\link{marks}} \tab
	extract marks, change marks or attach marks  \cr
	\code{\link{rotate}} \tab
	rotate pattern  \cr
	\code{\link{shift} } \tab
	translate pattern  \cr
	\code{\link{flipxy} } \tab
	swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link{reflect} } \tab
	reflect in the origin  \cr
	\code{\link{periodify} } \tab
	make several translated copies  \cr
	\code{\link{affine}} \tab
	apply affine transformation\cr
	\code{\link{scalardilate}} \tab
	apply scalar dilation\cr
	\code{\link{density.ppp}} \tab
	kernel estimation of point pattern intensity\cr
	\code{\link{Smooth.ppp}} \tab
	kernel smoothing of marks of point pattern\cr
	\code{\link{nnmark}} \tab
	mark value of nearest data point\cr
	\code{\link{sharpen.ppp}} \tab
	data sharpening\cr
	\code{\link{identify.ppp}} \tab
	interactively identify points \cr
	\code{\link{unique.ppp}} \tab
	remove duplicate points \cr
	\code{\link{duplicated.ppp}} \tab
	determine which points are duplicates \cr
	\code{\link{connected.ppp}} \tab find clumps of points \cr
	\code{\link{dirichlet}} \tab
	compute Dirichlet-Voronoi tessellation \cr
	\code{\link{delaunay}} \tab
	compute Delaunay triangulation \cr 
	\code{\link{delaunayDistance}} \tab
	graph distance in Delaunay triangulation \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{discretise}} \tab discretise coordinates \cr
	\code{\link{pixellate.ppp}} \tab approximate point pattern by 
	pixel image \cr
	\code{\link{as.im.ppp}} \tab approximate point pattern by 
	pixel image 
      }
      See \code{\link{spatstat.options}} to control plotting behaviour.
      
    \bold{To create a window:}

    An object of class \code{"owin"} describes a spatial region
      (a window of observation).

      \tabular{ll}{
	\code{\link{owin}}	\tab	Create a window object \cr
	\tab \code{owin(xlim, ylim)} for rectangular window \cr
	\tab \code{owin(poly)} for polygonal window \cr
	\tab \code{owin(mask)} for binary image window \cr
	\code{\link{Window}}	\tab
	Extract window of another object \cr
	\code{\link{Frame}}	\tab
	Extract the containing rectangle ('frame') of another object \cr
	\code{\link{as.owin}}	\tab
	Convert other data to a window object \cr
	\code{\link{square}}    \tab	make a square window \cr
	\code{\link{disc}}    \tab	make a circular window \cr
	\code{\link{ellipse}}    \tab	make an elliptical window \cr
	\code{\link{ripras}}    \tab
	Ripley-Rasson estimator of window, given only the points \cr
	\code{\link{convexhull}} \tab compute convex hull of something \cr 
	\code{\link{letterR}}    \tab
	polygonal window in the shape of the \R logo \cr
	\code{\link{clickpoly}}    \tab
	interactively draw a polygonal window  \cr
	\code{\link{clickbox}}    \tab
	interactively draw a rectangle  
      }

    \bold{To manipulate a window:}

    \tabular{ll}{
	\code{\link{plot.owin}}	\tab	plot a window. \cr
	\tab		\code{plot(W)}\cr
	\code{\link{boundingbox}} \tab
	Find a tight bounding box for the window \cr
	\code{\link{erosion}}	\tab
	erode window by a distance r\cr
	\code{\link{dilation}}	\tab
	dilate window by a distance r\cr
	\code{\link{closing}}	\tab
	close window by a distance r\cr
	\code{\link{opening}}	\tab
	open window by a distance r\cr
	\code{\link{border}}	\tab
	difference between window and its erosion/dilation \cr 
	\code{\link{complement.owin}}	\tab
	invert (swap inside and outside)\cr
	\code{\link{simplify.owin}}	\tab
	approximate a window by a simple polygon  \cr
	\code{\link{rotate}} \tab rotate window  \cr
	\code{\link{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link{shift} } \tab translate window  \cr
	\code{\link{periodify} } \tab make several translated copies  \cr
	\code{\link{affine}} \tab apply affine transformation \cr
	\code{\link{as.data.frame.owin}} \tab
	convert window to data frame 
      }

    \bold{Digital approximations:}

    \tabular{ll}{
	\code{\link{as.mask}}	\tab
	Make a discrete pixel approximation of a given window \cr
	\code{\link{as.im.owin}} \tab convert window to pixel image \cr
	\code{\link{pixellate.owin}} \tab convert window to pixel image \cr
	\code{\link{commonGrid}} \tab find common pixel grid for windows \cr
	\code{\link{nearest.raster.point}} \tab
	map continuous coordinates to raster locations\cr
	\code{\link{raster.x}} \tab
	raster x coordinates \cr
	\code{\link{raster.y}} \tab
	raster y coordinates \cr
	\code{\link{raster.xy}} \tab
	raster x and y coordinates \cr
	\code{\link{as.polygonal}} \tab
	convert pixel mask to polygonal window 
      }
      See \code{\link{spatstat.options}} to control the approximation

    \bold{Geometrical computations with windows:}

    \tabular{ll}{
	\code{\link{edges}}	\tab	extract boundary edges \cr
	\code{\link{intersect.owin}}	\tab	intersection of two windows\cr
	\code{\link{union.owin}}	\tab	union of two windows\cr
	\code{\link{setminus.owin}}	\tab	set subtraction of two windows\cr
	\code{\link{inside.owin}}	\tab	determine whether a point is inside a window\cr
	\code{\link{area.owin}}	\tab	compute area \cr
	\code{\link{perimeter}}	\tab	compute perimeter length \cr
	\code{\link{diameter.owin}}	\tab	compute diameter\cr
	\code{\link{incircle}}	\tab	find largest circle inside a window \cr
	\code{\link{inradius}}	\tab	radius of incircle \cr
	\code{\link{connected.owin}}    \tab find connected components of window \cr
	\code{\link{eroded.areas}}	\tab	compute areas of eroded windows\cr
	\code{\link{dilated.areas}}	\tab	compute areas of dilated windows\cr
	\code{\link{bdist.points}}	\tab	compute distances from data points to window boundary \cr
	\code{\link{bdist.pixels}}	\tab	compute distances from all pixels to window boundary \cr
	\code{\link{bdist.tiles}}	\tab
	boundary distance for each tile in tessellation \cr
	\code{\link{distmap.owin}}	\tab	distance transform image \cr
	\code{\link{distfun.owin}}	\tab	distance transform \cr
	\code{\link{centroid.owin}}	\tab compute centroid (centre of mass) of window\cr
	\code{\link{is.subset.owin}}    \tab determine whether one
	window contains another \cr
	\code{\link{is.convex}} \tab determine whether a window is convex \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{triangulate.owin}} \tab decompose into triangles \cr
	\code{\link{as.mask}} \tab pixel approximation of window \cr
	\code{\link{as.polygonal}} \tab polygonal approximation of window \cr
	\code{\link{is.rectangle}} \tab test whether window is a rectangle \cr 
	\code{\link{is.polygonal}} \tab test whether window is polygonal \cr 
	\code{\link{is.mask}} \tab test whether window is a mask \cr 
	\code{\link{setcov}} \tab spatial covariance function of window \cr
	\code{\link{pixelcentres}} \tab extract centres of pixels in mask \cr
	\code{\link{clickdist}}    \tab
	measure distance between two points clicked by user
      }

    \bold{Pixel images:}
    An object of class \code{"im"} represents a pixel image. 
    Such objects are returned by some of the functions in
      \pkg{spatstat} including \code{\link{Kmeasure}},
      \code{\link{setcov}} and \code{\link{density.ppp}}. 
      \tabular{ll}{
	\code{\link{im}} \tab create a pixel image\cr
	\code{\link{as.im}} \tab convert other data to a pixel image\cr
	\code{\link{pixellate}} \tab convert other data to a pixel image\cr
	\code{\link{as.matrix.im}} \tab convert pixel image to matrix\cr
	\code{\link{as.data.frame.im}} \tab convert pixel image to data frame\cr
	\code{\link{as.function.im}} \tab convert pixel image to function\cr
	\code{\link{plot.im}}	\tab	plot a pixel image on screen as a digital image\cr
	\code{\link{contour.im}}	\tab draw contours of a pixel image \cr
	\code{\link{persp.im}}	\tab draw perspective plot of a pixel image \cr
	\code{\link{rgbim}}	\tab create colour-valued pixel image \cr
	\code{\link{hsvim}}	\tab create colour-valued pixel image \cr
	\code{\link{[.im}} 	\tab extract a subset of a pixel image\cr
	\code{\link{[<-.im}} 	\tab replace a subset of a pixel image\cr
	\code{\link{rotate.im}} \tab rotate pixel image \cr
	\code{\link{shift.im}} \tab apply vector shift to pixel image \cr
	\code{\link{affine.im}} \tab apply affine transformation to image \cr
	\code{X}	\tab	print very basic information about image \code{X}\cr
	\code{\link{summary}(X)} \tab	summary of image \code{X} \cr
	\code{\link{hist.im}} \tab	histogram of image \cr
	\code{\link{mean.im}} \tab	mean pixel value of image  \cr
	\code{\link{integral.im}} \tab	integral of pixel values  \cr
	\code{\link{quantile.im}} \tab	quantiles of image  \cr
	\code{\link{cut.im}} \tab	convert numeric image to factor image \cr
	\code{\link{is.im}} \tab test whether an object is a pixel image\cr
	\code{\link{interp.im}} \tab interpolate a pixel image\cr
	\code{\link{blur}} \tab apply Gaussian blur to image\cr
	\code{\link{Smooth.im}} \tab apply Gaussian blur to image\cr
	\code{\link{connected.im}} \tab find connected components \cr
	\code{\link{compatible.im}} \tab test whether two images have
	compatible dimensions \cr
	\code{\link{harmonise.im}} \tab make images compatible \cr
	\code{\link{commonGrid}} \tab find a common pixel grid for images \cr
	\code{\link{eval.im}} \tab evaluate any expression involving images\cr
	\code{\link{scaletointerval}} \tab rescale pixel values \cr
	\code{\link{zapsmall.im}} \tab set very small pixel values to zero \cr
	\code{\link{levelset}} \tab level set of an image\cr 
	\code{\link{solutionset}} \tab region where an expression is true \cr
	\code{\link{imcov}} \tab spatial covariance function of image \cr
	\code{\link{convolve.im}} \tab spatial convolution of images \cr
	\code{\link{transect.im}} \tab line transect of image \cr
	\code{\link{pixelcentres}} \tab extract centres of pixels \cr
	\code{\link{transmat}} \tab convert matrix of pixel values \cr
	                       \tab to a different indexing convention \cr
	\code{\link{rnoise}} \tab random pixel noise
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.
    \tabular{ll}{
      \code{\link{psp}} \tab create a line segment pattern \cr
      \code{\link{as.psp}} \tab convert other data into a line segment pattern \cr
      \code{\link{edges}} \tab extract edges of a window \cr
      \code{\link{is.psp}} \tab determine whether a dataset has class \code{"psp"} \cr
      \code{\link{plot.psp}} \tab plot a line segment pattern \cr
      \code{\link{print.psp}} \tab print basic information \cr
      \code{\link{summary.psp}} \tab print summary information \cr
      \code{\link{[.psp}} \tab extract a subset of a line segment pattern \cr
      \code{\link{as.data.frame.psp}} \tab
      convert line segment pattern to data frame \cr
      \code{\link{marks.psp}} \tab extract marks of line segments \cr
      \code{\link{marks<-.psp}} \tab assign new marks to line segments \cr
      \code{\link{unmark.psp}} \tab delete marks from line segments \cr
      \code{\link{midpoints.psp}} \tab
      compute the midpoints of line segments \cr
      \code{\link{endpoints.psp}} \tab extract the endpoints of line segments \cr
      \code{\link{lengths.psp}} \tab compute the lengths of line segments \cr
      \code{\link{angles.psp}} \tab compute the orientation angles of line segments \cr
      \code{\link{superimpose}} \tab combine several line segment patterns  \cr
      \code{\link{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates \cr      
      \code{\link{rotate.psp}} \tab rotate a line segment pattern \cr      
      \code{\link{shift.psp}} \tab shift a line segment pattern \cr      
      \code{\link{periodify}} \tab make several shifted copies \cr      
      \code{\link{affine.psp}} \tab apply an affine transformation \cr      
      \code{\link{pixellate.psp}} \tab approximate line segment pattern
    by pixel image \cr      
      \code{\link{as.mask.psp}} \tab approximate line segment pattern
    by binary mask \cr      
      \code{\link{distmap.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{distfun.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{density.psp}} \tab kernel smoothing of line segments\cr
      \code{\link{selfcrossing.psp}} \tab find crossing points between
      line segments \cr
      \code{\link{selfcut.psp}} \tab cut segments where they cross \cr
      \code{\link{crossing.psp}} \tab find crossing points between
      two line segment patterns \cr
      \code{\link{nncross}} \tab find distance to nearest line segment
      from a given point\cr
      \code{\link{nearestsegment}} \tab find line segment closest to a
      given point \cr
      \code{\link{project2segment}} \tab find location along a line segment
      closest to a given point \cr
      \code{\link{pointsOnLines}} \tab generate points evenly spaced
      along line segment \cr
      \code{\link{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window\cr
      \code{\link{rlinegrid}} \tab generate a random array of parallel
      lines through a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link{tess}} \tab create a tessellation \cr
      \code{\link{quadrats}} \tab create a tessellation of rectangles\cr
      \code{\link{hextess}} \tab create a tessellation of hexagons \cr
      \code{\link{quantess}} \tab quantile tessellation \cr
      \code{\link{as.tess}} \tab convert other data to a tessellation \cr
      \code{\link{plot.tess}} \tab plot a tessellation \cr
      \code{\link{tiles}} \tab extract all the tiles of a tessellation \cr
      \code{\link{[.tess}} \tab extract some tiles of a tessellation \cr
      \code{\link{[<-.tess}} \tab change some tiles of a tessellation \cr
      \code{\link{intersect.tess}} \tab intersect two tessellations \cr
                              \tab or restrict a tessellation to a window \cr
      \code{\link{chop.tess}} \tab subdivide a tessellation by a line \cr
      \code{\link{dirichlet}} \tab compute Dirichlet-Voronoi tessellation of points\cr
      \code{\link{delaunay}} \tab compute Delaunay triangulation of points\cr
      \code{\link{rpoislinetess}} \tab generate tessellation using Poisson line
      process \cr
      \code{\link{tile.areas}}	\tab
      area of each tile in tessellation \cr
      \code{\link{bdist.tiles}}	\tab
      boundary distance for each tile in tessellation 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link{pp3}} \tab create a 3-D point pattern \cr
      \code{\link{plot.pp3}} \tab plot a 3-D point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{subset.pp3}} \tab extract subset of 3-D point pattern \cr
      \code{\link{unitname.pp3}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link{rpoispp3}} \tab generate Poisson random points in 3-D \cr
      \code{\link{envelope.pp3}} \tab generate simulation envelopes for
      3-D pattern \cr
      \code{\link{box3}} \tab create a 3-D rectangular box \cr
      \code{\link{as.box3}} \tab convert data to 3-D rectangular box \cr
      \code{\link{unitname.box3}} \tab name of unit of length \cr
      \code{\link{diameter.box3}} \tab diameter of box \cr
      \code{\link{volume.box3}} \tab volume of box \cr
      \code{\link{shortside.box3}} \tab shortest side of box \cr
      \code{\link{eroded.volumes}} \tab volumes of erosions of box 
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link{ppx}} \tab create a multidimensional space-time point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{subset.ppx}} \tab extract subset \cr
      \code{\link{unitname.ppx}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{runifpointx}} \tab generate uniform random points \cr
      \code{\link{rpoisppx}} \tab generate Poisson random points \cr
      \code{\link{boxx}} \tab define multidimensional box  \cr
      \code{\link{diameter.boxx}} \tab diameter of box \cr
      \code{\link{volume.boxx}} \tab volume of box \cr
      \code{\link{shortside.boxx}} \tab shortest side of box \cr
      \code{\link{eroded.volumes.boxx}} \tab volumes of erosions of box 
    }
    
    \bold{Point patterns on a linear network}

    An object of class \code{"linnet"} represents a linear network
    (for example, a road network).

    \tabular{ll}{
      \code{\link{linnet}} \tab create a linear network \cr
      \code{\link{clickjoin}} \tab interactively join vertices in network \cr
      \code{\link{iplot.linnet}} \tab interactively plot network \cr
      \code{\link{simplenet}} \tab simple example of network \cr
      \code{\link{lineardisc}} \tab disc in a linear network \cr
      \code{\link{delaunayNetwork}} \tab network of Delaunay triangulation \cr
      \code{\link{dirichletNetwork}} \tab network of Dirichlet edges \cr
      \code{\link{methods.linnet}} \tab methods for \code{linnet} objects\cr
      \code{\link{vertices.linnet}} \tab nodes of network \cr 
      \code{\link{pixellate.linnet}} \tab approximate by pixel image
    }
    
    An object of class \code{"lpp"} represents a 
    point pattern on a linear network (for example,
    road accidents on a road network).
    
    \tabular{ll}{
      \code{\link{lpp}} \tab create a point pattern on a linear network \cr
      \code{\link{methods.lpp}} \tab methods for \code{lpp} objects \cr
      \code{\link{subset.lpp}} \tab method for \code{subset} \cr
      \code{\link{rpoislpp}} \tab simulate Poisson points on linear network \cr
      \code{\link{runiflpp}} \tab simulate random points on a linear network \cr
      \code{\link{chicago}} \tab Chicago street crime data \cr
      \code{\link{dendrite}} \tab Dendritic spines data \cr
      \code{\link{spiders}} \tab Spider webs on mortar lines of brick wall 
    }
    
    \bold{Hyperframes}

    A hyperframe is like a data frame, except that the entries
    may be objects of any kind.

    \tabular{ll}{
      \code{\link{hyperframe}} \tab create a hyperframe \cr
      \code{\link{as.hyperframe}} \tab convert data to hyperframe \cr
      \code{\link{plot.hyperframe}} \tab plot hyperframe \cr
      \code{\link{with.hyperframe}} \tab evaluate expression using each row
      of hyperframe \cr
      \code{\link{cbind.hyperframe}} \tab combine hyperframes by columns\cr
      \code{\link{rbind.hyperframe}} \tab combine hyperframes by rows\cr
      \code{\link{as.data.frame.hyperframe}} \tab convert hyperframe to
      data frame \cr
      \code{\link{subset.hyperframe}} \tab method for \code{subset} \cr
      \code{\link{head.hyperframe}} \tab first few rows of hyperframe \cr
      \code{\link{tail.hyperframe}} \tab last few rows of hyperframe
    }
    
    \bold{Layered objects}

    A layered object represents data that should be plotted in
    successive layers, for example, a background and a foreground.

     \tabular{ll}{
      \code{\link{layered}} \tab create layered object \cr
      \code{\link{plot.layered}} \tab plot layered object\cr
      \code{\link{[.layered}} \tab extract subset of layered object
    }
    
    \bold{Colour maps}

    A colour map is a mechanism for associating colours with data.
    It can be regarded as a function, mapping data to colours.
    Using a \code{colourmap} object in a plot command
    ensures that the mapping from numbers to colours is
    the same in different plots. 

     \tabular{ll}{
      \code{\link{colourmap}} \tab create a colour map  \cr
      \code{\link{plot.colourmap}} \tab plot the colour map only\cr
      \code{\link{tweak.colourmap}} \tab alter individual colour values \cr
      \code{\link{interp.colourmap}} \tab make a smooth transition
      between colours \cr
      \code{\link{beachcolourmap}} \tab one special colour map
    }
}
  
\section{II. EXPLORATORY DATA ANALYSIS}{

  \bold{Inspection of data:}
  \tabular{ll}{
    \code{\link{summary}(X)} \tab
    print useful summary of point pattern \code{X}\cr
    \code{X} \tab
    print basic description of point pattern \code{X}  \cr
    \code{any(duplicated(X))} \tab
    check for duplicated points in pattern \code{X} \cr
    \code{\link{istat}(X)} \tab Interactive exploratory analysis \cr
    \code{\link{View}(X)} \tab spreadsheet-style viewer
  }

  \bold{Classical exploratory tools:}
  \tabular{ll}{
    \code{\link{clarkevans}} \tab Clark and Evans aggregation index \cr
    \code{\link{fryplot}} \tab Fry plot \cr
    \code{\link{miplot}} \tab Morisita Index plot
  }

  \bold{Smoothing:}
  \tabular{ll}{
    \code{\link{density.ppp}} \tab kernel smoothed density/intensity\cr
    \code{\link{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link{bw.diggle}} \tab cross-validated bandwidth selection
                               for \code{\link{density.ppp}}\cr
    \code{\link{bw.ppl}} \tab likelihood cross-validated bandwidth selection
                               for \code{\link{density.ppp}}\cr
    \code{\link{bw.scott}} \tab Scott's rule of thumb 
                                for density estimation\cr
    \code{\link{bw.relrisk}} \tab cross-validated bandwidth selection
                               for \code{\link{relrisk}} \cr
    \code{\link{bw.smoothppp}} \tab cross-validated bandwidth selection
                               for \code{\link{Smooth.ppp}} \cr
    \code{\link{bw.frac}} \tab bandwidth selection using window geometry\cr
    \code{\link{bw.stoyan}} \tab Stoyan's rule of thumb for bandwidth
                               for \code{\link{pcf}}
  }

  \bold{Modern exploratory tools:}
  \tabular{ll}{
    \code{\link{clusterset}} \tab Allard-Fraley feature detection  \cr
    \code{\link{nnclean}} \tab Byers-Raftery feature detection  \cr
    \code{\link{sharpen.ppp}} \tab Choi-Hall data sharpening \cr
    \code{\link{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link{rho2hat}} \tab Kernel estimate of effect of two covariates\cr
    \code{\link{spatialcdf}} \tab Spatial cumulative distribution function\cr
    \code{\link{roc}} \tab Receiver operating characteristic curve
  }

  \bold{Summary statistics for a point pattern:}
  Type \code{demo(sumfun)} for a demonstration of many
  of the summary statistics.
  \tabular{ll}{
    \code{\link{intensity}} \tab Mean intensity \cr
    \code{\link{quadratcount}} \tab Quadrat counts \cr
    \code{\link{intensity.quadratcount}} \tab Mean intensity in quadrats \cr
    \code{\link{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link{Lest}} \tab Besag \eqn{L}-function\cr
    \code{\link{Tstat}} \tab Third order \eqn{T}-function \cr
    \code{\link{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link{pcf}} \tab 	pair correlation function \cr
    \code{\link{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link{Linhom}} \tab \eqn{L} for inhomogeneous point patterns \cr
    \code{\link{pcfinhom}} \tab pair correlation for inhomogeneous patterns\cr
    \code{\link{Finhom}} \tab \eqn{F} for inhomogeneous point patterns \cr
    \code{\link{Ginhom}} \tab \eqn{G} for inhomogeneous point patterns \cr
    \code{\link{Jinhom}} \tab \eqn{J} for inhomogeneous point patterns \cr
    \code{\link{localL}} \tab Getis-Franklin neighbourhood density function\cr
    \code{\link{localK}} \tab neighbourhood K-function\cr
    \code{\link{localpcf}} \tab local pair correlation function\cr
    \code{\link{localKinhom}} \tab local \eqn{K} for inhomogeneous point patterns \cr
    \code{\link{localLinhom}} \tab local \eqn{L} for inhomogeneous point patterns \cr
    \code{\link{localpcfinhom}} \tab local pair correlation for inhomogeneous patterns\cr
    \code{\link{Ksector}} \tab Directional \eqn{K}-function\cr
    \code{\link{Kscaled}} \tab locally scaled \eqn{K}-function \cr
    \code{\link{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link{envelope}} \tab simulation envelopes for a summary
    function \cr
    \code{\link{varblock}} \tab variances and confidence intervals\cr
                   \tab for a summary function \cr
    \code{\link{lohboot}} \tab bootstrap for a summary function 
  }

  Related facilities:
  \tabular{ll}{
    \code{\link{plot.fv}} \tab plot a summary function\cr
    \code{\link{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link{harmonise.fv}} \tab make functions compatible \cr
    \code{\link{eval.fasp}} \tab evaluate any expression involving
    an array of functions\cr
    \code{\link{with.fv}} \tab evaluate an expression for a 
    summary function\cr
    \code{\link{Smooth.fv}} \tab apply smoothing to a summary function\cr
    \code{\link{deriv.fv}} \tab calculate derivative of a summary function\cr
    \code{\link{pool.fv}} \tab pool several estimates of a summary function\cr
    \code{\link{nndist}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich}} \tab find nearest neighbours \cr
    \code{\link{pairdist}} \tab distances between all pairs of points\cr
    \code{\link{crossdist}} \tab distances between points in two patterns\cr
    \code{\link{nncross}} \tab nearest neighbours between two point patterns \cr
    \code{\link{exactdt}} \tab distance from any location to nearest data point\cr 
    \code{\link{distmap}} \tab distance map image\cr
    \code{\link{distfun}} \tab distance map function\cr
    \code{\link{nnmap}} \tab nearest point image \cr
    \code{\link{nnfun}} \tab nearest point function \cr
    \code{\link{density.ppp}} \tab kernel smoothed density\cr
    \code{\link{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link{sharpen.ppp}} \tab data sharpening  \cr
    \code{\link{rknn}} \tab theoretical distribution of nearest
    neighbour distance
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} such that \code{marks(X)} is a factor. 
  \tabular{ll}{
    \code{\link{relrisk}} \tab kernel estimation of relative risk  \cr
    \code{\link{scan.test}} \tab spatial scan test of elevated risk  \cr
    \code{\link{Gcross},\link{Gdot},\link{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{G[i,j], G[i.]} \cr
    \code{\link{Kcross},\link{Kdot}, \link{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{K[i,j], K[i.]} \cr
    \code{\link{Lcross},\link{Ldot}} \tab
    multitype \eqn{L}-functions 
    \eqn{L_{ij}, L_{i\bullet}}{L[i,j], L[i.]} \cr
    \code{\link{Jcross},\link{Jdot},\link{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{J[i,j],J[i.]} \cr
    \code{\link{pcfcross}} \tab
    multitype pair correlation function \eqn{g_{ij}}{g[i,j]} \cr
    \code{\link{pcfdot}} \tab
    multitype pair correlation function \eqn{g_{i\bullet}}{g[i.]} \cr
    \code{\link{pcfmulti}} \tab
    general pair correlation function \cr
    \code{\link{markconnect}} \tab
    marked connection function \eqn{p_{ij}}{p[i,j]} \cr
    \code{\link{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link{Iest}} \tab  multitype \eqn{I}-function\cr
    \code{\link{Kcross.inhom},\link{Kdot.inhom}} \tab
    inhomogeneous counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link{Lcross.inhom},\link{Ldot.inhom}} \tab
    inhomogeneous counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link{pcfcross.inhom},\link{pcfdot.inhom}} \tab
    inhomogeneous counterparts of \code{pcfcross}, \code{pcfdot} 
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type. For numeric marks, there are the
  following functions:
  \tabular{ll}{
    \code{\link{markmean}} \tab smoothed local average of marks \cr
    \code{\link{markvar}} \tab smoothed local variance of marks \cr
    \code{\link{markcorr}} \tab mark correlation function \cr
    \code{\link{markcrosscorr}} \tab mark cross-correlation function \cr
    \code{\link{markvario}} \tab mark variogram \cr
    \code{\link{Kmark}} \tab mark-weighted \eqn{K} function \cr
    \code{\link{Emark}} \tab mark independence diagnostic \eqn{E(r)} \cr
    \code{\link{Vmark}} \tab mark independence diagnostic \eqn{V(r)} \cr
    \code{\link{nnmean}} \tab nearest neighbour mean index \cr
    \code{\link{nnvario}} \tab nearest neighbour mark variance index 
  }
  For marks of any type, there are the following:
  \tabular{ll}{
    \code{\link{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link{applynbd}} \tab apply function to every neighbourhood
    in a point pattern \cr
    \code{\link{markstat}} \tab apply function to the marks of neighbours
    in a point pattern \cr
    \code{\link{marktable}} \tab tabulate the marks of neighbours
    in a point pattern \cr
    \code{\link{pppdist}} \tab find the optimal match between two point
    patterns
  }

  \bold{Summary statistics for a point pattern on a linear network:}

  These are for point patterns on a linear network (class \code{lpp}).
  For unmarked patterns:
  
  \tabular{ll}{
    \code{\link{linearK}} \tab
    \eqn{K} function on linear network \cr
    \code{\link{linearKinhom}} \tab
    inhomogeneous \eqn{K} function on linear network \cr
    \code{\link{linearpcf}} \tab
    pair correlation function on linear network \cr
    \code{\link{linearpcfinhom}} \tab
    inhomogeneous pair correlation on linear network
  }

  For multitype patterns:
  \tabular{ll}{
    \code{\link{linearKcross}} \tab
    \eqn{K} function between two types of points \cr
    \code{\link{linearKdot}} \tab
    \eqn{K} function from one type to any type \cr
    \code{\link{linearKcross.inhom}} \tab
    Inhomogeneous version of \code{\link{linearKcross}} \cr
    \code{\link{linearKdot.inhom}} \tab
    Inhomogeneous version of \code{\link{linearKdot}} \cr
    \code{\link{linearmarkconnect}} \tab
    Mark connection function  on linear network \cr
    \code{\link{linearmarkequal}} \tab
    Mark equality function on linear network \cr
    \code{\link{linearpcfcross}} \tab
    Pair correlation between two types of points \cr
    \code{\link{linearpcfdot}} \tab
    Pair correlation from one type to any type \cr
    \code{\link{linearpcfcross.inhom}} \tab
    Inhomogeneous version of \code{\link{linearpcfcross}} \cr
    \code{\link{linearpcfdot.inhom}} \tab
    Inhomogeneous version of \code{\link{linearpcfdot}} 
  }

  Related facilities:
  
  \tabular{ll}{
    \code{\link{pairdist.lpp}} \tab distances between pairs  \cr
    \code{\link{crossdist.lpp}} \tab distances between pairs \cr
    \code{\link{nndist.lpp}} \tab nearest neighbour distances  \cr
    \code{\link{nncross.lpp}} \tab nearest neighbour distances  \cr
    \code{\link{nnwhich.lpp}} \tab find nearest neighbours  \cr
    \code{\link{nnfun.lpp}} \tab find nearest data point  \cr
    \code{\link{density.lpp}} \tab kernel smoothing estimator of intensity  \cr
    \code{\link{distfun.lpp}} \tab distance transform  \cr
    \code{\link{envelope.lpp}} \tab simulation envelopes  \cr
    \code{\link{rpoislpp}} \tab simulate Poisson points on linear network \cr
    \code{\link{runiflpp}} \tab simulate random points on a linear network 
  }
  
  It is also possible to fit point process models to \code{lpp} objects.
  See Section IV.
  
  \bold{Summary statistics for a three-dimensional point pattern:}

  These are for 3-dimensional point pattern objects (class \code{pp3}).

  \tabular{ll}{
    \code{\link{F3est}} \tab empty space function \eqn{F} \cr
    \code{\link{G3est}} \tab nearest neighbour function \eqn{G} \cr
    \code{\link{K3est}} \tab \eqn{K}-function \cr
    \code{\link{pcf3est}} \tab pair correlation function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link{envelope.pp3}} \tab simulation envelopes \cr
    \code{\link{pairdist.pp3}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.pp3}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.pp3}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.pp3}} \tab find nearest neighbours \cr
    \code{\link{nncross.pp3}} \tab find nearest neighbours in another pattern
  }

  \bold{Computations for multi-dimensional point pattern:}

  These are for multi-dimensional space-time
  point pattern objects (class \code{ppx}).

  \tabular{ll}{
    \code{\link{pairdist.ppx}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.ppx}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.ppx}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.ppx}} \tab find nearest neighbours
  }

  \bold{Summary statistics for random sets:}
  
  These work for point patterns (class \code{ppp}),
  line segment patterns (class \code{psp})
  or windows (class \code{owin}).
  
  \tabular{ll}{
    \code{\link{Hest}} \tab spherical contact distribution \eqn{H} \cr
    \code{\link{Gfox}} \tab Foxall \eqn{G}-function \cr
    \code{\link{Jfox}} \tab Foxall \eqn{J}-function
  }
  
}

\section{III. MODEL FITTING (COX AND CLUSTER MODELS)}{
  
  Cluster process models (with homogeneous or inhomogeneous intensity)
  and Cox processes can be fitted by the function \code{\link{kppm}}.
  Its result is an object of class \code{"kppm"}.
  The fitted model can be printed, plotted, predicted, simulated
  and updated.

  \tabular{ll}{
    \code{\link{kppm}} \tab  Fit model\cr
    \code{\link{plot.kppm}} \tab  Plot the fitted model\cr
    \code{\link{summary.kppm}} \tab  Summarise the fitted model\cr
    \code{\link{fitted.kppm}} \tab Compute fitted intensity \cr
    \code{\link{predict.kppm}} \tab Compute fitted intensity \cr
    \code{\link{update.kppm}} \tab Update the model \cr
    \code{\link{improve.kppm}} \tab Refine the estimate of trend \cr
    \code{\link{simulate.kppm}} \tab Generate simulated realisations \cr
    \code{\link{vcov.kppm}} \tab Variance-covariance matrix of coefficients \cr
    \code{\link[spatstat:methods.kppm]{coef.kppm}}
    \tab Extract trend coefficients \cr
    \code{\link[spatstat:methods.kppm]{formula.kppm}}
    \tab Extract trend formula \cr
    \code{\link{parameters}} \tab Extract all model parameters \cr
    \code{\link{clusterfield}} \tab Compute offspring density \cr
    \code{\link{clusterradius}} \tab Radius of support of offspring density \cr
    \code{\link{Kmodel.kppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link{pcfmodel.kppm}} \tab Pair correlation of fitted model 
  }
  
  For model selection, you can also use 
  the generic functions \code{\link{step}}, \code{\link{drop1}} 
  and \code{\link{AIC}} on fitted point process models.

  The theoretical models can also be simulated,
  for any choice of parameter values,
  using \code{\link{rThomas}}, \code{\link{rMatClust}},
  \code{\link{rCauchy}}, \code{\link{rVarGamma}},
  and \code{\link{rLGCP}}.
  
  Lower-level fitting functions include:

  \tabular{ll}{
    \code{\link{lgcp.estK}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link{lgcp.estpcf}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link{thomas.estK}} \tab fit the Thomas process model \cr
    \code{\link{thomas.estpcf}} \tab fit the Thomas process model \cr
    \code{\link{matclust.estK}} \tab fit the Matern Cluster process model \cr
    \code{\link{matclust.estpcf}} \tab fit the Matern Cluster process model \cr
    \code{\link{cauchy.estK}} \tab fit a Neyman-Scott Cauchy cluster process \cr
   \code{\link{cauchy.estpcf}} \tab fit a Neyman-Scott Cauchy cluster process\cr
   \code{\link{vargamma.estK}} \tab fit a Neyman-Scott Variance Gamma process\cr
   \code{\link{vargamma.estpcf}} \tab fit a Neyman-Scott Variance Gamma process\cr
    \code{\link{mincontrast}} \tab low-level algorithm for fitting models
    \cr \tab by the method of minimum contrast 
  }
}

\section{IV. MODEL FITTING (POISSON AND GIBBS MODELS)}{
  
  \bold{Types of models}
  
  Poisson point processes are the simplest models for point patterns.
  A Poisson model assumes that the points are stochastically
  independent. It may allow the points to have a non-uniform spatial
  density. The special case of a Poisson process with a uniform
  spatial density is often called Complete Spatial Randomness.
  
  Poisson point processes are included in the more general class of Gibbs point
  process models. In a Gibbs model, there is \emph{interaction}
  or dependence between points. Many different types of interaction
  can be specified.
  
  For a detailed explanation of how to fit Poisson or
  Gibbs point process models to point pattern data using \pkg{spatstat},
  see Baddeley and Turner (2005b) or Baddeley (2008).
  
  \bold{To fit a Poisson or Gibbs point process model:}

  Model fitting in \pkg{spatstat} is performed mainly by the function
  \code{\link{ppm}}. Its result is an object of class \code{"ppm"}.
  
  Here are some examples, where \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{ppm(X)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ x)} \tab Poisson process with \cr
                             \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{ppm(X ~ 1, Strauss(0.1))} \tab Stationary Strauss process \cr
    \code{ppm(X ~ x, Strauss(0.1))} \tab Strauss process with \cr
                              \tab conditional intensity loglinear in \eqn{x}
  }
  It is also possible to fit models that depend on
  other covariates.

  \bold{Manipulating the fitted model:}

  \tabular{ll}{
    \code{\link{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link{predict.ppm}}
    \tab   Compute the spatial trend and conditional intensity\cr
    \tab   of the fitted point process model \cr
    \code{\link{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link{parameters}} \tab Extract all model parameters\cr
    \code{\link{formula.ppm}} \tab Extract the trend formula\cr
    \code{\link{intensity.ppm}} \tab Compute fitted intensity \cr
    \code{\link{Kmodel.ppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link{pcfmodel.ppm}} \tab pair correlation of fitted model \cr
    \code{\link{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link{residuals.ppm}} \tab Compute point process residuals at quadrature points \cr
    \code{\link{update.ppm}} \tab Update the fit \cr
    \code{\link{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link{simulate.ppm}} \tab Simulate from fitted model  \cr
    \code{\link{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link{effectfun}} \tab Compute the fitted effect of one covariate\cr
    \code{\link{logLik.ppm}} \tab log-likelihood or log-pseudolikelihood\cr
    \code{\link{anova.ppm}} \tab Analysis of deviance \cr
    \code{\link{model.frame.ppm}} \tab Extract data frame used to fit model  \cr
    \code{\link{model.images}} \tab Extract spatial data used to fit model  \cr
    \code{\link{model.depends}} \tab Identify variables in the model \cr
    \code{\link{as.interact}} \tab Interpoint interaction component of model \cr
    \code{\link{fitin}} \tab Extract fitted interpoint interaction \cr
    \code{\link{is.hybrid}} \tab Determine whether the model is a hybrid \cr
    \code{\link{valid.ppm}} \tab Check the model is a valid point process \cr
    \code{\link{project.ppm}} \tab Ensure the model is a valid point process 
  }
  For model selection, you can also use 
  the generic functions \code{\link{step}}, \code{\link{drop1}} 
  and \code{\link{AIC}} on fitted point process models.
  
  See \code{\link{spatstat.options}} to control plotting of fitted model.
  
  \bold{To specify a point process model:}
  
  The first order ``trend'' of the model is determined by an \R 
  language formula. The formula specifies the form of the
  \emph{logarithm} of the trend.
  
  \tabular{ll}{
    \code{X ~ 1}  \tab No trend (stationary) \cr
    \code{X ~ x}   \tab Loglinear trend
      \eqn{\lambda(x,y) = \exp(\alpha + \beta x)}{lambda(x,y) =	exp(alpha + beta * x)} \cr
    \tab where \eqn{x,y} are Cartesian coordinates \cr
    \code{X ~ polynom(x,y,3)}  \tab Log-cubic polynomial trend  \cr
    \code{X ~ harmonic(x,y,2)}  \tab Log-harmonic polynomial trend \cr
    \code{X ~ Z}   \tab Loglinear function of covariate \code{Z} \cr
    \tab  \eqn{\lambda(x,y) = \exp(\alpha + \beta Z(x,y))}{lambda(x,y) =	exp(alpha + beta * Z(x,y))} 
  }

  The higher order (``interaction'') components are described by
  an object of class \code{"interact"}. Such objects are created by:
  \tabular{ll}{
    \code{\link{Poisson}()} \tab		the Poisson point process\cr
    \code{\link{AreaInter}()}	 \tab Area-interaction process\cr
    \code{\link{BadGey}()} \tab	multiscale Geyer process\cr
    \code{\link{Concom}()} \tab	connected component interaction\cr
    \code{\link{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link{DiggleGatesStibbard}() } \tab Diggle-Gates-Stibbard potential \cr
    \code{\link{Fiksel}()}	 \tab Fiksel pairwise interaction process\cr
    \code{\link{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link{Hardcore}()}	 \tab Hard core process\cr
    \code{\link{HierHard}()}	 \tab Hierarchical multiype hard core process\cr
    \code{\link{HierStrauss}()}	 \tab Hierarchical multiype Strauss process\cr
    \code{\link{HierStraussHard}()}	 \tab
    Hierarchical multiype Strauss-hard core process\cr
    \code{\link{Hybrid}()}	 \tab Hybrid of several interactions\cr
    \code{\link{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link{MultiHard}()} \tab 		multitype hard core process \cr
    \code{\link{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link{MultiStraussHard}()} \tab 	multitype Strauss/hard core process \cr
    \code{\link{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link{Penttinen}()} \tab	Penttinen pairwise interaction\cr
    \code{\link{SatPiece}()} \tab	Saturated pair model, piecewise  constant potential\cr
    \code{\link{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link{Softcore}()} \tab pairwise interaction, soft core potential\cr
    \code{\link{Strauss}()} \tab  Strauss process \cr
    \code{\link{StraussHard}()} \tab Strauss/hard core point process \cr
    \code{\link{Triplets}()} \tab Geyer triplets process
  }
  Note that it is also possible to combine several such interactions
  using \code{\link{Hybrid}}.
  
  \bold{Finer control over model fitting:}
  
  A quadrature scheme is represented by an object of
  class \code{"quad"}. To create a quadrature scheme, typically
  use \code{\link{quadscheme}}.
  
  \tabular{ll}{
    \code{\link{quadscheme}} \tab default quadrature scheme \cr
                             \tab using rectangular cells or Dirichlet cells\cr
    \code{\link{pixelquad}}  \tab quadrature scheme based on image pixels \cr
    \code{\link{quad}}       \tab create an object of class \code{"quad"}
  }
  
  To inspect a quadrature scheme:
  \tabular{ll}{
    \code{plot(Q)} \tab plot quadrature scheme \code{Q}\cr
    \code{print(Q)} \tab print basic information about quadrature scheme \code{Q}\cr
    \code{\link{summary}(Q)} \tab summary of quadrature scheme \code{Q}
  }

  A quadrature scheme consists of data points, dummy points, and
  weights. To generate dummy points:
  \tabular{ll}{
    \code{\link{default.dummy}} \tab default pattern of dummy points \cr
    \code{\link{gridcentres}} \tab dummy points in a rectangular grid \cr
    \code{\link{rstrat}} \tab stratified random dummy pattern \cr
    \code{\link{spokes}} \tab radial pattern of dummy points  \cr
    \code{\link{corners}} \tab dummy points at corners of the window
  }
  
  To compute weights:
  \tabular{ll}{
    \code{\link{gridweights}} \tab quadrature weights by the grid-counting rule  \cr
    \code{\link{dirichletWeights}} \tab quadrature weights are
    Dirichlet tile areas
  }

  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link{simulate.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }

  \bold{Point process models on a linear network:}

  An object of class \code{"lpp"} represents a pattern of points on
  a linear network. Point process models can also be fitted to these
  objects. Currently only Poisson models can be fitted.

  \tabular{ll}{
    \code{\link{lppm}} \tab point process model on linear network \cr
    \code{\link{anova.lppm}} \tab analysis of deviance for \cr
    \tab point process model on linear network \cr
    \code{\link{envelope.lppm}} \tab simulation envelopes for \cr
    \tab point process model on linear network \cr
    \code{\link{fitted.lppm}} \tab fitted intensity values \cr
    \code{\link{predict.lppm}} \tab model prediction on linear network \cr
    \code{\link{linim}} \tab pixel image on linear network \cr
    \code{\link{plot.linim}} \tab plot a pixel image on linear network \cr
    \code{\link{eval.linim}} \tab evaluate expression involving images \cr
    \code{\link{linfun}} \tab function defined on linear network \cr
    \code{\link{methods.linfun}} \tab conversion facilities
  }
}


\section{V. MODEL FITTING (DETERMINANTAL POINT PROCESS MODELS)}{

  Code for fitting \emph{determinantal point process models} has 
  recently been added to \pkg{spatstat}.

  For information, see the help file for \code{\link{dppm}}.
}

\section{VI. MODEL FITTING (SPATIAL LOGISTIC REGRESSION)}{
  
  \bold{Logistic regression}
  
  Pixel-based spatial logistic regression is an alternative
  technique for analysing spatial point patterns
  that is widely used in Geographical Information Systems.
  It is approximately equivalent to fitting a Poisson point process
  model.
  
  In pixel-based logistic regression, the spatial domain is
  divided into small pixels, the presence or absence of a
  data point in each pixel is recorded, and logistic regression
  is used to model the presence/absence indicators as a function
  of any covariates.
  
  Facilities for performing spatial logistic regression are
  provided in \pkg{spatstat} for comparison purposes.
  
  \bold{Fitting a spatial logistic regression}
  
  Spatial logistic regression is performed by the function
  \code{\link{slrm}}. Its result is an object of class \code{"slrm"}.
  There are many methods for this class, including methods for
  \code{print}, \code{fitted}, \code{predict}, \code{simulate},
  \code{anova}, \code{coef}, \code{logLik}, \code{terms},
  \code{update}, \code{formula} and \code{vcov}. 
  
  For example, if \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{slrm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{slrm(X ~ x)} \tab Poisson process with \cr
                \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{slrm(X ~ Z)} \tab Poisson process with \cr
                \tab intensity loglinear in covariate \code{Z}
  }

  \bold{Manipulating a fitted spatial logistic regression}
  
  \tabular{ll}{
    \code{\link{anova.slrm}} \tab Analysis of deviance \cr
    \code{\link{coef.slrm}}  \tab Extract fitted coefficients \cr
    \code{\link{vcov.slrm}}  \tab Variance-covariance matrix of fitted coefficients \cr
    \code{\link{fitted.slrm}} \tab Compute fitted probabilities or
    intensity \cr
    \code{\link{logLik.slrm}}   \tab Evaluate loglikelihood of fitted
    model \cr
    \code{\link{plot.slrm}}    \tab Plot fitted probabilities or
    intensity \cr
    \code{\link{predict.slrm}} \tab Compute predicted probabilities or
    intensity with new data \cr
    \code{\link{simulate.slrm}} \tab Simulate model
  }
  
  There are many other undocumented methods for this class,
  including methods for \code{print}, \code{update}, \code{formula}
  and \code{terms}. Stepwise model selection is
  possible using \code{step} or \code{stepAIC}.
}


\section{VII. SIMULATION}{

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}

  \tabular{ll}{
    \code{\link{runifpoint}} \tab
    generate \eqn{n} independent uniform random points \cr
    \code{\link{rpoint}} \tab
    generate \eqn{n} independent random points \cr
    \code{\link{rmpoint}} \tab
    generate \eqn{n} independent multitype random points \cr
    \code{\link{rpoispp}} \tab
    simulate the (in)homogeneous Poisson point process \cr
    \code{\link{rmpoispp}} \tab
    simulate the (in)homogeneous multitype Poisson point process \cr
    \code{\link{runifdisc}} \tab
    generate \eqn{n} independent uniform random points in disc\cr
    \code{\link{rstrat}} \tab
    stratified random sample of points \cr
    \code{\link{rsyst}} \tab
    systematic random sample (grid) of points \cr
    \code{\link{rMaternI}}  \tab
    simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Model I inhibition process\cr
    \code{\link{rMaternII}} \tab
    simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Model II inhibition process\cr
    \code{\link{rSSI}} \tab
    simulate Simple Sequential Inhibition process\cr
    \code{\link{rHardcore}} \tab
    simulate hard core process (perfect simulation)\cr
    \code{\link{rStrauss}} \tab
    simulate Strauss process (perfect simulation)\cr
    \code{\link{rStraussHard}} \tab
    simulate Strauss-hard core process (perfect simulation)\cr
    \code{\link{rDiggleGratton}} \tab
    simulate Diggle-Gratton process (perfect simulation)\cr
    \code{\link{rDGS}} \tab
    simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
    \code{\link{rPenttinen}} \tab
    simulate Penttinen process (perfect simulation)\cr
    \code{\link{rNeymanScott}} \tab
    simulate a general Neyman-Scott process\cr
    \code{\link{rMatClust}} \tab
    simulate the \ifelse{latex}{\out{Mat\'ern}}{Matern} Cluster process\cr
    \code{\link{rThomas}} \tab
    simulate the Thomas process  \cr
    \code{\link{rLGCP}} \tab
    simulate the log-Gaussian Cox process  \cr
    \code{\link{rGaussPoisson}}  \tab
    simulate the Gauss-Poisson cluster process\cr
    \code{\link{rCauchy}} \tab
    simulate Neyman-Scott process with Cauchy clusters \cr
    \code{\link{rVarGamma}} \tab
    simulate Neyman-Scott process with Variance Gamma clusters \cr
    \code{\link{rcell}} \tab
    simulate the Baddeley-Silverman cell process  \cr
    \code{\link{runifpointOnLines}} \tab
    generate \eqn{n} random points along specified line segments \cr
    \code{\link{rpoisppOnLines}} \tab
    generate Poisson random points along specified line segments 
  }
      
  \bold{Resampling a point pattern:}

  \tabular{ll}{
    \code{\link{quadratresample}} \tab block resampling \cr
    \code{\link{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link{rthin}} \tab  random thinning 
  }
  
  See also \code{\link{varblock}} for estimating the variance
  of a summary statistic by block resampling, and
  \code{\link{lohboot}} for another bootstrap technique.
  
  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.

   Cluster process models 
   are fitted by the function \code{\link{kppm}} yielding an
   object of class \code{"kppm"}. To generate one or more simulated
   realisations of this fitted model, use 
   \code{\link{simulate.kppm}}.

   Gibbs point process models 
   are fitted by the function \code{\link{ppm}} yielding an
   object of class \code{"ppm"}. To generate a simulated
   realisation of this fitted model, use \code{\link{rmh}}.
   To generate one or more simulated realisations of the fitted model,
   use \code{\link{simulate.ppm}}.

   \bold{Other random patterns:}

   \tabular{ll}{
     \code{\link{rlinegrid}} \tab
     generate a random array of parallel lines through a window \cr
     \code{\link{rpoisline}} \tab
     simulate the Poisson line process within a window \cr
     \code{\link{rpoislinetess}} \tab
     generate random tessellation using Poisson line process \cr
     \code{\link{rMosaicSet}} \tab
     generate random set by selecting some tiles of a tessellation \cr
     \code{\link{rMosaicField}} \tab
     generate random pixel image by assigning random values
     in each tile of a tessellation
   }

   \bold{Simulation-based inference}

   \tabular{ll}{
    \code{\link{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link{qqplot.ppm}} \tab diagnostic plot for interpoint
    interaction \cr
    \code{\link{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link{studpermu.test}} \tab studentised permutation test\cr
    \code{\link{segregation.test}} \tab test of segregation of types 
  }
}


\section{VIII. TESTS AND DIAGNOSTICS}{

  \bold{Classical hypothesis tests:}
  
  \tabular{ll}{
    \code{\link{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link{clarkevans.test}} \tab Clark and Evans test \cr
    \code{\link{cdf.test}} \tab Spatial distribution goodness-of-fit test\cr
    \code{\link{berman.test}} \tab Berman's goodness-of-fit tests\cr
    \code{\link{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link{dclf.test}} \tab Diggle-Cressie-Loosmore-Ford test \cr
    \code{\link{mad.test}} \tab Mean Absolute Deviation test \cr
    \code{\link{dclf.progress}} \tab Progress plot for DCLF test \cr
    \code{\link{mad.progress}} \tab Progress plot for MAD test \cr
    \code{\link{anova.ppm}} \tab Analysis of Deviance for
    point process models 
  }

  \bold{Sensitivity diagnostics:}

  Classical measures of model sensitivity such as leverage and influence
  have been adapted to point process models.
  
    \tabular{ll}{
      \code{\link{leverage.ppm}} \tab Leverage for point process model\cr
      \code{\link{influence.ppm}} \tab Influence for point process model\cr
      \code{\link{dfbetas.ppm}} \tab Parameter influence\cr
    }
  
  \bold{Diagnostics for covariate effect:}

  Classical diagnostics for covariate effects have been adapted to
  point process models.

  \tabular{ll}{
    \code{\link{parres}} \tab Partial residual plot\cr
    \code{\link{addvar}} \tab Added variable plot \cr
    \code{\link{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link{rho2hat}} \tab Kernel estimate of covariate effect
    (bivariate) 
  }
  
  \bold{Residual diagnostics:}
  
  Residuals for a fitted point process model, and diagnostic plots
  based on the residuals, were introduced in Baddeley et al (2005) and
  Baddeley, Rubak and \ifelse{latex}{\out{M\o ller}}{Moller} (2011).
  
  Type \code{demo(diagnose)}
  for a demonstration of the diagnostics features.

  \tabular{ll}{
    \code{\link{diagnose.ppm}} \tab diagnostic plots for spatial trend\cr
    \code{\link{qqplot.ppm}} \tab diagnostic Q-Q plot for interpoint interaction\cr
    \code{\link{residualspaper}} \tab examples from Baddeley et al (2005) \cr
    \code{\link{Kcom}} \tab model compensator of \eqn{K} function \cr
    \code{\link{Gcom}} \tab model compensator of \eqn{G} function \cr
    \code{\link{Kres}} \tab score residual of \eqn{K} function \cr
    \code{\link{Gres}} \tab score residual of \eqn{G} function \cr
    \code{\link{psst}} \tab pseudoscore residual of summary function \cr
    \code{\link{psstA}} \tab pseudoscore residual of empty space function \cr
    \code{\link{psstG}} \tab pseudoscore residual of \eqn{G} function \cr
    \code{\link{compareFit}} \tab compare compensators of several fitted
    models
  }


  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link{quadratresample}} \tab block resampling \cr
    \code{\link{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link{rthin}} \tab  random thinning  
  }
}


\section{IX. DOCUMENTATION}{
  The online manual entries are quite detailed and should be consulted
  first for information about a particular function.
  
  The book
  Baddeley, Rubak and Turner (2015) is a complete course
  on analysing spatial point patterns, with full details about
  \pkg{spatstat}.

  Older material (which is now out-of-date but is freely available)
  includes Baddeley and Turner (2005a), a brief overview of
  the package in its early development;
  Baddeley and Turner (2005b), a more detailed explanation of
  how to fit point process models to data; and 
  Baddeley (2010), a complete set of notes from a 2-day workshop
  on the use of \pkg{spatstat}. 

  Type \code{citation("spatstat")} to get a list of these references.
}
\references{
  Baddeley, A. (2010)
  \emph{Analysing spatial point patterns in R}.
  Workshop notes. Version 4.1. CSIRO online technical publication.
  \url{https://research.csiro.au/software/r-workshop-notes/}

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
  
  Baddeley, A. and Turner, R. (2005a)
  Spatstat: an R package for analyzing spatial point patterns.
  \emph{Journal of Statistical Software} \bold{12}:6, 1--42.
  URL: \code{www.jstatsoft.org}, ISSN: 1548-7660.

  Baddeley, A. and Turner, R. (2005b)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.

  Baddeley, A., Turner, R.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., Rubak, E. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.

  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  \emph{Journal of Statistical Software} \bold{55}:11, 1--43.
  \url{http://www.jstatsoft.org/v55/i11/}

  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.

  Diggle, P.J. (2014)
  \emph{Statistical Analysis of Spatial and Spatio-Temporal Point Patterns},
  Third edition. {Chapman and Hall/CRC}.

  Gelfand, A.E., Diggle, P.J., Fuentes, M. and Guttorp, P., editors (2010)
  \emph{Handbook of Spatial Statistics}.
  CRC Press.

  Huang, F. and Ogata, Y. (1999)
  Improvements of the maximum pseudo-likelihood
  estimators in various spatial statistical models.
  \emph{Journal of Computational and Graphical Statistics}
  \bold{8}, 510--530.

  Illian, J., Penttinen, A., Stoyan, H. and Stoyan, D. (2008)
  \emph{Statistical Analysis and Modelling of Spatial Point Patterns.}
  Wiley.
  
  Waagepetersen, R.
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63} (2007) 252--258.
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie, 
  Ute Hahn,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Tuomas Rajala,
  Dominic Schuhmacher and
  Rasmus Waagepetersen
  made substantial contributions of code.

  Additional contributions by
  Monsuru Adepeju,
  Ang Qi Wei,
  Sandro Azaele,
  Malissa Baddeley,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  Jean-Francois Coeurjolly,
  Robin Corria Ainslie,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Julian Gilbey,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \ifelse{latex}{\out{B\o gsted}}{Bogsted} Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Paul Hewson,
  Kassel Hingee,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Frederic Lavancier,
  Tom Lawrence,
  Robert Lamb,
  Jonathan Lee,
  George Leser,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \ifelse{latex}{\out{M\o ller}}{Moller},
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \ifelse{latex}{\out{Oehlschl\"{a}gel}}{Oehlschlaegel},
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Suman Rakshit,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \ifelse{latex}{\out{S\"{a}rkk\"{a}}}{Sarkka},
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

