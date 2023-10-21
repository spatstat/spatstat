\name{spatstat-package}
\alias{spatstat-package} 
\alias{spatstat} 
\docType{package}
\title{The Spatstat Package}
\description{
  This is a summary of the features of 
  \pkg{spatstat}, a family of \R packages
  for the statistical analysis of spatial point patterns.
}
\details{
  \pkg{spatstat} is a family of \R packages
  for the statistical analysis
  of spatial data. Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  \pkg{spatstat} is designed to support a complete statistical analysis
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
  a function analogous to \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
  Fitted models can be printed, plotted, predicted, simulated and so on.
}
\section{Getting Started}{
  For a quick introduction to \pkg{spatstat}, read
  the package vignette \emph{Getting started with spatstat}
  installed with \pkg{spatstat}. To read that document, you can either
  \itemize{
    \item visit \url{https://cran.r-project.org/package=spatstat}
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
  at the book companion website, \url{https://book.spatstat.org/}.

%  For information about handling data in \bold{shapefiles},
%  see Chapter 3, or the Vignette
%  \emph{Handling shapefiles in the spatstat package},
%  installed with \pkg{spatstat}, accessible as
%  \code{vignette('shapefiles')}.
}
\section{Structure of the spatstat family}{  
  The original \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing functionality for geometrical
    operations, and defining the main classes of spatial objects
    \item \pkg{spatstat.explore} containing the main functions for exploratory
    analysis of spatial data
    \item \pkg{spatstat.model} containing the main functions for
    parametric statistical modelling and analysis,
    and formal inference, for spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  The breakup has been done in such a way that the user should not notice
  any difference. Source code that worked with the old \pkg{spatstat}
  package should work with the new \pkg{spatstat} family.
  Code that is documented in our books, journal articles and vignettes
  should still work.
  
  When you install \pkg{spatstat}, the sub-packages listed above are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported. 

  This help file covers all the functionality and datasets
  that are provided in the sub-packages listed above.
}
\section{Extension packages}{
  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
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
    \code{\link[spatstat.geom]{ppp}} \tab point pattern \cr
    \code{\link[spatstat.geom]{owin}} \tab window (spatial region) \cr
    \code{\link[spatstat.geom]{im}} \tab pixel image \cr
    \code{\link[spatstat.geom]{psp}} \tab line segment pattern \cr
    \code{\link[spatstat.geom]{tess}} \tab tessellation \cr
    \code{\link[spatstat.geom]{pp3}} \tab three-dimensional point pattern \cr
    \code{\link[spatstat.geom]{ppx}} \tab point pattern in any number of dimensions \cr
    \code{\link[spatstat.linnet]{lpp}} \tab point pattern on a linear network
  }

  \bold{To create a point pattern:}
  
  \tabular{ll}{
    \code{\link[spatstat.geom]{ppp}} \tab
    create a point pattern from \eqn{(x,y)} and window information
    \cr
    \tab
    \code{ppp(x, y, xlim, ylim)} for rectangular window\cr
    \tab
    \code{ppp(x, y, poly)} for polygonal window \cr
    \tab
    \code{ppp(x, y, mask)} for binary image window \cr
    \code{\link[spatstat.geom]{as.ppp}} \tab
    convert other types of data to a \code{ppp} object \cr
    \code{\link[spatstat.geom]{clickppp}} \tab
    interactively add points to a plot \cr
    \code{\link[spatstat.geom]{marks<-}}, \code{\%mark\%}  \tab
    attach/reassign marks to a point pattern
  }
      
  \bold{To simulate a random point pattern:}
  
  \tabular{ll}{
	\code{\link[spatstat.random]{runifpoint}} \tab
	generate \eqn{n} independent uniform random points \cr
	\code{\link[spatstat.random]{rpoint}} \tab
	generate \eqn{n} independent random points \cr
	\code{\link[spatstat.random]{rmpoint}} \tab
	generate \eqn{n} independent multitype random points \cr
	\code{\link[spatstat.random]{rpoispp}} \tab
	simulate the (in)homogeneous Poisson point process \cr
	\code{\link[spatstat.random]{rmpoispp}} \tab
	simulate the (in)homogeneous multitype Poisson point process \cr
	\code{\link[spatstat.random]{runifdisc}} \tab
	generate \eqn{n} independent uniform random points in disc\cr
	\code{\link[spatstat.random]{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link[spatstat.geom]{rsyst}} \tab
	systematic random sample of points \cr
	\code{\link[spatstat.geom]{rjitter}} \tab
	apply random displacements to points in a pattern\cr
	\code{\link[spatstat.random]{rMaternI}}  \tab
	simulate the \Matern Model I inhibition process\cr
	\code{\link[spatstat.random]{rMaternII}} \tab
	simulate the \Matern Model II inhibition process\cr
	\code{\link[spatstat.random]{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link[spatstat.random]{rStrauss}} \tab
	simulate Strauss process (perfect simulation)\cr
	\code{\link[spatstat.random]{rHardcore}} \tab
	simulate Hard Core process (perfect simulation)\cr
	\code{\link[spatstat.random]{rStraussHard}} \tab
	simulate Strauss-hard core process (perfect simulation)\cr
	\code{\link[spatstat.random]{rDiggleGratton}} \tab
	simulate Diggle-Gratton process (perfect simulation)\cr
	\code{\link[spatstat.random]{rDGS}} \tab
	simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
	\code{\link[spatstat.random]{rPenttinen}} \tab
	simulate Penttinen process (perfect simulation)\cr
	\code{\link[spatstat.random]{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link[spatstat.random]{rPoissonCluster}} \tab
	simulate a general Poisson cluster process\cr
	\code{\link[spatstat.random]{rMatClust}} \tab
	simulate the \Matern Cluster process\cr
	\code{\link[spatstat.random]{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link[spatstat.random]{rGaussPoisson}}  \tab
	simulate the Gauss-Poisson cluster process\cr
	\code{\link[spatstat.random]{rCauchy}} \tab
	simulate Neyman-Scott Cauchy cluster process \cr
	\code{\link[spatstat.random]{rVarGamma}} \tab
	simulate Neyman-Scott Variance Gamma cluster process \cr
	\code{\link[spatstat.random]{rthin}} \tab  random thinning  \cr
	\code{\link[spatstat.random]{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link[spatstat.random]{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link[spatstat.model]{simulate.ppm}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link[spatstat.random]{runifpointOnLines}} \tab
	generate \eqn{n} random points along specified line segments \cr
	\code{\link[spatstat.random]{rpoisppOnLines}} \tab
	generate Poisson random points along specified line segments 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link[spatstat.random]{rshift}} \tab random shifting of points \cr
        \code{\link[spatstat.geom]{rjitter}} \tab
	apply random displacements to points in a pattern\cr
	\code{\link[spatstat.random]{rthin}} \tab  random thinning \cr
	\code{\link[spatstat.random]{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
        \code{\link[spatstat.random]{quadratresample}} \tab block resampling 
      }

      \bold{Standard point pattern datasets:}

      Datasets in \pkg{spatstat} are lazy-loaded, so you can simply
      type the name of the dataset to use it; there is no need
      to type \code{\link[utils]{data}(amacrine)} etc.

      Type \code{demo(data)} to see a display of all the datasets
      installed with the package.

      Type \code{vignette('datasets')} for a document giving an overview
      of all datasets, including background information, and plots.
      
      \tabular{ll}{
	\code{\link[spatstat.data]{amacrine}} \tab Austin Hughes' rabbit amacrine cells \cr
	\code{\link[spatstat.data]{anemones}} \tab Upton-Fingleton sea anemones data\cr
	\code{\link[spatstat.data]{ants}} \tab Harkness-Isham ant nests data\cr
	\code{\link[spatstat.data]{bdspots}} \tab Breakdown spots in microelectrodes \cr
	\code{\link[spatstat.data]{bei}} \tab Tropical rainforest trees \cr
	\code{\link[spatstat.data]{betacells}} \tab Waessle et al. cat retinal ganglia data \cr
	\code{\link[spatstat.data]{bramblecanes}} \tab Bramble Canes data \cr
	\code{\link[spatstat.data]{bronzefilter}} \tab Bronze Filter Section data \cr
	\code{\link[spatstat.data]{cells}} \tab Crick-Ripley biological cells data \cr
	\code{\link[spatstat.data]{chicago}} \tab Chicago crimes \cr
	\code{\link[spatstat.data]{chorley}} \tab Chorley-Ribble cancer data \cr
	\code{\link[spatstat.data]{clmfires}} \tab Castilla-La Mancha forest fires \cr
	\code{\link[spatstat.data]{copper}} \tab Berman-Huntington copper deposits data \cr
	\code{\link[spatstat.data]{dendrite}} \tab Dendritic spines \cr
	\code{\link[spatstat.data]{demohyper}} \tab Synthetic point patterns\cr
	\code{\link[spatstat.data]{demopat}} \tab Synthetic point pattern \cr
	\code{\link[spatstat.data]{finpines}} \tab Finnish Pines data \cr
	\code{\link[spatstat.data]{flu}} \tab Influenza virus proteins \cr
	\code{\link[spatstat.data]{gordon}} \tab People in Gordon Square, London \cr
	\code{\link[spatstat.data]{gorillas}} \tab Gorilla nest sites \cr
	\code{\link[spatstat.data]{hamster}} \tab Aherne's hamster tumour data \cr
	\code{\link[spatstat.data]{humberside}} \tab North Humberside childhood leukaemia data \cr
	\code{\link[spatstat.data]{hyytiala}} \tab {Mixed forest in
              \ifelse{latex}{\out{Hyyti{\"a}l{\"a}}}{Hyytiala}, Finland}\cr
	\code{\link[spatstat.data]{japanesepines}} \tab Japanese Pines data \cr
	\code{\link[spatstat.data]{lansing}} \tab Lansing Woods data \cr
	\code{\link[spatstat.data]{longleaf}} \tab Longleaf Pines data \cr
	\code{\link[spatstat.data]{mucosa}} \tab Cells in gastric mucosa \cr
	\code{\link[spatstat.data]{murchison}} \tab Murchison gold deposits \cr
	\code{\link[spatstat.data]{nbfires}} \tab New Brunswick fires data \cr
	\code{\link[spatstat.data]{nztrees}} \tab Mark-Esler-Ripley trees data \cr
	\code{\link[spatstat.data]{osteo}} \tab Osteocyte lacunae (3D, replicated) \cr
	\code{\link[spatstat.data]{paracou}} \tab Kimboto trees in Paracou, French Guiana \cr
	\code{\link[spatstat.data]{ponderosa}} \tab Getis-Franklin ponderosa pine trees data \cr
	\code{\link[spatstat.data]{pyramidal}} \tab Pyramidal neurons from 31 brains \cr
	\code{\link[spatstat.data]{redwood}} \tab Strauss-Ripley redwood saplings data \cr
	\code{\link[spatstat.data]{redwoodfull}} \tab Strauss redwood saplings data (full set) \cr
	\code{\link[spatstat.data]{residualspaper}} \tab Data from Baddeley et al (2005) \cr
	\code{\link[spatstat.data]{shapley}} \tab Galaxies in an astronomical survey \cr
	\code{\link[spatstat.data]{simdat}} \tab Simulated point pattern (inhomogeneous, with interaction) \cr
	\code{\link[spatstat.data]{spiders}} \tab Spider webs on mortar lines of brick wall \cr
	\code{\link[spatstat.data]{sporophores}} \tab Mycorrhizal fungi around a tree \cr
	\code{\link[spatstat.data]{spruces}} \tab Spruce trees in Saxonia \cr
	\code{\link[spatstat.data]{swedishpines}} \tab Strand-Ripley Swedish pines data \cr
	\code{\link[spatstat.data]{urkiola}} \tab Urkiola Woods data \cr
	\code{\link[spatstat.data]{waka}} \tab Trees in Waka national park \cr
	\code{\link[spatstat.data]{waterstriders}} \tab Insects on water surface 
      }

      \bold{To manipulate a point pattern:}

      \tabular{ll}{
	\code{\link[spatstat.geom]{plot.ppp}} \tab
	plot a point pattern (e.g. \code{plot(X)}) \cr
	\code{spatstat.gui::iplot} \tab
	plot a point pattern interactively \cr
	\code{\link[spatstat.geom]{edit.ppp}} \tab interactive text editor \cr
	\code{\link[spatstat.geom]{[.ppp}} \tab
	extract or replace a subset of a point pattern \cr
        \tab \code{pp[subset]} or \code{pp[subwindow]} \cr
	\code{\link[spatstat.geom]{subset.ppp}} \tab
	extract subset of point pattern satisfying a condition \cr
	\code{\link[spatstat.geom]{superimpose}} \tab
	combine several point patterns  \cr
	\code{\link[spatstat.geom]{by.ppp}} \tab
	apply a function to sub-patterns of a point pattern \cr
	\code{\link[spatstat.geom]{cut.ppp}} \tab
	classify the points in a point pattern \cr
	\code{\link[spatstat.geom]{split.ppp}} \tab
	divide pattern into sub-patterns \cr
	\code{\link[spatstat.geom]{unmark}} \tab
	remove marks  \cr
	\code{\link[spatstat.geom]{npoints}} \tab
	count the number of points  \cr
	\code{\link[spatstat.geom]{coords}} \tab
	extract coordinates, change coordinates  \cr
	\code{\link[spatstat.geom]{marks}} \tab
	extract marks, change marks or attach marks  \cr
	\code{\link[spatstat.geom]{rotate}} \tab
	rotate pattern  \cr
	\code{\link[spatstat.geom]{shift} } \tab
	translate pattern  \cr
	\code{\link[spatstat.geom]{flipxy} } \tab
	swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link[spatstat.geom]{reflect} } \tab
	reflect in the origin  \cr
	\code{\link[spatstat.geom]{periodify} } \tab
	make several translated copies  \cr
	\code{\link[spatstat.geom]{affine}} \tab
	apply affine transformation\cr
	\code{\link[spatstat.geom]{scalardilate}} \tab
	apply scalar dilation\cr
	\code{\link[spatstat.explore]{density.ppp}} \tab
	kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.explore]{densityHeat.ppp}} \tab
	diffusion kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.explore]{Smooth.ppp}} \tab
	kernel smoothing of marks of point pattern\cr
	\code{\link[spatstat.geom]{nnmark}} \tab
	mark value of nearest data point\cr
	\code{\link[spatstat.explore]{sharpen.ppp}} \tab
	data sharpening\cr
	\code{\link[spatstat.geom]{identify.ppp}} \tab
	interactively identify points \cr
	\code{\link[spatstat.geom]{unique.ppp}} \tab
	remove duplicate points \cr
	\code{\link[spatstat.geom]{duplicated.ppp}} \tab
	determine which points are duplicates \cr
	\code{\link[spatstat.geom]{uniquemap.ppp}} \tab
	map duplicated points to unique points \cr
	\code{\link[spatstat.geom]{connected.ppp}} \tab find clumps of points \cr
	\code{\link[spatstat.geom]{dirichlet}} \tab
	compute Dirichlet-Voronoi tessellation \cr
	\code{\link[spatstat.geom]{delaunay}} \tab
	compute Delaunay triangulation \cr 
	\code{\link[spatstat.geom]{delaunayDistance}} \tab
	graph distance in Delaunay triangulation \cr 
	\code{\link[spatstat.geom]{convexhull}} \tab compute convex hull \cr
	\code{\link[spatstat.geom]{discretise}} \tab discretise coordinates \cr
	\code{\link[spatstat.geom]{pixellate.ppp}} \tab approximate point pattern by 
	pixel image \cr
	\code{\link[spatstat.geom]{as.im.ppp}} \tab approximate point pattern by 
	pixel image 
      }
      See \code{\link[spatstat.geom]{spatstat.options}} to control plotting behaviour.
      
    \bold{To create a window:}

    An object of class \code{"owin"} describes a spatial region
      (a window of observation).

      \tabular{ll}{
	\code{\link[spatstat.geom]{owin}}	\tab	Create a window object \cr
	\tab \code{owin(xlim, ylim)} for rectangular window \cr
	\tab \code{owin(poly)} for polygonal window \cr
	\tab \code{owin(mask)} for binary image window \cr
	\code{\link[spatstat.geom]{Window}}	\tab
	Extract window of another object \cr
	\code{\link[spatstat.geom]{Frame}}	\tab
	Extract the containing rectangle ('frame') of another object \cr
	\code{\link[spatstat.geom]{as.owin}}	\tab
	Convert other data to a window object \cr
	\code{\link[spatstat.geom]{square}}    \tab	make a square window \cr
	\code{\link[spatstat.geom]{disc}}    \tab	make a circular window \cr
	\code{\link[spatstat.geom]{ellipse}}    \tab	make an elliptical window \cr
	\code{\link[spatstat.geom]{ripras}}    \tab
	Ripley-Rasson estimator of window, given only the points \cr
	\code{\link[spatstat.geom]{convexhull}} \tab compute convex hull of something \cr 
	\code{\link[spatstat.data]{letterR}}    \tab
	polygonal window in the shape of the \R logo \cr
	\code{\link[spatstat.geom]{clickpoly}}    \tab
	interactively draw a polygonal window  \cr
	\code{\link[spatstat.geom]{clickbox}}    \tab
	interactively draw a rectangle  
      }

    \bold{To manipulate a window:}

    \tabular{ll}{
	\code{\link[spatstat.geom]{plot.owin}}	\tab	plot a window. \cr
	\tab		\code{plot(W)}\cr
	\code{\link[spatstat.geom]{boundingbox}} \tab
	Find a tight bounding box for the window \cr
	\code{\link[spatstat.geom]{erosion}}	\tab
	erode window by a distance r\cr
	\code{\link[spatstat.geom]{dilation}}	\tab
	dilate window by a distance r\cr
	\code{\link[spatstat.geom]{closing}}	\tab
	close window by a distance r\cr
	\code{\link[spatstat.geom]{opening}}	\tab
	open window by a distance r\cr
	\code{\link[spatstat.geom]{border}}	\tab
	difference between window and its erosion/dilation \cr 
	\code{\link[spatstat.geom]{complement.owin}}	\tab
	invert (swap inside and outside)\cr
	\code{\link[spatstat.geom]{simplify.owin}}	\tab
	approximate a window by a simple polygon  \cr
	\code{\link[spatstat.geom]{rotate}} \tab rotate window  \cr
	\code{\link[spatstat.geom]{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link[spatstat.geom]{shift} } \tab translate window  \cr
	\code{\link[spatstat.geom]{periodify} } \tab make several translated copies  \cr
	\code{\link[spatstat.geom]{affine}} \tab apply affine transformation \cr
	\code{\link[spatstat.geom]{as.data.frame.owin}} \tab
	convert window to data frame 
      }

    \bold{Digital approximations:}

    \tabular{ll}{
	\code{\link[spatstat.geom]{as.mask}}	\tab
	Make a discrete pixel approximation of a given window \cr
	\code{\link[spatstat.geom]{as.im.owin}} \tab convert window to pixel image \cr
	\code{\link[spatstat.geom]{pixellate.owin}} \tab convert window to pixel image \cr
	\code{\link[spatstat.geom]{commonGrid}} \tab find common pixel grid for windows \cr
	\code{\link[spatstat.geom]{nearest.raster.point}} \tab
	map continuous coordinates to raster locations\cr
	\code{\link[spatstat.geom]{raster.x}} \tab
	raster x coordinates \cr
	\code{\link[spatstat.geom]{raster.y}} \tab
	raster y coordinates \cr
	\code{\link[spatstat.geom]{raster.xy}} \tab
	raster x and y coordinates \cr
	\code{\link[spatstat.geom]{as.polygonal}} \tab
	convert pixel mask to polygonal window 
      }
      See \code{\link[spatstat.geom]{spatstat.options}} to control the approximation

    \bold{Geometrical computations with windows:}

    \tabular{ll}{
	\code{\link[spatstat.geom]{edges}}	\tab	extract boundary edges \cr
	\code{\link[spatstat.geom]{intersect.owin}}	\tab	intersection of two windows\cr
	\code{\link[spatstat.geom]{union.owin}}	\tab	union of two windows\cr
	\code{\link[spatstat.geom]{setminus.owin}}	\tab	set subtraction of two windows\cr
	\code{\link[spatstat.geom]{inside.owin}}	\tab	determine whether a point is inside a window\cr
	\code{\link[spatstat.geom]{area.owin}}	\tab	compute area \cr
	\code{\link[spatstat.geom]{perimeter}}	\tab	compute perimeter length \cr
	\code{\link[spatstat.geom]{diameter.owin}}	\tab	compute diameter\cr
	\code{\link[spatstat.geom]{incircle}}	\tab	find largest circle inside a window \cr
	\code{\link[spatstat.geom]{inradius}}	\tab	radius of incircle \cr
	\code{\link[spatstat.geom]{connected.owin}}    \tab find connected components of window \cr
	\code{\link[spatstat.geom]{eroded.areas}}	\tab	compute areas of eroded windows\cr
	\code{\link[spatstat.geom]{dilated.areas}}	\tab	compute areas of dilated windows\cr
	\code{\link[spatstat.geom]{bdist.points}}	\tab	compute distances from data points to window boundary \cr
	\code{\link[spatstat.geom]{bdist.pixels}}	\tab	compute distances from all pixels to window boundary \cr
	\code{\link[spatstat.geom]{bdist.tiles}}	\tab
	boundary distance for each tile in tessellation \cr
	\code{\link[spatstat.geom]{distmap.owin}}	\tab	distance transform image \cr
	\code{\link[spatstat.geom]{distfun.owin}}	\tab	distance transform \cr
	\code{\link[spatstat.geom]{centroid.owin}}	\tab compute centroid (centre of mass) of window\cr
	\code{\link[spatstat.geom]{is.subset.owin}}    \tab determine whether one
	window contains another \cr
	\code{\link[spatstat.geom]{is.convex}} \tab determine whether a window is convex \cr 
	\code{\link[spatstat.geom]{convexhull}} \tab compute convex hull \cr
	\code{\link[spatstat.geom]{triangulate.owin}} \tab decompose into triangles \cr
	\code{\link[spatstat.geom]{as.mask}} \tab pixel approximation of window \cr
	\code{\link[spatstat.geom]{as.polygonal}} \tab polygonal approximation of window \cr
	\code{\link[spatstat.geom]{is.rectangle}} \tab test whether window is a rectangle \cr 
	\code{\link[spatstat.geom]{is.polygonal}} \tab test whether window is polygonal \cr 
	\code{\link[spatstat.geom]{is.mask}} \tab test whether window is a mask \cr 
	\code{\link[spatstat.geom]{setcov}} \tab spatial covariance function of window \cr
	\code{\link[spatstat.geom]{pixelcentres}} \tab extract centres of pixels in mask \cr
	\code{\link[spatstat.geom]{clickdist}}    \tab
	measure distance between two points clicked by user
      }

    \bold{Pixel images:}
    An object of class \code{"im"} represents a pixel image. 
    Such objects are returned by some of the functions in
      \pkg{spatstat} including \code{\link[spatstat.explore]{Kmeasure}},
      \code{\link[spatstat.geom]{setcov}} and \code{\link[spatstat.explore]{density.ppp}}. 
      \tabular{ll}{
	\code{\link[spatstat.geom]{im}} \tab create a pixel image\cr
	\code{\link[spatstat.geom]{as.im}} \tab convert other data to a pixel image\cr
	\code{\link[spatstat.geom]{pixellate}} \tab convert other data to a pixel image\cr
	\code{\link[spatstat.geom]{as.matrix.im}} \tab convert pixel image to matrix\cr
	\code{\link[spatstat.geom]{as.data.frame.im}} \tab convert pixel image to data frame\cr
	\code{\link[spatstat.geom]{as.function.im}} \tab convert pixel image to function\cr
	\code{\link[spatstat.geom]{plot.im}}	\tab	plot a pixel image on screen as a digital image\cr
	\code{\link[spatstat.geom]{contour.im}}	\tab draw contours of a pixel image \cr
	\code{\link[spatstat.geom]{persp.im}}	\tab draw perspective plot of a pixel image \cr
	\code{\link[spatstat.geom]{rgbim}}	\tab create colour-valued pixel image \cr
	\code{\link[spatstat.geom]{hsvim}}	\tab create colour-valued pixel image \cr
	\code{\link[spatstat.geom]{[.im}} 	\tab extract a subset of a pixel image\cr
	\code{\link[spatstat.geom]{[<-.im}} 	\tab replace a subset of a pixel image\cr
	\code{\link[spatstat.geom]{rotate.im}} \tab rotate pixel image \cr
	\code{\link[spatstat.geom]{shift.im}} \tab apply vector shift to pixel image \cr
	\code{\link[spatstat.geom]{affine.im}} \tab apply affine transformation to image \cr
	\code{X}	\tab	print very basic information about image \code{X}\cr
	\code{\link[spatstat.geom:summary.im]{summary}(X)} \tab	summary of image \code{X} \cr
	\code{\link[spatstat.geom]{hist.im}} \tab	histogram of image \cr
	\code{\link[spatstat.geom]{mean.im}} \tab	mean pixel value of image  \cr
	\code{\link[spatstat.geom]{integral.im}} \tab	integral of pixel values  \cr
	\code{\link[spatstat.geom]{quantile.im}} \tab	quantiles of image  \cr
	\code{\link[spatstat.geom]{cut.im}} \tab	convert numeric image to factor image \cr
	\code{\link[spatstat.geom]{is.im}} \tab test whether an object is a pixel image\cr
	\code{\link[spatstat.geom]{interp.im}} \tab interpolate a pixel image\cr
	\code{\link[spatstat.explore]{blur}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.explore]{Smooth.im}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.geom]{connected.im}} \tab find connected components \cr
	\code{\link[spatstat.geom]{compatible.im}} \tab test whether two images have
	compatible dimensions \cr
	\code{\link[spatstat.geom]{harmonise.im}} \tab make images compatible \cr
	\code{\link[spatstat.geom]{commonGrid}} \tab find a common pixel grid for images \cr
	\code{\link[spatstat.geom]{eval.im}} \tab evaluate any expression involving images\cr
	\code{\link[spatstat.geom]{im.apply}} \tab evaluate a function of several images \cr
	\code{\link[spatstat.geom]{scaletointerval}} \tab rescale pixel values \cr
	\code{\link[spatstat.geom]{zapsmall.im}} \tab set very small pixel values to zero \cr
	\code{\link[spatstat.geom]{levelset}} \tab level set of an image\cr 
	\code{\link[spatstat.geom]{solutionset}} \tab region where an expression is true \cr
	\code{\link[spatstat.geom]{imcov}} \tab spatial covariance function of image \cr
	\code{\link[spatstat.geom]{convolve.im}} \tab spatial convolution of images \cr
	\code{\link[spatstat.explore]{transect.im}} \tab line transect of image \cr
	\code{\link[spatstat.geom]{pixelcentres}} \tab extract centres of pixels \cr
	\code{\link[spatstat.geom]{transmat}} \tab convert matrix of pixel values \cr
	                       \tab to a different indexing convention \cr
	\code{\link[spatstat.random]{rnoise}} \tab random pixel noise
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.
    \tabular{ll}{
      \code{\link[spatstat.geom]{psp}} \tab create a line segment pattern \cr
      \code{\link[spatstat.geom]{as.psp}} \tab convert other data into a line segment pattern \cr
      \code{\link[spatstat.geom]{edges}} \tab extract edges of a window \cr
      \code{\link[spatstat.geom]{is.psp}} \tab determine whether a dataset has class \code{"psp"} \cr
      \code{\link[spatstat.geom]{plot.psp}} \tab plot a line segment pattern \cr
      \code{\link[spatstat.geom]{print.psp}} \tab print basic information \cr
      \code{\link[spatstat.geom]{summary.psp}} \tab print summary information \cr
      \code{\link[spatstat.geom]{[.psp}} \tab extract a subset of a line segment pattern \cr
      \code{\link[spatstat.geom]{subset.psp}} \tab extract subset of line segment pattern \cr
      \code{\link[spatstat.geom]{as.data.frame.psp}} \tab
      convert line segment pattern to data frame \cr
      \code{\link[spatstat.geom]{marks.psp}} \tab extract marks of line segments \cr
      \code{\link[spatstat.geom]{marks<-.psp}} \tab assign new marks to line segments \cr
      \code{\link[spatstat.geom]{unmark.psp}} \tab delete marks from line segments \cr
      \code{\link[spatstat.geom]{midpoints.psp}} \tab
      compute the midpoints of line segments \cr
      \code{\link[spatstat.geom]{endpoints.psp}} \tab extract the endpoints of line segments \cr
      \code{\link[spatstat.geom]{lengths_psp}} \tab compute the lengths of line segments \cr
      \code{\link[spatstat.geom]{angles.psp}} \tab compute the orientation angles of line segments \cr
      \code{\link[spatstat.geom]{superimpose}} \tab combine several line segment patterns  \cr
      \code{\link[spatstat.geom]{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates \cr      
      \code{\link[spatstat.geom]{rotate.psp}} \tab rotate a line segment pattern \cr      
      \code{\link[spatstat.geom]{shift.psp}} \tab shift a line segment pattern \cr      
      \code{\link[spatstat.geom]{periodify}} \tab make several shifted copies \cr      
      \code{\link[spatstat.geom]{affine.psp}} \tab apply an affine transformation \cr      
      \code{\link[spatstat.geom]{pixellate.psp}} \tab approximate line segment pattern
    by pixel image \cr      
      \code{\link[spatstat.geom]{psp2mask}} \tab approximate line segment pattern
    by binary mask \cr      
      \code{\link[spatstat.geom]{distmap.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link[spatstat.geom]{distfun.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link[spatstat.explore]{density.psp}} \tab kernel smoothing of line segments\cr
      \code{\link[spatstat.geom]{selfcrossing.psp}} \tab find crossing points between
      line segments \cr
      \code{\link[spatstat.geom]{selfcut.psp}} \tab cut segments where they cross \cr
      \code{\link[spatstat.geom]{crossing.psp}} \tab find crossing points between
      two line segment patterns \cr
      \code{\link[spatstat.geom]{extrapolate.psp}} \tab extrapolate line segments to
      infinite lines \cr
      \code{\link[spatstat.geom]{nncross}} \tab find distance to nearest line segment
      from a given point\cr
      \code{\link[spatstat.geom]{nearestsegment}} \tab find line segment closest to a
      given point \cr
      \code{\link[spatstat.geom]{project2segment}} \tab find location along a line segment
      closest to a given point \cr
      \code{\link[spatstat.geom]{pointsOnLines}} \tab generate points evenly spaced
      along line segment \cr
      \code{\link[spatstat.random]{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window\cr
      \code{\link[spatstat.geom]{rlinegrid}} \tab generate a random array of parallel
      lines through a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link[spatstat.geom]{tess}} \tab create a tessellation \cr
      \code{\link[spatstat.geom]{quadrats}} \tab create a tessellation of rectangles\cr
      \code{\link[spatstat.geom]{hextess}} \tab create a tessellation of hexagons \cr
      \code{\link[spatstat.geom]{polartess}} \tab tessellation using polar coordinates \cr
      \code{\link[spatstat.geom]{quantess}} \tab quantile tessellation \cr
      \code{\link[spatstat.geom]{venn.tess}} \tab Venn diagram tessellation \cr 
      \code{\link[spatstat.geom]{dirichlet}} \tab compute Dirichlet-Voronoi tessellation of points\cr
      \code{\link[spatstat.geom]{delaunay}} \tab compute Delaunay triangulation of points\cr
      \code{\link[spatstat.geom]{as.tess}} \tab convert other data to a tessellation \cr
      \code{\link[spatstat.geom]{plot.tess}} \tab plot a tessellation \cr
      \code{\link[spatstat.geom]{tiles}} \tab extract all the tiles of a tessellation \cr
      \code{\link[spatstat.geom]{[.tess}} \tab extract some tiles of a tessellation \cr
      \code{\link[spatstat.geom]{[<-.tess}} \tab change some tiles of a tessellation \cr
      \code{\link[spatstat.geom]{intersect.tess}} \tab intersect two tessellations \cr
                              \tab or restrict a tessellation to a window \cr
      \code{\link[spatstat.geom]{chop.tess}} \tab subdivide a tessellation by a line \cr
      \code{\link[spatstat.random]{rpoislinetess}} \tab generate tessellation using Poisson line
      process \cr
      \code{\link[spatstat.geom]{tile.areas}}	\tab
      area of each tile in tessellation \cr
      \code{\link[spatstat.geom]{bdist.tiles}}	\tab
      boundary distance for each tile in tessellation \cr
      \code{\link[spatstat.geom]{connected.tess}} \tab find connected components of tiles \cr
      \code{\link[spatstat.geom]{shift.tess}} \tab shift a tessellation \cr
      \code{\link[spatstat.geom]{rotate.tess}} \tab rotate a tessellation \cr
      \code{\link[spatstat.geom]{reflect.tess}} \tab reflect about the origin \cr
      \code{\link[spatstat.geom]{flipxy.tess}} \tab reflect about the diagonal \cr
      \code{\link[spatstat.geom]{affine.tess}} \tab apply affine transformation 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link[spatstat.geom]{pp3}} \tab create a 3-D point pattern \cr
      \code{\link[spatstat.geom]{plot.pp3}} \tab plot a 3-D point pattern \cr
      \code{\link[spatstat.geom]{coords}} \tab extract coordinates \cr
      \code{\link[spatstat.geom]{as.hyperframe}} \tab extract coordinates \cr
      \code{\link[spatstat.geom]{subset.pp3}} \tab extract subset of 3-D point pattern \cr
      \code{\link[spatstat.geom]{unitname.pp3}} \tab name of unit of length \cr
      \code{\link[spatstat.geom]{npoints}} \tab count the number of points  \cr
      \code{\link[spatstat.random]{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link[spatstat.random]{rpoispp3}} \tab generate Poisson random points in 3-D \cr
      \code{\link[spatstat.explore]{envelope.pp3}} \tab generate simulation envelopes for
      3-D pattern \cr
      \code{\link[spatstat.geom]{box3}} \tab create a 3-D rectangular box \cr
      \code{\link[spatstat.geom]{as.box3}} \tab convert data to 3-D rectangular box \cr
      \code{\link[spatstat.geom]{unitname.box3}} \tab name of unit of length \cr
      \code{\link[spatstat.geom]{diameter.box3}} \tab diameter of box \cr
      \code{\link[spatstat.geom]{volume.box3}} \tab volume of box \cr
      \code{\link[spatstat.geom]{shortside.box3}} \tab shortest side of box \cr
      \code{\link[spatstat.geom]{eroded.volumes}} \tab volumes of erosions of box 
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link[spatstat.geom]{ppx}} \tab create a multidimensional space-time point pattern \cr
      \code{\link[spatstat.geom]{coords}} \tab extract coordinates \cr
      \code{\link[spatstat.geom]{as.hyperframe}} \tab extract coordinates \cr
      \code{\link[spatstat.geom]{subset.ppx}} \tab extract subset \cr
      \code{\link[spatstat.geom]{unitname.ppx}} \tab name of unit of length \cr
      \code{\link[spatstat.geom]{npoints}} \tab count the number of points  \cr
      \code{\link[spatstat.random]{runifpointx}} \tab generate uniform random points \cr
      \code{\link[spatstat.random]{rpoisppx}} \tab generate Poisson random points \cr
      \code{\link[spatstat.geom]{boxx}} \tab define multidimensional box  \cr
      \code{\link[spatstat.geom]{diameter.boxx}} \tab diameter of box \cr
      \code{\link[spatstat.geom]{volume.boxx}} \tab volume of box \cr
      \code{\link[spatstat.geom]{shortside.boxx}} \tab shortest side of box \cr
      \code{\link[spatstat.geom]{eroded.volumes.boxx}} \tab volumes of erosions of box 
    }
    
    \bold{Point patterns on a linear network}

    An object of class \code{"linnet"} represents a linear network
    (for example, a road network).

    \tabular{ll}{
      \code{\link[spatstat.linnet]{linnet}} \tab create a linear network \cr
      \code{\link[spatstat.linnet]{clickjoin}} \tab interactively join vertices in network \cr
      \code{spatstat.gui::iplot.linnet} \tab interactively plot network \cr
      \code{\link[spatstat.data]{simplenet}} \tab simple example of network \cr
      \code{\link[spatstat.linnet]{lineardisc}} \tab disc in a linear network \cr
      \code{\link[spatstat.linnet]{delaunayNetwork}} \tab network of Delaunay triangulation \cr
      \code{\link[spatstat.linnet]{dirichletNetwork}} \tab network of Dirichlet edges \cr
      \code{\link[spatstat.linnet]{methods.linnet}} \tab methods for \code{linnet} objects\cr
      \code{\link[spatstat.linnet]{vertices.linnet}} \tab nodes of network \cr
      \code{\link[spatstat.linnet]{joinVertices}} \tab join existing vertices in a network \cr
      \code{\link[spatstat.linnet]{insertVertices}} \tab insert new vertices at positions
      along a network \cr
      \code{\link[spatstat.linnet]{addVertices}} \tab add new vertices, extending a network \cr
      \code{\link[spatstat.linnet]{thinNetwork}} \tab remove vertices or lines from a network \cr
      \code{\link[spatstat.linnet]{repairNetwork}} \tab repair internal format \cr
      \code{\link[spatstat.linnet]{pixellate.linnet}} \tab approximate by pixel image
    }
    
    An object of class \code{"lpp"} represents a 
    point pattern on a linear network (for example,
    road accidents on a road network).
    
    \tabular{ll}{
      \code{\link[spatstat.linnet]{lpp}} \tab create a point pattern on a linear network \cr
      \code{\link[spatstat.linnet]{methods.lpp}} \tab methods for \code{lpp} objects \cr
      \code{\link[spatstat.linnet]{subset.lpp}} \tab method for \code{subset} \cr
      \code{\link[spatstat.linnet]{rpoislpp}} \tab simulate Poisson points on linear network \cr
      \code{\link[spatstat.linnet]{runiflpp}} \tab simulate random points on a linear network \cr
      \code{\link[spatstat.data]{chicago}} \tab Chicago crime data \cr
      \code{\link[spatstat.data]{dendrite}} \tab Dendritic spines data \cr
      \code{\link[spatstat.data]{spiders}} \tab Spider webs on mortar lines of brick wall 
    }
    
    \bold{Hyperframes}

    A hyperframe is like a data frame, except that the entries
    may be objects of any kind.

    \tabular{ll}{
      \code{\link[spatstat.geom]{hyperframe}} \tab create a hyperframe \cr
      \code{\link[spatstat.geom]{as.hyperframe}} \tab convert data to hyperframe \cr
      \code{\link[spatstat.geom]{plot.hyperframe}} \tab plot hyperframe \cr
      \code{\link[spatstat.geom]{with.hyperframe}} \tab evaluate expression using each row
      of hyperframe \cr
      \code{\link[spatstat.geom]{cbind.hyperframe}} \tab combine hyperframes by columns\cr
      \code{\link[spatstat.geom]{rbind.hyperframe}} \tab combine hyperframes by rows\cr
      \code{\link[spatstat.geom]{as.data.frame.hyperframe}} \tab convert hyperframe to
      data frame \cr
      \code{\link[spatstat.geom]{subset.hyperframe}} \tab method for \code{subset} \cr
      \code{\link[spatstat.geom]{head.hyperframe}} \tab first few rows of hyperframe \cr
      \code{\link[spatstat.geom]{tail.hyperframe}} \tab last few rows of hyperframe
    }
    
    \bold{Layered objects}

    A layered object represents data that should be plotted in
    successive layers, for example, a background and a foreground.

     \tabular{ll}{
      \code{\link[spatstat.geom]{layered}} \tab create layered object \cr
      \code{\link[spatstat.geom]{plot.layered}} \tab plot layered object\cr
      \code{\link[spatstat.geom]{[.layered}} \tab extract subset of layered object
    }
    
    \bold{Colour maps}

    A colour map is a mechanism for associating colours with data.
    It can be regarded as a function, mapping data to colours.
    Using a \code{colourmap} object in a plot command
    ensures that the mapping from numbers to colours is
    the same in different plots. 

     \tabular{ll}{
      \code{\link[spatstat.geom]{colourmap}} \tab create a colour map  \cr
      \code{\link[spatstat.geom]{plot.colourmap}} \tab plot the colour map only\cr
      \code{\link[spatstat.geom]{tweak.colourmap}} \tab alter individual colour values \cr
      \code{\link[spatstat.geom]{interp.colourmap}} \tab make a smooth transition
      between colours \cr
      \code{\link[spatstat.geom]{beachcolourmap}} \tab one special colour map
    }
}
  
\section{II. EXPLORATORY DATA ANALYSIS}{

  \bold{Inspection of data:}
  \tabular{ll}{
    \code{\link[spatstat.geom:summary.ppp]{summary}(X)} \tab
    print useful summary of point pattern \code{X}\cr
    \code{X} \tab
    print basic description of point pattern \code{X}  \cr
    \code{any(duplicated(X))} \tab
    check for duplicated points in pattern \code{X} \cr
    \code{spatstat.gui::istat(X)} \tab Interactive exploratory analysis \cr
    \code{spatstat.gui::View.ppp(X)} \tab spreadsheet-style viewer
  }

  \bold{Classical exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{clarkevans}} \tab Clark and Evans aggregation index \cr
    \code{\link[spatstat.explore]{fryplot}} \tab Fry plot \cr
    \code{\link[spatstat.explore]{miplot}} \tab Morisita Index plot
  }

  \bold{Smoothing:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{density.ppp}} \tab kernel smoothed density/intensity\cr
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.explore]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.explore]{bw.diggle}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{density.ppp}}\cr
    \code{\link[spatstat.explore]{bw.ppl}} \tab likelihood cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{density.ppp}}\cr
    \code{\link[spatstat.explore]{bw.CvL}} \tab Cronie-Van Lieshout bandwidth selection
                                for density estimation\cr
    \code{\link[spatstat.explore]{bw.scott}} \tab Scott's rule of thumb 
                                for density estimation\cr
    \code{\link[spatstat.explore]{bw.abram}} \tab Abramson's rule for adaptive bandwidths\cr
    \code{\link[spatstat.explore]{bw.relrisk}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{relrisk}} \cr
    \code{\link[spatstat.explore]{bw.smoothppp}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.explore]{Smooth.ppp}} \cr
    \code{\link[spatstat.explore]{bw.frac}} \tab bandwidth selection using window geometry\cr
    \code{\link[spatstat.explore]{bw.stoyan}} \tab Stoyan's rule of thumb for bandwidth
                               for \code{\link[spatstat.explore]{pcf}}
  }

  \bold{Modern exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.explore]{clusterset}} \tab Allard-Fraley feature detection  \cr
    \code{\link[spatstat.explore]{nnclean}} \tab Byers-Raftery feature detection  \cr
    \code{\link[spatstat.explore]{sharpen.ppp}} \tab Choi-Hall data sharpening \cr
    \code{\link[spatstat.explore]{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.explore]{rho2hat}} \tab Kernel estimate of effect of two covariates\cr
    \code{\link[spatstat.explore]{spatialcdf}} \tab Spatial cumulative distribution function\cr
    \code{\link[spatstat.explore]{roc}} \tab Receiver operating characteristic curve
  }

  \bold{Summary statistics for a point pattern:}
  Type \code{demo(sumfun)} for a demonstration of many
  of the summary statistics.
  \tabular{ll}{
    \code{\link[spatstat.geom]{intensity}} \tab Mean intensity \cr
    \code{\link[spatstat.geom]{quadratcount}} \tab Quadrat counts \cr
    \code{\link[spatstat.geom]{intensity.quadratcount}} \tab Mean intensity in quadrats \cr
    \code{\link[spatstat.explore]{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.explore]{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link[spatstat.explore]{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link[spatstat.explore]{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link[spatstat.explore]{Lest}} \tab Besag \eqn{L}-function\cr
    \code{\link[spatstat.explore]{Tstat}} \tab Third order \eqn{T}-function \cr
    \code{\link[spatstat.explore]{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link[spatstat.explore]{pcf}} \tab 	pair correlation function \cr
    \code{\link[spatstat.explore]{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Linhom}} \tab \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{pcfinhom}} \tab pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.explore]{Finhom}} \tab \eqn{F} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Ginhom}} \tab \eqn{G} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{Jinhom}} \tab \eqn{J} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localL}} \tab Getis-Franklin neighbourhood density function\cr
    \code{\link[spatstat.explore]{localK}} \tab neighbourhood K-function\cr
    \code{\link[spatstat.explore]{localpcf}} \tab local pair correlation function\cr
    \code{\link[spatstat.explore]{localKinhom}} \tab local \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localLinhom}} \tab local \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.explore]{localpcfinhom}} \tab local pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.explore]{Ksector}} \tab Directional \eqn{K}-function\cr
    \code{\link[spatstat.explore]{Kscaled}} \tab locally scaled \eqn{K}-function \cr
    \code{\link[spatstat.explore]{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link[spatstat.explore]{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link[spatstat.explore]{envelope}} \tab simulation envelopes for a summary
    function \cr
    \code{\link[spatstat.explore]{varblock}} \tab variances and confidence intervals\cr
                   \tab for a summary function \cr
    \code{\link[spatstat.explore]{lohboot}} \tab bootstrap for a summary function 
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.explore]{plot.fv}} \tab plot a summary function\cr
    \code{\link[spatstat.explore]{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link[spatstat.explore]{harmonise.fv}} \tab make functions compatible \cr
    \code{\link[spatstat.explore]{eval.fasp}} \tab evaluate any expression involving
    an array of functions\cr
    \code{\link[spatstat.explore]{with.fv}} \tab evaluate an expression for a 
    summary function\cr
    \code{\link[spatstat.explore]{Smooth.fv}} \tab apply smoothing to a summary function\cr
    \code{\link[spatstat.explore]{deriv.fv}} \tab calculate derivative of a summary function\cr
    \code{\link[spatstat.explore]{pool.fv}} \tab pool several estimates of a summary function\cr
    \code{\link[spatstat.geom]{nndist}} \tab nearest neighbour distances \cr
    \code{\link[spatstat.geom]{nnwhich}} \tab find nearest neighbours \cr
    \code{\link[spatstat.geom]{pairdist}} \tab distances between all pairs of points\cr
    \code{\link[spatstat.geom]{crossdist}} \tab distances between points in two patterns\cr
    \code{\link[spatstat.geom]{nncross}} \tab nearest neighbours between two point patterns \cr
    \code{\link[spatstat.geom]{exactdt}} \tab distance from any location to nearest data point\cr 
    \code{\link[spatstat.geom]{distmap}} \tab distance map image\cr
    \code{\link[spatstat.geom]{distfun}} \tab distance map function\cr
    \code{\link[spatstat.geom]{nnmap}} \tab nearest point image \cr
    \code{\link[spatstat.geom]{nnfun}} \tab nearest point function \cr
    \code{\link[spatstat.explore]{density.ppp}} \tab kernel smoothed density\cr
    \code{\link[spatstat.explore]{densityHeat.ppp}} \tab diffusion kernel smoothed density\cr
    \code{\link[spatstat.explore]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.explore]{sharpen.ppp}} \tab data sharpening  \cr
    \code{\link[spatstat.random]{rknn}} \tab theoretical distribution of nearest
    neighbour distance
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} such that \code{marks(X)} is a factor. 
  \tabular{ll}{
    \code{\link[spatstat.explore]{relrisk}} \tab kernel estimation of relative risk  \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan test of elevated risk  \cr
    \code{\link[spatstat.explore]{Gcross},\link[spatstat.explore]{Gdot},\link[spatstat.explore]{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{G[i,j], G[i.]} \cr
    \code{\link[spatstat.explore]{Kcross},\link[spatstat.explore]{Kdot}, \link[spatstat.explore]{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{K[i,j], K[i.]} \cr
    \code{\link[spatstat.explore]{Lcross},\link[spatstat.explore]{Ldot}} \tab
    multitype \eqn{L}-functions 
    \eqn{L_{ij}, L_{i\bullet}}{L[i,j], L[i.]} \cr
    \code{\link[spatstat.explore]{Jcross},\link[spatstat.explore]{Jdot},\link[spatstat.explore]{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{J[i,j],J[i.]} \cr
    \code{\link[spatstat.explore]{pcfcross}} \tab
    multitype pair correlation function \eqn{g_{ij}}{g[i,j]} \cr
    \code{\link[spatstat.explore]{pcfdot}} \tab
    multitype pair correlation function \eqn{g_{i\bullet}}{g[i.]} \cr
    \code{\link[spatstat.explore]{pcfmulti}} \tab
    general pair correlation function \cr
    \code{\link[spatstat.explore]{markconnect}} \tab
    marked connection function \eqn{p_{ij}}{p[i,j]} \cr
    \code{\link[spatstat.explore]{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link[spatstat.explore]{Iest}} \tab  multitype \eqn{I}-function\cr
    \code{\link[spatstat.explore]{Kcross.inhom},\link[spatstat.explore]{Kdot.inhom}} \tab
    inhomogeneous counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.explore]{Lcross.inhom},\link[spatstat.explore]{Ldot.inhom}} \tab
    inhomogeneous counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.explore]{pcfcross.inhom},\link[spatstat.explore]{pcfdot.inhom}} \tab
    inhomogeneous counterparts of \code{pcfcross}, \code{pcfdot} \cr
    \code{\link[spatstat.explore]{localKcross},\link[spatstat.explore]{localKdot}} \tab
    local counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.explore]{localLcross},\link[spatstat.explore]{localLdot}} \tab
    local counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.explore]{localKcross.inhom},\link[spatstat.explore]{localLcross.inhom}} \tab
    local counterparts of \code{Kcross.inhom}, \code{Lcross.inhom} 
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type. For numeric marks, there are the
  following functions:
  \tabular{ll}{
    \code{\link[spatstat.explore]{markmean}} \tab smoothed local average of marks \cr
    \code{\link[spatstat.explore]{markvar}} \tab smoothed local variance of marks \cr
    \code{\link[spatstat.explore]{markcorr}} \tab mark correlation function \cr
    \code{\link[spatstat.explore]{markcrosscorr}} \tab mark cross-correlation function \cr
    \code{\link[spatstat.explore]{markvario}} \tab mark variogram \cr
    \code{\link[spatstat.explore]{markmarkscatter}} \tab mark-mark scatterplot \cr
    \code{\link[spatstat.explore]{Kmark}} \tab mark-weighted \eqn{K} function \cr
    \code{\link[spatstat.explore]{Emark}} \tab mark independence diagnostic \eqn{E(r)} \cr
    \code{\link[spatstat.explore]{Vmark}} \tab mark independence diagnostic \eqn{V(r)} \cr
    \code{\link[spatstat.explore]{nnmean}} \tab nearest neighbour mean index \cr
    \code{\link[spatstat.explore]{nnvario}} \tab nearest neighbour mark variance index 
  }
  For marks of any type, there are the following:
  \tabular{ll}{
    \code{\link[spatstat.explore]{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link[spatstat.explore]{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link[spatstat.explore]{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link[spatstat.geom]{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link[spatstat.geom]{applynbd}} \tab apply function to every neighbourhood
    in a point pattern \cr
    \code{\link[spatstat.geom]{markstat}} \tab apply function to the marks of neighbours
    in a point pattern \cr
    \code{\link[spatstat.explore]{marktable}} \tab tabulate the marks of neighbours
    in a point pattern \cr
    \code{\link[spatstat.geom]{pppdist}} \tab find the optimal match between two point
    patterns
  }

  \bold{Summary statistics for a point pattern on a linear network:}

  These are for point patterns on a linear network (class \code{lpp}).
  For unmarked patterns:
  
  \tabular{ll}{
    \code{\link[spatstat.linnet]{linearK}} \tab
    \eqn{K} function on linear network \cr
    \code{\link[spatstat.linnet]{linearKinhom}} \tab
    inhomogeneous \eqn{K} function on linear network \cr
    \code{\link[spatstat.linnet]{linearpcf}} \tab
    pair correlation function on linear network \cr
    \code{\link[spatstat.linnet]{linearpcfinhom}} \tab
    inhomogeneous pair correlation on linear network
  }

  For multitype patterns:
  \tabular{ll}{
    \code{\link[spatstat.linnet]{linearKcross}} \tab
    \eqn{K} function between two types of points \cr
    \code{\link[spatstat.linnet]{linearKdot}} \tab
    \eqn{K} function from one type to any type \cr
    \code{\link[spatstat.linnet]{linearKcross.inhom}} \tab
    Inhomogeneous version of \code{\link[spatstat.linnet]{linearKcross}} \cr
    \code{\link[spatstat.linnet]{linearKdot.inhom}} \tab
    Inhomogeneous version of \code{\link[spatstat.linnet]{linearKdot}} \cr
    \code{\link[spatstat.linnet]{linearmarkconnect}} \tab
    Mark connection function  on linear network \cr
    \code{\link[spatstat.linnet]{linearmarkequal}} \tab
    Mark equality function on linear network \cr
    \code{\link[spatstat.linnet]{linearpcfcross}} \tab
    Pair correlation between two types of points \cr
    \code{\link[spatstat.linnet]{linearpcfdot}} \tab
    Pair correlation from one type to any type \cr
    \code{\link[spatstat.linnet]{linearpcfcross.inhom}} \tab
    Inhomogeneous version of \code{\link[spatstat.linnet]{linearpcfcross}} \cr
    \code{\link[spatstat.linnet]{linearpcfdot.inhom}} \tab
    Inhomogeneous version of \code{\link[spatstat.linnet]{linearpcfdot}} 
  }

  Related facilities:
  
  \tabular{ll}{
    \code{\link[spatstat.linnet]{pairdist.lpp}} \tab distances between pairs  \cr
    \code{\link[spatstat.linnet]{crossdist.lpp}} \tab distances between pairs \cr
    \code{\link[spatstat.linnet]{nndist.lpp}} \tab nearest neighbour distances  \cr
    \code{\link[spatstat.linnet]{nncross.lpp}} \tab nearest neighbour distances  \cr
    \code{\link[spatstat.linnet]{nnwhich.lpp}} \tab find nearest neighbours  \cr
    \code{\link[spatstat.linnet]{nnfun.lpp}} \tab find nearest data point  \cr
    \code{\link[spatstat.linnet]{density.lpp}} \tab kernel smoothing estimator of intensity  \cr
    \code{\link[spatstat.linnet]{densityHeat.lpp}} \tab diffusion kernel estimate  \cr
    \code{\link[spatstat.linnet]{distfun.lpp}} \tab distance transform  \cr
    \code{\link[spatstat.linnet]{envelope.lpp}} \tab simulation envelopes  \cr
    \code{\link[spatstat.linnet]{rpoislpp}} \tab simulate Poisson points on linear network \cr
    \code{\link[spatstat.linnet]{runiflpp}} \tab simulate random points on a linear network 
  }
  
  It is also possible to fit point process models to \code{lpp} objects.
  See Section IV.
  
  \bold{Summary statistics for a three-dimensional point pattern:}

  These are for 3-dimensional point pattern objects (class \code{pp3}).

  \tabular{ll}{
    \code{\link[spatstat.explore]{F3est}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.explore]{G3est}} \tab nearest neighbour function \eqn{G} \cr
    \code{\link[spatstat.explore]{K3est}} \tab \eqn{K}-function \cr
    \code{\link[spatstat.explore]{pcf3est}} \tab pair correlation function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.explore]{envelope.pp3}} \tab simulation envelopes \cr
    \code{\link[spatstat.geom]{pairdist.pp3}} \tab distances between all pairs of
    points \cr
    \code{\link[spatstat.geom]{crossdist.pp3}} \tab distances between points in
    two patterns \cr
    \code{\link[spatstat.geom]{nndist.pp3}} \tab nearest neighbour distances \cr
    \code{\link[spatstat.geom]{nnwhich.pp3}} \tab find nearest neighbours \cr
    \code{\link[spatstat.geom]{nncross.pp3}} \tab find nearest neighbours in another pattern
  }

  \bold{Computations for multi-dimensional point pattern:}

  These are for multi-dimensional space-time
  point pattern objects (class \code{ppx}).

  \tabular{ll}{
    \code{\link[spatstat.geom]{pairdist.ppx}} \tab distances between all pairs of
    points \cr
    \code{\link[spatstat.geom]{crossdist.ppx}} \tab distances between points in
    two patterns \cr
    \code{\link[spatstat.geom]{nndist.ppx}} \tab nearest neighbour distances \cr
    \code{\link[spatstat.geom]{nnwhich.ppx}} \tab find nearest neighbours
  }

  \bold{Summary statistics for random sets:}
  
  These work for point patterns (class \code{ppp}),
  line segment patterns (class \code{psp})
  or windows (class \code{owin}).
  
  \tabular{ll}{
    \code{\link[spatstat.explore]{Hest}} \tab spherical contact distribution \eqn{H} \cr
    \code{\link[spatstat.explore]{Gfox}} \tab Foxall \eqn{G}-function \cr
    \code{\link[spatstat.explore]{Jfox}} \tab Foxall \eqn{J}-function
  }
  
}

\section{III. MODEL FITTING (COX AND CLUSTER MODELS)}{
  
  Cluster process models (with homogeneous or inhomogeneous intensity)
  and Cox processes can be fitted by the function \code{\link[spatstat.model]{kppm}}.
  Its result is an object of class \code{"kppm"}.
  The fitted model can be printed, plotted, predicted, simulated
  and updated.

  \tabular{ll}{
    \code{\link[spatstat.model]{kppm}} \tab  Fit model\cr
    \code{\link[spatstat.model]{plot.kppm}} \tab  Plot the fitted model\cr
    \code{\link[spatstat.model]{summary.kppm}} \tab  Summarise the fitted model\cr
    \code{\link[spatstat.model]{fitted.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{predict.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{update.kppm}} \tab Update the model \cr
    \code{\link[spatstat.model]{improve.kppm}} \tab Refine the estimate of trend \cr
    \code{\link[spatstat.model]{simulate.kppm}} \tab Generate simulated realisations \cr
    \code{\link[spatstat.model]{vcov.kppm}} \tab Variance-covariance matrix of coefficients \cr
    \code{\link[spatstat.model:methods.kppm]{coef.kppm}}
    \tab Extract trend coefficients \cr
    \code{\link[spatstat.model:methods.kppm]{formula.kppm}}
    \tab Extract trend formula \cr
    \code{\link[spatstat.model]{parameters}} \tab Extract all model parameters \cr
    \code{\link[spatstat.model]{clusterfield.kppm}} \tab Compute offspring density \cr
    \code{\link[spatstat.model]{clusterradius.kppm}} \tab Radius of support of offspring density \cr
    \code{\link[spatstat.model]{Kmodel.kppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.model]{pcfmodel.kppm}} \tab Pair correlation of fitted model 
  }
  
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.
  
  The theoretical models can also be simulated,
  for any choice of parameter values,
  using \code{\link[spatstat.random]{rThomas}}, \code{\link[spatstat.random]{rMatClust}},
  \code{\link[spatstat.random]{rCauchy}}, \code{\link[spatstat.random]{rVarGamma}},
  and \code{\link[spatstat.random]{rLGCP}}.
  
  Lower-level fitting functions include:

  \tabular{ll}{
    \code{\link[spatstat.model]{lgcp.estK}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.model]{lgcp.estpcf}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.model]{thomas.estK}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.model]{thomas.estpcf}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.model]{matclust.estK}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.model]{matclust.estpcf}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.model]{cauchy.estK}} \tab fit a Neyman-Scott Cauchy cluster process \cr
   \code{\link[spatstat.model]{cauchy.estpcf}} \tab fit a Neyman-Scott Cauchy cluster process\cr
   \code{\link[spatstat.model]{vargamma.estK}} \tab fit a Neyman-Scott Variance Gamma process\cr
   \code{\link[spatstat.model]{vargamma.estpcf}} \tab fit a Neyman-Scott Variance Gamma process\cr
    \code{\link[spatstat.model]{mincontrast}} \tab low-level algorithm for fitting models
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
  \code{\link[spatstat.model]{ppm}}. Its result is an object of class \code{"ppm"}.
  
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
    \code{\link[spatstat.model]{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link[spatstat.model]{predict.ppm}}
    \tab   Compute the spatial trend and conditional intensity\cr
    \tab   of the fitted point process model \cr
    \code{\link[spatstat.model]{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link[spatstat.model]{parameters}} \tab Extract all model parameters\cr
    \code{\link[spatstat.model]{formula.ppm}} \tab Extract the trend formula\cr
    \code{\link[spatstat.model]{intensity.ppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{Kmodel.ppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.model]{pcfmodel.ppm}} \tab pair correlation of fitted model \cr
    \code{\link[spatstat.model]{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link[spatstat.model]{residuals.ppm}} \tab Compute point process residuals at quadrature points \cr
    \code{\link[spatstat.model]{update.ppm}} \tab Update the fit \cr
    \code{\link[spatstat.model]{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link[spatstat.model]{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.model]{simulate.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.model]{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link[spatstat.model]{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link[spatstat.model]{effectfun}} \tab Compute the fitted effect of one covariate\cr
    \code{\link[spatstat.model]{logLik.ppm}} \tab log-likelihood or log-pseudolikelihood\cr
    \code{\link[spatstat.model]{anova.ppm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.model]{model.frame.ppm}} \tab Extract data frame used to fit model  \cr
    \code{\link[spatstat.model]{model.images}} \tab Extract spatial data used to fit model  \cr
    \code{\link[spatstat.model]{model.depends}} \tab Identify variables in the model \cr
    \code{\link[spatstat.model]{as.interact}} \tab Interpoint interaction component of model \cr
    \code{\link[spatstat.model]{fitin}} \tab Extract fitted interpoint interaction \cr
    \code{\link[spatstat.model]{is.hybrid}} \tab Determine whether the model is a hybrid \cr
    \code{\link[spatstat.model]{valid.ppm}} \tab Check the model is a valid point process \cr
    \code{\link[spatstat.model]{project.ppm}} \tab Ensure the model is a valid point process 
  }
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.
  
  See \code{\link[spatstat.geom]{spatstat.options}} to control plotting of fitted model.
  
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
    \code{\link[spatstat.model]{Poisson}()} \tab		the Poisson point process\cr
    \code{\link[spatstat.model]{AreaInter}()}	 \tab Area-interaction process\cr
    \code{\link[spatstat.model]{BadGey}()} \tab	multiscale Geyer process\cr
    \code{\link[spatstat.model]{Concom}()} \tab	connected component interaction\cr
    \code{\link[spatstat.model]{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link[spatstat.model]{DiggleGatesStibbard}() } \tab Diggle-Gates-Stibbard potential \cr
    \code{\link[spatstat.model]{Fiksel}()}	 \tab Fiksel pairwise interaction process\cr
    \code{\link[spatstat.model]{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link[spatstat.model]{Hardcore}()}	 \tab Hard core process\cr
    \code{\link[spatstat.model]{HierHard}()}	 \tab Hierarchical multiype hard core process\cr
    \code{\link[spatstat.model]{HierStrauss}()}	 \tab Hierarchical multiype Strauss process\cr
    \code{\link[spatstat.model]{HierStraussHard}()}	 \tab
    Hierarchical multiype Strauss-hard core process\cr
    \code{\link[spatstat.model]{Hybrid}()}	 \tab Hybrid of several interactions\cr
    \code{\link[spatstat.model]{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link[spatstat.model]{MultiHard}()} \tab 		multitype hard core process \cr
    \code{\link[spatstat.model]{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link[spatstat.model]{MultiStraussHard}()} \tab 	multitype Strauss/hard core process \cr
    \code{\link[spatstat.model]{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link[spatstat.model]{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link[spatstat.model]{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link[spatstat.model]{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link[spatstat.model]{Penttinen}()} \tab	Penttinen pairwise interaction\cr
    \code{\link[spatstat.model]{SatPiece}()} \tab	Saturated pair model, piecewise  constant potential\cr
    \code{\link[spatstat.model]{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link[spatstat.model]{Softcore}()} \tab pairwise interaction, soft core potential\cr
    \code{\link[spatstat.model]{Strauss}()} \tab  Strauss process \cr
    \code{\link[spatstat.model]{StraussHard}()} \tab Strauss/hard core point process \cr
    \code{\link[spatstat.model]{Triplets}()} \tab Geyer triplets process
  }
  Note that it is also possible to combine several such interactions
  using \code{\link[spatstat.model]{Hybrid}}.
  
  \bold{Finer control over model fitting:}
  
  A quadrature scheme is represented by an object of
  class \code{"quad"}. To create a quadrature scheme, typically
  use \code{\link[spatstat.geom]{quadscheme}}.
  
  \tabular{ll}{
    \code{\link[spatstat.geom]{quadscheme}} \tab default quadrature scheme \cr
                             \tab using rectangular cells or Dirichlet cells\cr
    \code{\link[spatstat.geom]{pixelquad}}  \tab quadrature scheme based on image pixels \cr
    \code{\link[spatstat.geom]{quad}}       \tab create an object of class \code{"quad"}
  }
  
  To inspect a quadrature scheme:
  \tabular{ll}{
    \code{plot(Q)} \tab plot quadrature scheme \code{Q}\cr
    \code{print(Q)} \tab print basic information about quadrature scheme \code{Q}\cr
    \code{\link[spatstat.geom:summary.quad]{summary}(Q)} \tab summary of quadrature scheme \code{Q}
  }

  A quadrature scheme consists of data points, dummy points, and
  weights. To generate dummy points:
  \tabular{ll}{
    \code{\link[spatstat.geom]{default.dummy}} \tab default pattern of dummy points \cr
    \code{\link[spatstat.geom]{gridcentres}} \tab dummy points in a rectangular grid \cr
    \code{\link[spatstat.random]{rstrat}} \tab stratified random dummy pattern \cr
    \code{\link[spatstat.geom]{spokes}} \tab radial pattern of dummy points  \cr
    \code{\link[spatstat.geom]{corners}} \tab dummy points at corners of the window
  }
  
  To compute weights:
  \tabular{ll}{
    \code{\link[spatstat.geom]{gridweights}} \tab quadrature weights by the grid-counting rule  \cr
    \code{\link[spatstat.geom]{dirichletWeights}} \tab quadrature weights are
    Dirichlet tile areas
  }

  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link[spatstat.model]{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.model]{simulate.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.explore]{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }

  \bold{Point process models on a linear network:}

  An object of class \code{"lpp"} represents a pattern of points on
  a linear network. Point process models can also be fitted to these
  objects. Currently only Poisson models can be fitted.

  \tabular{ll}{
    \code{\link[spatstat.linnet]{lppm}} \tab point process model on linear network \cr
    \code{\link[spatstat.linnet]{anova.lppm}} \tab analysis of deviance for \cr
    \tab point process model on linear network \cr
    \code{\link[spatstat.linnet]{envelope.lppm}} \tab simulation envelopes for \cr
    \tab point process model on linear network \cr
    \code{\link[spatstat.linnet]{fitted.lppm}} \tab fitted intensity values \cr
    \code{\link[spatstat.linnet]{predict.lppm}} \tab model prediction on linear network \cr
    \code{\link[spatstat.linnet]{linim}} \tab pixel image on linear network \cr
    \code{\link[spatstat.linnet]{plot.linim}} \tab plot a pixel image on linear network \cr
    \code{\link[spatstat.linnet]{eval.linim}} \tab evaluate expression involving images \cr
    \code{\link[spatstat.linnet]{linfun}} \tab function defined on linear network \cr
    \code{\link[spatstat.linnet]{methods.linfun}} \tab conversion facilities
  }
}


\section{V. MODEL FITTING (DETERMINANTAL POINT PROCESS MODELS)}{

  Code for fitting \emph{determinantal point process models} has 
  recently been added to \pkg{spatstat}.

  For information, see the help file for \code{\link[spatstat.model]{dppm}}.
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
  \code{\link[spatstat.model]{slrm}}. Its result is an object of class \code{"slrm"}.
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
    \code{\link[spatstat.model]{anova.slrm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.model]{coef.slrm}}  \tab Extract fitted coefficients \cr
    \code{\link[spatstat.model]{vcov.slrm}}  \tab Variance-covariance matrix of fitted coefficients \cr
    \code{\link[spatstat.model]{fitted.slrm}} \tab Compute fitted probabilities or
    intensity \cr
    \code{\link[spatstat.model]{logLik.slrm}}   \tab Evaluate loglikelihood of fitted
    model \cr
    \code{\link[spatstat.model]{plot.slrm}}    \tab Plot fitted probabilities or
    intensity \cr
    \code{\link[spatstat.model]{predict.slrm}} \tab Compute predicted probabilities or
    intensity with new data \cr
    \code{\link[spatstat.model]{simulate.slrm}} \tab Simulate model
  }
  
  There are many other undocumented methods for this class,
  including methods for \code{print}, \code{update}, \code{formula}
  and \code{terms}. Stepwise model selection is
  possible using \code{step} or \code{stepAIC}.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.
}


\section{VII. SIMULATION}{

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}

  \tabular{ll}{
    \code{\link[spatstat.random]{runifpoint}} \tab
    generate \eqn{n} independent uniform random points \cr
    \code{\link[spatstat.random]{rpoint}} \tab
    generate \eqn{n} independent random points \cr
    \code{\link[spatstat.random]{rmpoint}} \tab
    generate \eqn{n} independent multitype random points \cr
    \code{\link[spatstat.random]{rpoispp}} \tab
    simulate the (in)homogeneous Poisson point process \cr
    \code{\link[spatstat.random]{rmpoispp}} \tab
    simulate the (in)homogeneous multitype Poisson point process \cr
    \code{\link[spatstat.random]{runifdisc}} \tab
    generate \eqn{n} independent uniform random points in disc\cr
    \code{\link[spatstat.random]{rstrat}} \tab
    stratified random sample of points \cr
    \code{\link[spatstat.geom]{rsyst}} \tab
    systematic random sample (grid) of points \cr
    \code{\link[spatstat.random]{rMaternI}}  \tab
    simulate the \Matern Model I inhibition process\cr
    \code{\link[spatstat.random]{rMaternII}} \tab
    simulate the \Matern Model II inhibition process\cr
    \code{\link[spatstat.random]{rSSI}} \tab
    simulate Simple Sequential Inhibition process\cr
    \code{\link[spatstat.random]{rHardcore}} \tab
    simulate hard core process (perfect simulation)\cr
    \code{\link[spatstat.random]{rStrauss}} \tab
    simulate Strauss process (perfect simulation)\cr
    \code{\link[spatstat.random]{rStraussHard}} \tab
    simulate Strauss-hard core process (perfect simulation)\cr
    \code{\link[spatstat.random]{rDiggleGratton}} \tab
    simulate Diggle-Gratton process (perfect simulation)\cr
    \code{\link[spatstat.random]{rDGS}} \tab
    simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
    \code{\link[spatstat.random]{rPenttinen}} \tab
    simulate Penttinen process (perfect simulation)\cr
    \code{\link[spatstat.random]{rNeymanScott}} \tab
    simulate a general Neyman-Scott process\cr
    \code{\link[spatstat.random]{rMatClust}} \tab
    simulate the \Matern Cluster process\cr
    \code{\link[spatstat.random]{rThomas}} \tab
    simulate the Thomas process  \cr
    \code{\link[spatstat.random]{rLGCP}} \tab
    simulate the log-Gaussian Cox process  \cr
    \code{\link[spatstat.random]{rGaussPoisson}}  \tab
    simulate the Gauss-Poisson cluster process\cr
    \code{\link[spatstat.random]{rCauchy}} \tab
    simulate Neyman-Scott process with Cauchy clusters \cr
    \code{\link[spatstat.random]{rVarGamma}} \tab
    simulate Neyman-Scott process with Variance Gamma clusters \cr
    \code{\link[spatstat.random]{rcell}} \tab
    simulate the Baddeley-Silverman cell process  \cr
    \code{\link[spatstat.random]{runifpointOnLines}} \tab
    generate \eqn{n} random points along specified line segments \cr
    \code{\link[spatstat.random]{rpoisppOnLines}} \tab
    generate Poisson random points along specified line segments 
  }
      
  \bold{Resampling a point pattern:}

  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.geom]{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning 
  }
  
  See also \code{\link[spatstat.explore]{varblock}} for estimating the variance
  of a summary statistic by block resampling, and
  \code{\link[spatstat.explore]{lohboot}} for another bootstrap technique.
  
  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.

   Cluster process models 
   are fitted by the function \code{\link[spatstat.model]{kppm}} yielding an
   object of class \code{"kppm"}. To generate one or more simulated
   realisations of this fitted model, use 
   \code{\link[spatstat.model]{simulate.kppm}}.

   Gibbs point process models 
   are fitted by the function \code{\link[spatstat.model]{ppm}} yielding an
   object of class \code{"ppm"}. To generate a simulated
   realisation of this fitted model, use \code{\link[spatstat.random]{rmh}}.
   To generate one or more simulated realisations of the fitted model,
   use \code{\link[spatstat.model]{simulate.ppm}}.

   \bold{Other random patterns:}

   \tabular{ll}{
     \code{\link[spatstat.geom]{rlinegrid}} \tab
     generate a random array of parallel lines through a window \cr
     \code{\link[spatstat.random]{rpoisline}} \tab
     simulate the Poisson line process within a window \cr
     \code{\link[spatstat.random]{rpoislinetess}} \tab
     generate random tessellation using Poisson line process \cr
     \code{\link[spatstat.random]{rMosaicSet}} \tab
     generate random set by selecting some tiles of a tessellation \cr
     \code{\link[spatstat.random]{rMosaicField}} \tab
     generate random pixel image by assigning random values
     in each tile of a tessellation
   }

   \bold{Simulation-based inference}

   \tabular{ll}{
    \code{\link[spatstat.explore]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.explore]{bits.envelope}} \tab critical envelope for balanced
    two-stage Monte Carlo test \cr
    \code{\link[spatstat.model]{qqplot.ppm}} \tab diagnostic plot for interpoint
    interaction \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.explore]{studpermu.test}} \tab studentised permutation test\cr
    \code{\link[spatstat.explore]{segregation.test}} \tab test of segregation of types 
  }
}


\section{VIII. TESTS AND DIAGNOSTICS}{

  \bold{Hypothesis tests:}

  \tabular{ll}{
    \code{\link[spatstat.explore]{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link[spatstat.explore]{clarkevans.test}} \tab Clark and Evans test \cr
    \code{\link[spatstat.explore]{cdf.test}} \tab Spatial distribution goodness-of-fit test\cr
    \code{\link[spatstat.explore]{berman.test}} \tab Berman's goodness-of-fit tests\cr
    \code{\link[spatstat.explore]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.explore]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.explore]{dclf.test}} \tab Diggle-Cressie-Loosmore-Ford test \cr
    \code{\link[spatstat.explore]{mad.test}} \tab Mean Absolute Deviation test \cr
    \code{\link[spatstat.model]{anova.ppm}} \tab Analysis of Deviance for
    point process models 
  }

  More recently-developed tests:
  
  \tabular{ll}{
    \code{\link[spatstat.explore]{dg.test}} \tab Dao-Genton test \cr
    \code{\link[spatstat.explore]{bits.test}} \tab Balanced independent two-stage test \cr
    \code{\link[spatstat.explore]{dclf.progress}} \tab Progress plot for DCLF test \cr
    \code{\link[spatstat.explore]{mad.progress}} \tab Progress plot for MAD test \cr
  }

\bold{Sensitivity diagnostics:}

  Classical measures of model sensitivity such as leverage and influence
  have been adapted to point process models.
  
    \tabular{ll}{
      \code{\link[spatstat.model]{leverage.ppm}} \tab Leverage for point process model\cr
      \code{\link[spatstat.model]{influence.ppm}} \tab Influence for point process model\cr
      \code{\link[spatstat.model]{dfbetas.ppm}} \tab Parameter influence\cr
      \code{\link[spatstat.model]{dffit.ppm}} \tab Effect change diagnostic
    }
  
  \bold{Diagnostics for covariate effect:}

  Classical diagnostics for covariate effects have been adapted to
  point process models.

  \tabular{ll}{
    \code{\link[spatstat.model]{parres}} \tab Partial residual plot\cr
    \code{\link[spatstat.model]{addvar}} \tab Added variable plot \cr
    \code{\link[spatstat.explore]{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.explore]{rho2hat}} \tab Kernel estimate of covariate effect
    (bivariate) 
  }
  
  \bold{Residual diagnostics:}
  
  Residuals for a fitted point process model, and diagnostic plots
  based on the residuals, were introduced in Baddeley et al (2005) and
  Baddeley, Rubak and \Moller (2011).
  
  Type \code{demo(diagnose)}
  for a demonstration of the diagnostics features.

  \tabular{ll}{
    \code{\link[spatstat.model]{diagnose.ppm}} \tab diagnostic plots for spatial trend\cr
    \code{\link[spatstat.model]{qqplot.ppm}} \tab diagnostic Q-Q plot for interpoint interaction\cr
    \code{\link[spatstat.data]{residualspaper}} \tab examples from Baddeley et al (2005) \cr
    \code{\link[spatstat.model]{Kcom}} \tab model compensator of \eqn{K} function \cr
    \code{\link[spatstat.model]{Gcom}} \tab model compensator of \eqn{G} function \cr
    \code{\link[spatstat.model]{Kres}} \tab score residual of \eqn{K} function \cr
    \code{\link[spatstat.model]{Gres}} \tab score residual of \eqn{G} function \cr
    \code{\link[spatstat.model]{psst}} \tab pseudoscore residual of summary function \cr
    \code{\link[spatstat.model]{psstA}} \tab pseudoscore residual of empty space function \cr
    \code{\link[spatstat.model]{psstG}} \tab pseudoscore residual of \eqn{G} function \cr
    \code{\link[spatstat.model]{compareFit}} \tab compare compensators of several fitted
    models
  }


  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.geom]{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning  
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
  Workshop notes, Version 4.1.
  Online technical publication, CSIRO.
  \url{https://research.csiro.au/software/wp-content/uploads/sites/6/2015/02/Rspatialcourse_CMIS_PDF-Standard.pdf}

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
  
  Baddeley, A. and Turner, R. (2005a)
  Spatstat: an R package for analyzing spatial point patterns.
  \emph{Journal of Statistical Software} \bold{12}:6, 1--42.
  \code{DOI: 10.18637/jss.v012.i06}. 

  Baddeley, A. and Turner, R. (2005b)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.

  Baddeley, A., Turner, R.,
  \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., Rubak, E. and \Moller, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.

  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  \emph{Journal of Statistical Software} \bold{55}:11, 1--43.
  \url{https://www.jstatsoft.org/v55/i11/}

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
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Yongtao Guan,
  Ute Hahn,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  Additional contributions and suggestions from
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
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
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \Lucia Cobo Sanchez,
  Jean-\Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
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
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kassel Hingee,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \Frederic Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Bethany Macdonald,
  Ben Madin,
  Martin Maechler,
  Daniel Manrique-\ifelse{latex}{\out{Casta{\~n}o}}{Castano},
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Suman Rakshit,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Yonatan Rosen,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \Francois \ifelse{latex}{\out{S\'{e}m\'{e}curbe}}{Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

