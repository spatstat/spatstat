spatstat
========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat)](http://cran.r-project.org/web/packages/spatstat)
[![Travis-CI Build Status](https://travis-ci.org/spatstat/spatstat.png?branch=master)](https://travis-ci.org/spatstat/spatstat)
[![codecov.io](https://codecov.io/github/spatstat/spatstat/coverage.svg)](https://codecov.io/github/spatstat/spatstat?branch=covr)
[![appveyor build status](https://ci.appveyor.com/api/projects/status/github/spatstat/spatstat)](https://ci.appveyor.com/api/projects/status/github/spatstat/spatstat)

## What is spatstat?

`spatstat` is a family of R packages for analysing 
spatial point pattern data (and other kinds of spatial data).
See the website [www.spatstat.org](http://www.spatstat.org)
or read the [book](http://book.spatstat.org).

## spatstat has been split into a family of packages

Originally there was a single package called `spatstat`.
It grew so large (150,000 lines of code) that CRAN required us
to split it into pieces. 

Almost all of the code in the original `spatstat` has been 
 placed into a family of sub-packages:

| Sub-package | CRAN page | GitHub repository | Description |
| ----------  | --------- | ----------------- | ----------  |
| `spatstat.utils` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.utils)](http://cran.r-project.org/web/packages/spatstat.utils) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.utils)](https://github.com/spatstat/spatstat.utils) | Basic utilities |
| `spatstat.data` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.data)](http://cran.r-project.org/web/packages/spatstat.data) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.data)](https://github.com/spatstat/spatstat.data) | Datasets |
| `spatstat.sparse` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.sparse)](http://cran.r-project.org/web/packages/spatstat.sparse) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.sparse)](https://github.com/spatstat/spatstat.sparse) | Sparse arrays |
| `spatstat.geom` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.geom)](http://cran.r-project.org/web/packages/spatstat.geom) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.geom)](https://github.com/spatstat/spatstat.geom) | Spatial data classes; geometrical operations |
| `spatstat.core` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.core)](http://cran.r-project.org/web/packages/spatstat.core) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.core)](https://github.com/spatstat/spatstat.core) | Data analysis and modelling of spatial data |
| `spatstat.linnet` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.linnet)](http://cran.r-project.org/web/packages/spatstat.linnet) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.linnet)](https://github.com/spatstat/spatstat.linnet) | Spatial analysis on a linear network |
| `spatstat` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat)](http://cran.r-project.org/web/packages/spatstat) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat)](https://github.com/spatstat/spatstat) | Umbrella package |

Click the green badge to visit the CRAN page which contains the current
release of each sub-package. Click the blue badge to visit the GitHub repository
for the current development version of the sub-package. 

There still exists a package called `spatstat`, which is now an
**umbrella package** that requires all the sub-packages listed above.
When you install the new `spatstat`, all the sub-packages listed above will
be installed. When you load the `spatstat` package in an R session,
all the sub-packages listed above will be loaded or imported.

## Extension packages

Additionally there are **extension packages** which contain additional
functionality. These packages are not automatically installed or loaded;
the user must do that if these extra features are desired.

| Extension package | CRAN page | GitHub repository | Description |
| ----------------  | --------- | ----------------- | ----------  |
| `spatstat.gui` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.gui)](http://cran.r-project.org/web/packages/spatstat.gui)  | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.gui)](https://github.com/spatstat/spatstat.gui) | Graphical interface |
| `spatstat.Knet` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.Knet)](http://cran.r-project.org/web/packages/spatstat.Knet) | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.Knet)](https://github.com/spatstat/spatstat.Knet) | linear networks |
| `spatstat.local` | [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.local)](http://cran.r-project.org/web/packages/spatstat.local) | [![GitHub R package version](https://img.shields.io/github/r-package/v/baddstats/spatstat.local)](https://github.com/baddstats/spatstat.local) | Local (geographically weighted) models |
| `spatstat.sphere` | Not yet published | [![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.sphere)](https://github.com/spatstat/spatstat.sphere) | Spherical data |

## Family portrait 

The pink box marked `spatstat` contains all the code that will be
installed when you install the `spatstat` umbrella package, and loaded
or imported when you load the `spatstat` umbrella package.

The blue boxes are extension packages which must be installed and loaded
separately.

![Spatstat pieces](RepoStuff/newspatstat.jpg)
