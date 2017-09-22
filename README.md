spatstat
========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat)](http://cran.r-project.org/web/packages/spatstat)
[![Research software impact](http://depsy.org/api/package/cran/spatstat/badge.svg)](http://depsy.org/package/r/spatstat)
[![Travis-CI Build Status](https://travis-ci.org/spatstat/spatstat.png?branch=master)](https://travis-ci.org/spatstat/spatstat)
[![codecov.io](https://codecov.io/github/spatstat/spatstat/coverage.svg?branch=covr)](https://codecov.io/github/spatstat/spatstat?branch=covr)

## This is the development version

This repository holds a copy of the _current development version_ 
of the contributed R-package `spatstat`.

This development version is more recent than the official *release* 
of `spatstat` on CRAN. Each official release of `spatstat` has a
version number like `1.2-3` while the development version has a 
version number like `1.2-3.004`. Official releases occur every 8 weeks
(the minimum time permitted by CRAN policies) while the development code
is updated almost every day. 

## Where is the official release?

For the most recent **official release** of 
`spatstat`, see the [CRAN page](https://cran.r-project.org/web/packages/spatstat). 

## spatstat is now split into several packages

Recently we have started the process of splitting `spatstat` into several
packages (as mandated by CRAN, because `spatstat` is very large).

Your existing code will still work:
typing `library(spatstat)` will still give you
access to all the functions in `spatstat` that you know from previous versions.

However, messages from `R` about the installation and loading of the package
will now show that `spatstat` consists of several pieces.
Currently there are three pieces:

  . `spatstat`: contains the main functionality of the `spatstat` family.

  . `spatstat.data`: contains the datasets for the `spatstat` family.
  The current development version of `spatstat.data` is
  [here](https://github.com/spatstat/spatstat.data).

  . `spatstat.utils`: utility functions originally included in `spatstat`
  which are now accessible as a separate package.
  The current development version of `spatstat.utils` is
  [here](https://github.com/spatstat/spatstat.utils).

When you type `library(spatstat)` this will load
the main `spatstat` library
and the `spatstat.data` library, 
and will also *import* the `spatstat.utils` library. This means that
`spatstat.utils` functions can be used by `spatstat` but cannot be accessed by
the user. To access these utility functions directly, you need to type
`library(spatstat.utils)`.

### Extension packages

There are also *extension packages* which provide additional capabilities
and must be loaded explicitly when you need them. 
Currently there are two extension packages:

   . [spatstat.local](https://github.com/baddstats/spatstat.local)
   for local model-fitting, 

   . [spatstat.sphere](https://github.com/spatstat/spatstat.sphere)
   for analysing point patterns on a sphere.

## Installation

### Installing the official release

To install the official release of `spatstat` from CRAN, start `R` and type

```R
install.packages('spatstat')
```

### Installing the development version

The easiest way to install the development version of `spatstat` 
from github is through the `remotes` package. Start `R` and type

```R
require(remotes)
install_github('spatstat/spatstat.utils')
install_github('spatstat/spatstat.data')
install_github('spatstat/spatstat')
```

If you don't have `remotes` installed you should first run

```R
install.packages('remotes')
```

## Bug reports 

Users are encouraged to report bugs here.
Go to 
[issues](https://github.com/spatstat/spatstat/issues) in the menu above, 
and press *new issue* to start a new bug report, documentation correction
or feature request.

Please do not post *questions* on the Issues page
because it's too clunky for correspondence.

## Questions about spatstat

For questions about `spatstat`, first check 
the question-and-answer website [stackoverflow](http://stackoverflow.com/questions/tagged/spatstat).
If your question is not listed,
you can either post your question at stackoverflow, or
email the authors.

## Making your own changes

Feel free to fork `spatstat`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

However, this repository is only a copy of the development code, so 
your pull request may not be implemented until the next official release.

