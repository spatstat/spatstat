spatstat
========

This repository holds a *copy* of the _current development version_ 
of the contributed R-package `spatstat`.

This development version is more recent than the current *release* of `spatstat` on CRAN.
A release of `spatstat` has a version number like 1.2-3 while a development version has a 
version number like 1.2-3.002.

Users of `spatstat` are encouraged to report bugs and make feature
requests here (press *issue* in the menu on the right to start a new bug report or feature request).

Feel free to fork `spatstat` and make pull requests. 
However, we are only in the process of moving to git and github, 
so this is not the actual development code at the moment and
your pull request may not be implemented until the next official release.

## Installation

The easiest way to install the github version of `spatstat` is through the `devtools` package:

```R
require(devtools)
install_github('spatstat/spatstat')
```

If you don't have `devtools` installed you should first run

```R
install.packages('devtools')
```
