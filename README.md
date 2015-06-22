spatstat
========

This repository holds a copy of the _current development version_ 
of the contributed R-package `spatstat`.

This development version is more recent than the official *release* 
of `spatstat` on CRAN. Each official release of `spatstat` has a
version number like `1.2-3` while the development version has a 
version number like `1.2-3.004`. Official releases occur every 8 weeks
(the minimum time permitted by CRAN policies) while the development code
is updated almost every day.

## Installation

The easiest way to install the development version of `spatstat` 
from github is through the `devtools` package:

```R
require(devtools)
install_github('spatstat/spatstat')
```

If you don't have `devtools` installed you should first run

```R
install.packages('devtools')
```

## Bug reports 

Users of `spatstat` are encouraged to report bugs here 
(press *issue* in the menu on the right to start a new bug report
or feature request).

## Making your own changes

Feel free to fork `spatstat`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

However, this repository is only a copy of the development code, so 
your pull request may not be implemented until the next official release.

