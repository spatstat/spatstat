spatstat
========

This repository holds the currently released source code for the
R-package `spatstat`, but we aim at making it contain a development
version of `spatstat` in the "near" future.

Users of `spatstat` are encouraged to report bugs and make feature
requests here (press issue in the menu on the right to start a new bug report or feature request).

Feel free to fork `spatstat` and make pull
requests. However, we are only in the process of moving to git and
github, so this is not the actual development code at the moment and
your pull request may have to wait until the next official release.

## Installation

The easiest way to install the github version of `spatstat` is through the `devtools` package:

```R
require(devtools)
install_github('spatstat', username = 'spatstat')
```

If you don't have `devtools` installed you should first run

```R
install.packages('devtools')
```
