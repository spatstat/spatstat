#
#  beginner.R
#
# Helpful information for beginners
#
#  $Revision: 1.1 $  $Date: 2013/08/02 05:40:16 $
#

beginner <- function(package="spatstat") {
  package <- as.character(substitute(package))
  RShowDoc("BEGINNER.txt", type="txt", package=package)
  return(invisible(NULL))
}

class(beginner) <- "autoexec"

print.autoexec <- function(x, ...) { x() }
