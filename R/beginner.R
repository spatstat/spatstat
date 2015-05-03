#
#  beginner.R
#
# Helpful information for beginners
#
#  $Revision: 1.2 $  $Date: 2014/10/24 00:22:30 $
#

print.autoexec <- function(x, ...) { x() }

beginner <- function(package="spatstat") {
  package <- as.character(substitute(package))
  RShowDoc("BEGINNER.txt", type="txt", package=package)
  return(invisible(NULL))
}

class(beginner) <- "autoexec"

foo <- local({
  fooText <- paste0("Error: object 'foo' not found.\n\n",
                    "'foo' is not a defined variable or function.\n",
                    "It is a placeholder name, which serves only to ",
                    "demonstrate a concept. It represents the name of ",
                    "any desired object or function. ", 
                    "Other placeholder names popular with computer scientists ",
                    "are 'bar', 'foobar', 'qux' and 'mork'.")

  foo <- function() {
    splat(fooText) 
    return(invisible(NULL))
  }
  class(foo) <- "autoexec"
  foo
})

plot.foo <- function(x, ...) foo()
