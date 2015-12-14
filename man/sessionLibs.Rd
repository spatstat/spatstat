\name{sessionLibs}
\alias{sessionLibs}
\title{
  Print Names and Version Numbers of Libraries Loaded
}
\description{
  Prints the names and version numbers
  of libraries currently loaded by the user.
}
\usage{
  sessionLibs()
}
\details{
  This function prints a list of the libraries loaded by the user
  in the current session, giving just their name and version number.
  It obtains this information from \code{\link[utils]{sessionInfo}}.

  This function is not needed in an interactive \R session
  because the package startup messages will usually provide this
  information.
  
  Its main use is in an \code{\link{Sweave}} script, where it is needed
  because the package startup messages are not printed. 
}
\value{
  Null.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link[utils]{sessionInfo}}
}
\keyword{data}
