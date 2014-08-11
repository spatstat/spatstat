#
#     options.R
#
#     Spatstat options and other internal states
#
#    $Revision: 1.53 $   $Date: 2014/04/25 07:15:10 $
#
#

.spEnv <- new.env()

putSpatstatVariable <- function(name, value) {
  assign(name, value, envir=.spEnv)
}
getSpatstatVariable <- function(name) {
  get(name, envir=.spEnv)
}

putSpatstatVariable("Spatstat.Options", list())
putSpatstatVariable("Spatstat.ProgressBar", NULL)
putSpatstatVariable("Spatstat.ProgressData", NULL)
putSpatstatVariable("warnedkeys", character(0))

warn.once <- function(key, ...) {
  warned <- getSpatstatVariable("warnedkeys")
  if(!(key %in% warned)) {
    warning(paste(...), call.=FALSE)
    putSpatstatVariable("warnedkeys", c(warned, key))
  }
  return(invisible(NULL))
}

".Spat.Stat.Opt.Table" <-
  list(
       allow.logi.influence=list(
         ## whether leverage/influence calculations are permitted
         ## on a fitted model with method='logi'
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       checkpolygons = list(
         ## superseded
         default=FALSE,
         check=function(x) {
           warning("spatstat.options('checkpolygons') will be ignored in future versions of spatstat", call.=FALSE)
           return(is.logical(x) && length(x) == 1)
         },
         valid="a single logical value"
         ),
       checksegments = list(
         ## default value of 'check' for psp objects
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       closepairs.newcode=list(
         ## use new code for 'closepairs'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       crossing.psp.useCall=list(
         ## use new code for 'crossing.psp'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       crosspairs.newcode=list(
         ## use new code for 'crosspairs'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       densityC=list(
         ## use C routines for 'density.ppp'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       dupC = list(
         ## value of 'DUP' for .C calls
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       exactdt.checks.data=list(
         ## whether 'exactdt' checks validity of return value
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       expand=list(
         ## default area expansion factor
         default=2,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && x > 1
         },
         valid="a single numeric value, greater than 1"
       ),
       fasteval=list(
         ## whether to use 'fasteval' code if available
         default="on",
         check=function(x) { x %in% c("off", "on", "test") },
         valid="one of the strings \'off\', \'on\' or \'test\'"
       ),
       fastK.lgcp=list(
         ## whether to cut a few corners in 'lgcp.estK'
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       fixpolygons = list(
         ## whether to repair polygons automatically
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       gpclib=list(
         ## defunct!
         default=FALSE,
         check=function(x) {
           message("gpclib is no longer needed")
           return(TRUE)
         },
         valid="a single logical value"
         ),
       huge.npoints=list(
         ## threshold to trigger a warning from rpoispp 
         default=1e6,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       image.colfun=list(
         ## default colour scheme for plot.im
         default=function(n){topo.colors(n)},
         check=function(x) {
           if(!is.function(x) || length(formals(x)) == 0) return(FALSE)
           y <- x(42)
           if(length(y) != 42 || !is.character(y)) return(FALSE)
           z <- try(col2rgb(y), silent=TRUE)
           return(!inherits(z, "try-error"))
         },
         valid="a function f(n) that returns character strings, interpretable as colours"
         ),
       Kcom.remove.zeroes=list(
         ## whether Kcom removes zero distances
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       maxedgewt=list(
         ## maximum edge correction weight 
         default=100,
         check=function(x){
           is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1
         },
         valid="a finite numeric value, not less than 1"
       ),
       maxmatrix=list(
         ## maximum size of matrix of pairs of points in mpl.R
         default=2^24, # 16,777,216
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       monochrome = list(
         ## switch for monochrome colour scheme
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       n.bandwidth=list(
         ## number of values of bandwidth to try in bandwidth selection
         default=32,
         check=function(x) {
           is.numeric(x) && (length(x) == 1) && (x == ceiling(x)) && (x > 2)
         },
         valid="a single integer, greater than 2"
       ),
       ndummy.min=list(
         ## minimum grid size for dummy points
         default=32,
         check=function(x) {
           is.numeric(x) && length(x) <= 2 && all(x == ceiling(x)) && all(x > 1)
         },
         valid="a single integer or a pair of integers, greater than 1"
       ),
       ngrid.disc=list(
         ## number of grid points used to calculate area in area-interaction
         default=128,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1
         },
         valid="a single integer, greater than 1"
       ),
       npixel=list(
         ## default pixel dimensions
         default=128,
         check=function(x){
           is.numeric(x) && (length(x) %in% c(1,2)) && is.finite(x) &&
           all(x == ceiling(x)) && all(x > 1) 
         },
         valid="an integer, or a pair of integers, greater than 1"
        ),
       nvoxel=list(
         ## default total number of voxels
         default=2^22,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 2^12
         },
         valid="a single integer, greater than 2^12"
       ),
       old.morpho.psp=list(
         ## use old code for morphological operations
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       par.binary=list(
         ## default graphics parameters for masks
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.persp=list(
         ## default graphics parameters for 'persp' plots
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.points=list(
         ## default graphics parameters for 'points'
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.contour=list(
         ## default graphics parameters for 'contour'
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.fv=list(
         ## default graphics parameters for 'plot.fv'
         default=list(),
         check=is.list,
         valid="a list"
         ),
       print.ppm.SE=list(
         ## under what conditions to print estimated SE in print.ppm
         default="poisson",
         check=function(x) { is.character(x) && length(x) == 1 &&
                             x %in% c("always", "poisson", "never") },
         valid="one of the strings \'always\', \'poisson\' or \'never\'"
       ),
       progress = list(
         ## how to display progress reports
         default="tty",
         check=function(x){ x %in% c("tty", "txtbar") },
         valid=paste("one of the strings", dQuote("tty"),
           "or", dQuote("txtbar"))
         ),
       project.fast=list(
         ## whether to cut corners when projecting an invalid ppm object
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       psstG.remove.zeroes=list(
         ## whether to remove zero distances in psstG
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       psstA.ngrid=list(
         ## size of point grid for computing areas in psstA
         default=32,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x >= 8
         },
         valid="a single integer, greater than or equal to 8"
       ),
       psstA.nr=list(
         ## number of 'r' values to consider in psstA
         default=30,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x >= 4
         },
         valid="a single integer, greater than or equal to 4"
       ),
       rmh.p=list(
         ## default value of parameter 'p' in rmh
         default=0.9,
         check=function(x) { is.numeric(x) && length(x) == 1 &&
                             x >= 0 && x <= 1 },
         valid="a single numerical value, between 0 and 1"
       ),
       rmh.q=list(
         ## default value of parameter 'q' in rmh
         default=0.9,
         check=function(x) { is.numeric(x) && length(x) == 1 &&
                             x > 0 && x < 1 },
         valid="a single numerical value, strictly between 0 and 1"
       ),
       rmh.nrep=list(
         ## default value of parameter 'nrep' in rmh
         default=5e5, 
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 0
         },
         valid="a single integer, greater than 0"
       ),
       scalable = list(
         ## whether certain calculations in ppm should be scalable
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       selfcrossing.psp.useCall=list(
         ## whether to use new code in selfcrossing.psp
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       terse = list(
         ## Level of terseness in printed output (higher => more terse)
         default=0,
         check=function(x) { length(x) == 1 && (x %in% 0:3) },
         valid="an integer between 0 and 3"
       ),
       use.Krect=list(
         ## whether to use function Krect in Kest(X) when window is rectangle
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       expand.polynom=list(
         ## whether to expand polynom() in ppm formulae
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       )
       )
# end of options list

reset.spatstat.options <- function() {
  Spatstat.Options <- lapply(.Spat.Stat.Opt.Table,
                               function(z) { z$default })
  putSpatstatVariable("Spatstat.Options", Spatstat.Options)
  invisible(Spatstat.Options)  
}

reset.spatstat.options()

"spatstat.options" <-
function (...) 
{
    Spatstat.Options <- getSpatstatVariable("Spatstat.Options")
    called <- list(...)    

    if(length(called) == 0)
    	return(Spatstat.Options)

    if(is.null(names(called)) && length(called)==1) {
      # spatstat.options(x) 
      x <- called[[1]]
      if(is.null(x))
        return(Spatstat.Options)  # spatstat.options(NULL)
      if(is.list(x))
        called <- x 
    }
    
    if(is.null(names(called))) {
        # spatstat.options("par1", "par2", ...)
	ischar <- unlist(lapply(called, is.character))
	if(all(ischar)) {
          choices <- unlist(called)
          ok <- choices %in% names(Spatstat.Options)
          if(!all(ok))
            stop(paste("Unrecognised option(s):", called[!ok]))
          if(length(called) == 1)
            return(Spatstat.Options[[choices]])
          else
            return(Spatstat.Options[choices])
	} else {
	   wrong <- called[!ischar]
	   offending <- unlist(lapply(wrong,
	   		function(x) { y <- x;
	     		short.deparse(substitute(y)) }))
	   offending <- paste(offending, collapse=",")
           stop(paste("Unrecognised mode of argument(s) [",
		offending,
	   "]: should be character string or name=value pair"))
    	}
    }
# spatstat.options(name=value, name2=value2,...)
    assignto <- names(called)
    if (is.null(assignto) || !all(nzchar(assignto)))
        stop("options must all be identified by name=value")
    recog <- assignto %in% names(.Spat.Stat.Opt.Table)
    if(!all(recog))
	stop(paste("Unrecognised option(s):", assignto[!recog]))
# validate new values
    for(i in seq_along(assignto)) {
      nama <- assignto[i]
      valo <- called[[i]]
      entry <- .Spat.Stat.Opt.Table[[nama]]
      ok <- do.call(entry$check, list(valo))
      if(!ok)
        stop(paste("Parameter", dQuote(nama), "should be",
                   entry$valid))
    }
# reassign
  changed <- Spatstat.Options[assignto]
  Spatstat.Options[assignto] <- called
  putSpatstatVariable("Spatstat.Options", Spatstat.Options)
  
# return 
    invisible(changed)
}

