#
#     options.R
#
#     Spatstat Options
#
#    $Revision: 1.47 $   $Date: 2013/07/17 03:07:17 $
#
#

.spEnv <- new.env()
assign(".Spatstat.Options", list(), envir = .spEnv)

".Spat.Stat.Opt.Table" <-
  list(
       scalable = list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       npixel=list(
         default=128,
         check=function(x){
           is.numeric(x) && (length(x) %in% c(1,2)) && is.finite(x) &&
           all(x == ceiling(x)) && all(x > 1) 
         },
         valid="an integer, or a pair of integers, greater than 1"
        ),
       maxedgewt=list(
         default=100,
         check=function(x){
           is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1
         },
         valid="a finite numeric value, not less than 1"
       ),
       par.binary=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.persp=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.points=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.contour=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.fv=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       image.colfun=list(
         default=function(n){topo.colors(n)},
         check=function(x) {
           is.function(x) && length(formals(x)) > 0 && all(is.character(x(42)))
         },
         valid="a function f(n) that returns character values"
         ),
       ndummy.min=list(
         default=32,
         check=function(x) {
           is.numeric(x) && length(x) <= 2 && all(x == ceiling(x)) && all(x > 1)
         },
         valid="a single integer or a pair of integers, greater than 1"
       ),
       dupC = list(
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       progress = list(
         default="tty",
         check=function(x){ x %in% c("tty", "txtbar") },
         valid=paste("one of the strings", dQuote("tty"),
           "or", dQuote("txtbar"))
         ),
       checkpolygons = list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       checksegments = list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       ngrid.disc=list(
         default=128,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1
         },
         valid="a single integer, greater than 1"
       ),
       gpclib=list(
         default=FALSE,
         check=function(x) {
           message("gpclib is no longer needed")
           return(TRUE)
         },
         valid="a single logical value"
         ),
       maxmatrix=list(
         default=2^24, # 16,777,216
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       huge.npoints=list(
         default=1e6,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       expand=list(
         default=2,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && x > 1
         },
         valid="a single numeric value, greater than 1"
       ),
       fasteval=list(
         default="on",
         check=function(x) { x %in% c("off", "on", "test") },
         valid="one of the strings \'off\', \'on\' or \'test\'"
       ),
       densityC=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       exactdt.checks.data=list(
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       closepairs.newcode=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       crosspairs.newcode=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       psstG.remove.zeroes=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Kcom.remove.zeroes=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       psstA.ngrid=list(
         default=32,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x >= 8
         },
         valid="a single integer, greater than or equal to 8"
       ),
       psstA.nr=list(
         default=30,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x >= 4
         },
         valid="a single integer, greater than or equal to 4"
       ),
       crossing.psp.useCall=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       selfcrossing.psp.useCall=list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       n.bandwidth=list(
         default=32,
         check=function(x) {
           is.numeric(x) && (length(x) == 1) && (x == ceiling(x)) && (x > 2)
         },
         valid="a single integer, greater than 2"
       ),
       project.fast=list(
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       rmh.p=list(
         default=0.9,
         check=function(x) { is.numeric(x) && length(x) == 1 &&
                             x >= 0 && x <= 1 },
         valid="a single numerical value, between 0 and 1"
       ),
       rmh.q=list(
         default=0.9,
         check=function(x) { is.numeric(x) && length(x) == 1 &&
                             x > 0 && x < 1 },
         valid="a single numerical value, strictly between 0 and 1"
       ),
       rmh.nrep=list(
         default=5e5, 
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 0
         },
         valid="a single integer, greater than 0"
       ),
       print.ppm.SE=list(
         default="poisson",
         check=function(x) { is.character(x) && length(x) == 1 &&
                             x %in% c("always", "poisson", "never") },
         valid="one of the strings \'always\', \'poisson\' or \'never\'"
       ),
       nvoxel=list(
         default=2^22,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 2^12
         },
         valid="a single integer, greater than 2^12"
       ),
       fastK.lgcp=list(
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       )
       )
# end of options list

reset.spatstat.options <- function() {
  .Spatstat.Options <- lapply(.Spat.Stat.Opt.Table,
                               function(z) { z$default })
  assign(".Spatstat.Options", .Spatstat.Options, envir = .spEnv)
  invisible(.Spatstat.Options)  
}

reset.spatstat.options()

"spatstat.options" <-
function (...) 
{
    .Spatstat.Options <- get(".Spatstat.Options", envir = .spEnv)
    called <- list(...)    

    if(length(called) == 0)
    	return(.Spatstat.Options)

    if(is.null(names(called)) && length(called)==1) {
      # spatstat.options(x) 
      x <- called[[1]]
      if(is.null(x))
        return(.Spatstat.Options)  # spatstat.options(NULL)
      if(is.list(x))
        called <- x 
    }
    
    if(is.null(names(called))) {
        # spatstat.options("par1", "par2", ...)
	ischar <- unlist(lapply(called, is.character))
	if(all(ischar)) {
          choices <- unlist(called)
          ok <- choices %in% names(.Spatstat.Options)
          if(!all(ok))
            stop(paste("Unrecognised option(s):", called[!ok]))
          if(length(called) == 1)
            return(.Spatstat.Options[[choices]])
          else
            return(.Spatstat.Options[choices])
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
  changed <- .Spatstat.Options[assignto]
  .Spatstat.Options[assignto] <- called
  assign(".Spatstat.Options", .Spatstat.Options, envir = .spEnv)
  
# return 
    invisible(changed)
}

