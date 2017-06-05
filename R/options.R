#
#     options.R
#
#     Spatstat options and other internal states
#
#    $Revision: 1.80 $   $Date: 2017/06/05 10:31:58 $
#
#

.spEnv <- new.env()

putSpatstatVariable <- function(name, value) {
  assign(name, value, envir=.spEnv)
}
getSpatstatVariable <- function(name) {
  get(name, envir=.spEnv)
}
existsSpatstatVariable <- function(name) {
  exists(name, envir=.spEnv)
}

putSpatstatVariable("Spatstat.Options", list())
putSpatstatVariable("Spatstat.ProgressBar", NULL)
putSpatstatVariable("Spatstat.ProgressData", NULL)
putSpatstatVariable("warnedkeys", character(0))

## Kovesi's uniform colour map, row 29, linear 'bmy'
putSpatstatVariable("DefaultImageColours", 
c("#000C7D", "#000D7E", "#000D80", "#000E81", "#000E83", "#000E85", 
"#000F86", "#000F88", "#00108A", "#00108B", "#00118D", "#00118F", 
"#001190", "#001292", "#001293", "#001295", "#001396", "#001398", 
"#001399", "#00149A", "#00149C", "#00149D", "#00149E", "#00159F", 
"#0015A0", "#0015A1", "#0015A2", "#0015A3", "#0015A4", "#0016A5", 
"#0016A6", "#0016A6", "#0016A7", "#0016A8", "#0016A8", "#0016A8", 
"#0A16A9", "#1516A9", "#1D15A9", "#2315A9", "#2915A9", "#2F15A8", 
"#3414A8", "#3914A7", "#3E13A6", "#4313A5", "#4712A4", "#4C12A3", 
"#5011A2", "#5311A1", "#5710A0", "#5A0F9F", "#5E0F9E", "#610E9E", 
"#640E9D", "#670D9C", "#6A0D9B", "#6C0C9A", "#6F0B99", "#720B98", 
"#740A98", "#770A97", "#790996", "#7C0896", "#7E0895", "#800794", 
"#810794", "#840693", "#860692", "#880692", "#8A0591", "#8C0591", 
"#8E0490", "#900490", "#92048F", "#94038F", "#96038E", "#98038E", 
"#9A028D", "#9C028D", "#9E028D", "#A0018C", "#A2018C", "#A4018B", 
"#A6018B", "#A8008A", "#AA008A", "#AB0089", "#AD0089", "#AF0088", 
"#B10088", "#B30087", "#B50087", "#B70086", "#B80086", "#BA0086", 
"#BC0085", "#BE0085", "#C00084", "#C20084", "#C30083", "#C50083", 
"#C70082", "#C90082", "#CB0081", "#CD0081", "#CE0080", "#D00080", 
"#D20080", "#D40080", "#D5007F", "#D7007F", "#D9007E", "#DA007E", 
"#DC007D", "#DD007C", "#DF017C", "#E1027B", "#E2047B", "#E4067A", 
"#E5087A", "#E70B79", "#E80D78", "#E91078", "#EB1277", "#EC1477", 
"#ED1676", "#EF1875", "#F01A75", "#F11C74", "#F31E73", "#F42073", 
"#F52272", "#F62471", "#F72671", "#F82870", "#FA2A6F", "#FB2C6F", 
"#FC2E6E", "#FD306D", "#FE326C", "#FE346C", "#FE366B", "#FE386A", 
"#FE3A6A", "#FE3D69", "#FE3F68", "#FE4167", "#FE4366", "#FE4566", 
"#FE4765", "#FE4964", "#FE4B63", "#FE4D62", "#FE5062", "#FE5261", 
"#FE5460", "#FE565F", "#FE585E", "#FE5A5D", "#FE5D5C", "#FE5F5B", 
"#FE615B", "#FE635A", "#FE6559", "#FE6758", "#FE6A57", "#FE6C56", 
"#FE6E55", "#FE7054", "#FE7253", "#FE7452", "#FE7651", "#FE7850", 
"#FE7A4E", "#FE7C4D", "#FE7E4C", "#FE7F4B", "#FE804A", "#FE8249", 
"#FE8448", "#FE8647", "#FE8745", "#FE8944", "#FE8B43", "#FE8D42", 
"#FE8E40", "#FE903F", "#FE923E", "#FE943C", "#FE953B", "#FE9739", 
"#FE9938", "#FE9A36", "#FE9C35", "#FE9E33", "#FE9F32", "#FEA130", 
"#FEA22F", "#FEA42E", "#FEA52C", "#FEA72B", "#FEA82A", "#FEAA29", 
"#FEAB28", "#FEAD27", "#FEAE26", "#FEB026", "#FEB125", "#FEB324", 
"#FEB423", "#FEB523", "#FEB722", "#FEB822", "#FEBA21", "#FEBB20", 
"#FEBC20", "#FEBE1F", "#FEBF1F", "#FEC11F", "#FEC21E", "#FEC31E", 
"#FEC51E", "#FEC61D", "#FEC71D", "#FEC91D", "#FECA1D", "#FECB1D", 
"#FECD1D", "#FECE1C", "#FECF1C", "#FED11C", "#FED21C", "#FED31C", 
"#FED51C", "#FED61D", "#FED71D", "#FED91D", "#FEDA1D", "#FEDB1D", 
"#FEDD1D", "#FEDE1E", "#FEDF1E", "#FEE11E", "#FEE21E", "#FEE31F", 
"#FEE51F", "#FEE61F", "#FEE720", "#FEE820", "#FEEA21", "#FEEB21", 
"#FEEC22", "#FEEE22", "#FEEF23", "#FEF023"))

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
       checkpolygons = list(
         ## superseded
         superseded=TRUE,
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
       densityTransform=list(
         ## use experimental new C routines for 'density.ppp'
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
       dpp.maxmatrix=list(
         ## maximum size of matrix in dppeigen
         default=2^24, # 16,777,216
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
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
       expand.polynom=list(
         ## whether to expand polynom() in ppm formulae
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       fasteval=list(
         ## whether to use 'fasteval' code if available
         default="on",
         check=function(x) { x %in% c("off", "on", "test") },
         valid="one of the strings \'off\', \'on\' or \'test\'"
       ),
       fastpois=list(
         # whether to use fast algorithm for rpoispp() when lambda is an image
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       fastthin=list(
         # whether to use fast C algorithm for rthin() when P is constant
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
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
         superseded=TRUE, 
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
#         default=function(n){topo.colors(n)},
         default=function(n) {
           z <- getSpatstatVariable("DefaultImageColours")
           interp.colours(z, n)
         },
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
       par.pp3=list(
         ## default graphics parameters for 'plot.pp3'
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
         check=function(x){ x %in% c("tty", "tk", "txtbar") },
         valid="one of the strings 'tty', 'tk' or 'txtbar'"
         ),
       project.fast=list(
         ## whether to cut corners when projecting an invalid ppm object
         default=FALSE,
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
       psstG.remove.zeroes=list(
         ## whether to remove zero distances in psstG
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
      eroded.intensity=list(
         ## whether to compute intensity estimate in eroded window
         ## e.g. for Kcom, Gcom
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       rmh.nrep=list(
         ## default value of parameter 'nrep' in rmh
         default=5e5, 
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 0
         },
         valid="a single integer, greater than 0"
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
         check=function(x) { length(x) == 1 && (x %in% 0:4) },
         valid="an integer between 0 and 4"
       ),
       transparent=list(
         ## whether to allow transparent colours in default colour maps
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       units.paren=list(
         default="(",
         check=function(x) {
           is.character(x) && (length(x) == 1) &&
             (x %in% c("(", "[", "{", ""))
         },
         valid="one of the strings '(', '[', '{' or '' "
       ),
       use.Krect=list(
         ## whether to use function Krect in Kest(X) when window is rectangle
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Cwhist=list(
         ## whether to use C code for whist
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Cbdrymask=list(
         ## whether to use C code for bdry.mask
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       kppm.canonical=list(
         ## whether to use 'canonical' parameters in kppm
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       kppm.adjusted=list(
         ## experimental
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       check.rpanel.loaded=list(
         # internal debugging
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       check.RandomFields.loaded=list(
         # this is working OK so no need to check unless debugging
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       check.RandomFieldsUtils.loaded=list(
         # this is working OK so no need to check unless debugging
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Clinequad = list(
         # use C code for 'linequad'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Ccountends = list(
         # use C code for 'countends'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Clinearradius = list(
         # use C code for 'boundingradius.linnet'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Cnndistlpp = list(
         # use C code for 'nndist.lpp'/'nnwhich.lpp'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       Cnncrosslpp = list(
         # use C code for 'nncross.lpp'
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       ),
       developer = list(
         # general purpose; user is a developer; use experimental code, etc
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1 },
         valid="a single logical value"
       )
    )
# end of options list

reset.spatstat.options <- function() {
  Spatstat.Options <- lapply(.Spat.Stat.Opt.Table, getElement, name="default")
  putSpatstatVariable("Spatstat.Options", Spatstat.Options)
  invisible(Spatstat.Options)  
}

reset.spatstat.options()

spatstat.options <- local({

  spatstat.options <- function (...) {
    Spatstat.Options <- getSpatstatVariable("Spatstat.Options")
    called <- list(...)    

    if(length(called) == 0) {
      # return all options, except superseded ones
      allofem <- .Spat.Stat.Opt.Table[names(Spatstat.Options)]
      retain <- sapply(lapply(allofem, getElement, name="superseded"), is.null)
      return(Spatstat.Options[retain])
    }
    
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
	   offending <- sapply(wrong, ShortDeparse)
	   offending <- paste(offending, collapse=",")
           stop(paste("Unrecognised mode of argument(s) [",
		offending,
	   "]: should be character string or name=value pair"))
    	}
    }
    ## spatstat.options(name=value, name2=value2,...)
    assignto <- names(called)
    if (is.null(assignto) || !all(nzchar(assignto)))
        stop("options must all be identified by name=value")
    recog <- assignto %in% names(.Spat.Stat.Opt.Table)
    if(!all(recog))
	stop(paste("Unrecognised option(s):", assignto[!recog]))
    ## validate new values
    for(i in seq_along(assignto)) {
      nama <- assignto[i]
      valo <- called[[i]]
      entry <- .Spat.Stat.Opt.Table[[nama]]
      ok <- entry$check(valo)
      if(!ok)
        stop(paste("Parameter", dQuote(nama), "should be",
                   entry$valid))
    }
    ## reassign
    changed <- Spatstat.Options[assignto]
    Spatstat.Options[assignto] <- called
    putSpatstatVariable("Spatstat.Options", Spatstat.Options)
  
    ## return 
    invisible(changed)
  }

  ShortDeparse <- function(x) {
    y <- x
    dont.complain.about(y)
    short.deparse(substitute(y))
  }
    
  spatstat.options
})

