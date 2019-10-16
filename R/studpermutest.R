#'
#'  studpermtest.R
#' 
#'  Original by Ute Hahn 2014
#'
#' $Revision: 1.10 $ $Date: 2019/10/16 02:36:54 $
#' 
#' Studentized permutation test for comparison of grouped point patterns;
#' functions to generate these grouped point patterns;
#' wrapper for test of reweighted second order stationarity.
#' 


#' studpermu.test
#' studentized permutation test for grouped point patterns
#' interpreted version, random permutations only. 
#' A group needs to contain at least two point patterns with at least minpoints each.
#
#' X               the data, may be a list of lists of point patterns, or a hyperframe
#' formula         if X is a hyperframe, relates point patterns to factor variables that
#'                 determine the groups. If missing, the first column of X that contains
#'                 a factor variable is used.
#' summaryfunction the function used in the test
#' ...             additional arguments for summaryfunction
#' rinterval       r-interval where summaryfunction is evaluated. If NULL, the
#'                 interval is calculated from spatstat defaults 
#'                 (intersection for all patterns)    
#' nperm           number of random permutations
#' use.Tbar        use the alternative test statistic, for summary functions with
#'                 roughly constant variance, such as K/r or L
#' minpoints       the minimum number of points a pattern needs to have. Patterns
#'                 with fewer points are not used.
#' rsteps          discretization steps of the r-interval
#' r               arguments at which to evaluate summaryfunction, overrides rinterval
#'                 Should normally not be given, replace by rinterval instead,
#'                 this allows r_0 > 0. Also, there is no plausibility check for r so far
#' arguments.in.data if TRUE, individual extra arguments to summary function that 
#'                 change are taken from X (which has to be a hyperframe then).
#'                 Assumes that the first argument of summaryfunction always is the
#'                 point pattern.
#'                 This is meant for internal purposes (automatisation)
#
#' returns an object of classes htest and studpermutest, that can be plotted. The
#' plot shows the summary functions for the groups (and the means if requested)

studpermu.test <- local({

studpermu.test <-
  function (X, formula, summaryfunction = Kest,
            ...,
            rinterval = NULL, nperm = 999,
            use.Tbar = FALSE, # the alternative statistic, use with K/r or L  
            minpoints = 20, 
            rsteps = 128, r = NULL,
            arguments.in.data = FALSE) {
    #' ---- the loooong preliminaries -------  

    #' ---- argument checking paranoia ----
    if (arguments.in.data & !is.hyperframe(X))
      stop(paste("X needs to be a hyperframe",
                 "if arguments for summary function are to be retrieved"),
           call.=FALSE)
    stopifnot(is.function(summaryfunction))
    #' there could be more...
  
    #' first prepare the data
    if(is.hyperframe(X)) {
      if(dim(X)[2] < 2) 
        stop(paste("Hyperframe X needs to contain at least 2 columns,",
                   "one for patterns, one indicating groups"),
             call.=FALSE)
      data <- X # renaming for later. 
      Xclass <- unclass(X)$vclass
      factorcandidate <-
        Xclass %in% c("integer", "numeric", "character", "factor")
      ppcandidate <- Xclass == "ppp"
      names(factorcandidate) <- names(ppcandidate) <-
        names(Xclass) <- Xnames <- names(X)
      if(all(!factorcandidate) || all(!ppcandidate))
        stop(paste("Hyperframe X needs to contain at least a column",
                   "with point patterns, and one indicating groups"),
             call.=FALSE)
      if(!missing(formula)){
        #' safety precautions ;-)
        if(!inherits(formula, "formula")) 
          stop(paste("Argument", dQuote("formula"), "should be a formula"))
        if (length(formula) < 3) 
          stop(paste("Argument", sQuote("formula"),
                     "must have a left hand side"))
        rhs <- rhs.of.formula(formula)
        ppname <- formula[[2]]
        if (!is.name(ppname)) 
          stop("Left hand side of formula should be a single name")
        ppname <- paste(ppname)
        if(!ppcandidate[ppname])
          stop(paste("Left hand side of formula",
                     "should be the name of a column of point patterns"),
               call.=FALSE)
        groupvars <- all.vars(as.expression(rhs))
        if(!all(groupvars %in% Xnames) || any(!factorcandidate[groupvars]))
          stop(paste("Not all variables on right hand side of formula",
                     "can be interpreted as factors"),
               call.=FALSE)
        #' make the groups to be compared
        group <-
          interaction(lapply(as.data.frame(data[ , groupvars, drop=FALSE]),
                             factor))
        #' rename the point patterns, needs the patch      
        newnames <- Xnames
        newnames[Xnames == ppname] <- "pp"
        names(data) <- newnames
        data$group <- group
      } else {
        #' No formula supplied.
        #' Choose first ppp column and first factor column to make pp and groups
        thepp <- which.max(ppcandidate)
        thegroup <- which.max(factorcandidate)
        #' fake formula for output of test result
        formula <- as.formula(paste( Xnames[thepp],"~", Xnames[thegroup]))
        newnames <- Xnames
        newnames[thepp] <- "pp"
        newnames[thegroup] <- "group"
        names(data) <- newnames
        data$group <- as.factor(data$group)
      }
    } else {
      #' X is not a hyperframe, but hopefully a list of ppp
      if(!is.list(X))
        stop("X should be a hyperframe or a list of lists of point patterns")
      if (!is.list(X[[1]]) || !is.ppp(X[[1]][[1]]))
        stop("X is a list, but not a list of lists of point patterns")
      nams <- names(X)
      if(is.null(nams)) nams <- paste("group", seq_along(X))
    
      pp <- list()
      group <- NULL
      for (i in seq_along(X)) {
        pp <- c(pp, X[[i]])
        group <- c(group, rep(nams[i], length(X[[i]])))
      }
      group <- as.factor(group) 
      data <-  hyperframe(pp = pp, group = group)
      ppname <- "pp"
    }
  
    framename <- short.deparse(substitute(X))
    fooname <- short.deparse(substitute(summaryfunction))
  
    #' sorting out the patterns that contain too few points
  
    OK <- sapply(data$pp, npoints) >= minpoints
    if((nbad <- sum(!OK)) > 0)
      warning(paste(nbad,
                    "patterns have been discarded",
                    "because they contained fewer than",
                    minpoints, "points"),
              call.=FALSE)
    data <- data[OK, ,drop=FALSE]
    pp <- data$pp

    #' ---- the groups,
    #' or what remains after discarding the poor patterns with few points -----
  
    #' check if at least two observations in each group
    groupi <- as.integer(data$group)
    ngroups <- max(groupi)
    if(ngroups < 2)
      stop(paste("Sorry, after discarding patterns with fewer than",
                 minpoints,
                 "points,",
                 if(ngroups < 1) "nothing" else "only one group",
                 "is left over.",
                 "\n- nothing to compare, take a break!"),
           call.=FALSE)
 
    lev <- 1:ngroups
    m <- as.vector(table(groupi))
    if (any(m < 3))
      stop(paste("Data groups need to contain at least two patterns;",
                 "\nafter discarding those with fewer than",
                 minpoints,
                 "points, the remaining group sizes are",
                 commasep(m)),
           call.=FALSE)
    #' check if number of possible outcomes is small
    #' WAS: npossible <- factorial(sum(m))/prod(factorial(m))/prod(factorial(table(m)))
    lognpossible <- lgamma(sum(m)+1)-sum(lgamma(m+1))-sum(lgamma(table(m)+1))
    if (lognpossible < log(max(100, nperm)))
      warning("Don't expect exact results - group sizes are too small")
  
    #' --------- real preliminaries now ------
    
    #' get interval for arguments

    if(!is.null(r)){
      rinterval <- range(r)
      rsteps <- length(r)
    } else if (is.null(rinterval)) {
      foochar <- substr(fooname, 1, 1)
      if (foochar %in% c("p", "L")) foochar <- "K"
      if (fooname %in%  c("Kscaled", "Lscaled")) foochar <- "Kscaled"
      rinterval <-
        c(0, min(with(data, rmax.rule(foochar, Window(pp), intensity(pp)))))
    }  
    
    ranger <- diff(range(rinterval))
    #' r sequence needs to start at 0 for Kest and such
    rr <- r %orifnull% seq(0, rinterval[2], length.out = rsteps + 1) 
    taker <- rr >= rinterval[1] & rr <= rinterval[2] # used for testing
 
    #' now estimate the summary function, finally...
    #' TO DO!!!! Match function call of summary function with data columns!
    #' use arguments.in.data, if applicable. This is for inhomogeneous summary 
    #' functions

    #' Force all calls to summaryfunction to use the same edge correction,
    #' rather than allowing correction to depend on npoints
    needcorx <- "correction" %in% names(formals(summaryfunction))
    gavecorx <- "correction" %in% names(list(...))
    corx <- if(needcorx && !gavecorx) "best" else NULL

    #' --------- retrieve arguments for summary function from data, hvis det er
    fvlist <-
      if(arguments.in.data) {
        #' use arguments in hyperframe 'data' as well as explicit arguments
        if(is.null(corx)) {
          multicall(summaryfunction, pp, data, r = rr, ...)
        } else {
          multicall(summaryfunction, pp, data, r = rr, ..., correction=corx)
        }
      } else {
        #' use explicit arguments only
        if(is.null(corx)) {
          with(data, summaryfunction(pp, r = rr, ...))
        } else {
          with(data, summaryfunction(pp, r = rr, ..., correction=corx))
        }
      }

    fvtemplate <- fvlist[[1]]
 
    valu <- attr(fvtemplate, "valu")
    argu <- attr(fvtemplate, "argu")

    foar <- sapply(lapply(fvlist, "[[", valu),
                   "[", taker)

    #' --------- the real stuff --------------
  
    #' function that calculates the discrepancy
    #' slow version
    combs <- combn(lev, 2)

    #' --------- now do the real real stuff :-)  -------------
  
    #' generate "simulated values" from random permutations. 
    #' possible improvement for future work:
    #' If the number of all permutations (combis) is small,
    #' first generate all permutations and then
    #' sample from them to improve precision
  
    predigested <- list(lev=lev,
                        foar=foar, m=m, combs=combs, rrr=rr[taker],
                        ranger=ranger)
    if(use.Tbar) {
      Tobs <- Tbarstat(groupi, predigested)
      Tsim <- replicate(nperm, Tbarstat(sample(groupi), predigested))    
    } else {
      Tobs <- Tstat(groupi, predigested)
      Tsim <- replicate(nperm, Tstat(sample(groupi), predigested))
    }  
    names(Tobs) <- if(use.Tbar) "Tbar" else "T"

    pval <- (1 + sum(Tobs < Tsim))/(1 + nperm)
  
    #' ----- making a test object -----
    method <- c("Studentized permutation test for grouped point patterns",
                if(is.hyperframe(X)) pasteFormula(formula) else NULL,
                choptext(ngroups, "groups:",
                         paste(levels(data$group), collapse=", ")),
                choptext("summary function:",
                         paste0(fooname, ","), 
                         "evaluated on r in",
                         prange(rinterval)),
                choptext("test statistic:",
                         if(use.Tbar) "Tbar," else "T,", 
                         nperm, "random permutations"))
    fooshort <- switch(fooname, pcf = "pair correlation ",
                       Kinhom = "inhomogeneous K-",
                       Linhom = "inhomogeneous L-",
                       Kscaled = "locally scaled K-",
                       Lscaled = "locally scaled L-",
                       paste(substr(fooname, 1, 1),"-",sep=""))
    alternative <- c(paste("not the same ",fooshort,"function", sep=""))
  
    testerg <- list(statistic = Tobs, 
                    p.value = pval, 
                    alternative = alternative,            
                    method = method, 
                    data.name = framename)
    class(testerg) <- c("studpermutest", "htest")
    #' Add things for plotting
    
    #' prepare the fvlist, so that it only contains the estimates used,
    fvs <- lapply(fvlist, "[.fv", j=c(argu, valu))
    #' with rinterval as alim
    fvs <- lapply(fvs, "attr<-", which="alim", value=rinterval)

    testerg$curves <- hyperframe(fvs = fvs, groups = data$group)
  
    fvtheo <- fvlist[[1]]
    fvnames(fvtheo, ".y") <- "theo"
    attr(fvtheo, "alim") <- rinterval
    testerg$curvtheo <- fvtheo[ , c(argu, "theo")]
  
    #' group means
    grmn <- lapply(lev, splitmean, ind=groupi, f=foar)
    testerg$groupmeans <- lapply(grmn, makefv, xvals=rr[taker],
                                 template=fvtheo)
  
    return(testerg)
  }  

  splitmean <- function(l, ind, f) {
    apply(f[ , ind == l], 1, mean)
  }
  splitvarn <- function(l, ind, f, m) {
    apply(f[ , ind == l], 1, var) / m[l]
  }
  studentstat <- function(i, grmean, grvar) {
    (grmean[, i[1]] - grmean[, i[2]])^2 / (grvar[i[1],] + grvar[i[2], ])
  }
    
  Tstat <- function (ind = groupi, predigested) {
    #' predigested should be a list with entries lev, foar, m, combs, rrr
    with(predigested, {
      grmean <- sapply(lev, splitmean, ind=ind, f=foar)
      grvar <- t(sapply(lev, splitvarn, ind=ind, f=foar, m=m))
      y <- apply(combs, 2, studentstat, grmean=grmean, grvar=grvar)
      sum(apply(y, 2, trapint, x = rrr))
    })
  } 

  intstudent <- function(i, rrr, grmean, meangrvar) {
    trapint(rrr, (grmean[, i[1]] - grmean[, i[2]])^2 / 
            (meangrvar[i[1]] + meangrvar[i[2]]))
  }
    
  Tbarstat <- function (ind = groupi, predigested) {
    #' predigested should be a list
    #' with entries lev, foar, m, combs, rrr, ranger
    with(predigested, {
      grmean <- sapply(lev, splitmean, ind=ind, f=foar)
      grvar <- t(sapply(lev, splitvarn, ind=ind, f=foar, m=m))
      meangrvar <- apply(grvar, 1, trapint, x = rrr)/ranger
      sum(apply(combs, 2, intstudent, 
                rrr=rrr, grmean=grmean, meangrvar=meangrvar))
      #' trapint(rr[taker], grvar[i[1],] + grvar[i[2], ]))))
    })
  }

  makefv <- function(yvals, xvals, template) {
    fdf <- data.frame(r = xvals, y = yvals)
    argu <- fvnames(template, ".x")
    valu <- fvnames(template, ".y")
    names(fdf) <- c(argu,valu)
    fv(fdf, argu = argu, ylab = attr(template, "ylab"), valu = valu, 
       fmla = attr(template,"fmla"), alim = attr(template, "alim"))
  }

  #' Trapezoidal rule approximation to integral
  #' ------- Trapezregel, mit Behandlung von NAns:
  #'                  die werden einfach ignoriert ----
  trapint <- function(x, y) {
    nonan <- !is.na(y)
    nn <- sum(nonan)
    if(nn < 2L) return(0)
    Y <- y[nonan]
    X <- x[nonan]
    0.5 * sum( (Y[-1] + Y[-nn]) * diff(X))
  }

  #' call foo(x, further arguments) repeatedly
  #' further arguments are taken from hyperframe H and ...
  multicall <- function(foo, x, H, ...){
    stopifnot(is.hyperframe(H))
    if (is.hyperframe(x)) {
      x <- as.list(x)[[1]] 
    } else if(!is.list(x))
      stop("in multicall: x should be a hyperframe or list", call.=FALSE)
  
    #' check if same length
    nrows <- dim(H)[1]
    if (length(x) != nrows)
      stop(paste("in multicall: x and H need to have",
                 "the same number of rows or list elements"),
           call.=FALSE)
    dotargs <- list(...)
    hnames <- names(H)
    argnames <- names(formals(foo))#' always assume first argument is given
  
    ppname <- argnames[1]
    argnames <- argnames[-1]
    dotmatch <- pmatch(names(dotargs), argnames)
    dotmatched <- dotmatch[!is.na(dotmatch)]
    dotuseargs <- dotargs[!is.na(dotmatch)]
    restargs <- if(length(dotmatched) >0) argnames[-dotmatched] else argnames
    hmatch <- pmatch(hnames, restargs)
    huse <- !is.na(hmatch)
    lapply(seq_len(nrows), function (i) 
           do.call(foo, c(list(x[[i]]),
                          as.list(H[i, huse, drop=TRUE, strip=FALSE]),
                          dotargs))) 
  }

  studpermu.test
})


#' ------------------- plot studpermutest ---------------------------------------
#
#' plot.studpermutest
#' plot the functions that were used in studperm.test
#' also plot group means, if requested
#
#' x               a studpermtest object, the test result
#' fmla            a plot formula as in plot.fv, should be generic, using "." for values
#' ...             further plot parameters
#' col, lty, lwd   parameter (vectors) for plotting the individual summary functions,
#'                 according to group, if vectors  
#' col.theo, lty.theo, lwd.theo if not all are NULL, the "theo" curve is also plotted
#' lwd.mean        a multiplyer for the line width of the group means. 
#'                 if NULL, group means are not plotted, defaults to NULL
#' lty.mean, col.mean  selbsterklaerend
#' separately      generate a separate plot for each group (then no legends are plotted)
#' meanonly        do not plot individual summary functions
#' legend          if TRUE, and plots are not separate, plot a legend 
#' legendpos       ...
#' lbox            if TRUE, draw box around legend. Defaults to FALSE
#' add             ...

plot.studpermutest <- local({

  plot.studpermutest <-
    function(x, fmla, ..., 
             lty = NULL, col = NULL, lwd = NULL,
             lty.theo = NULL, col.theo = NULL, lwd.theo = NULL,
             lwd.mean = if(meanonly) 1 else NULL,
             lty.mean = lty, col.mean = col,
             separately = FALSE, meanonly = FALSE,
             main = if(meanonly) "group means" else NULL,
             xlim = NULL, ylim = NULL, ylab = NULL,
             legend = !add, legendpos = "topleft", lbox=FALSE,
             add = FALSE) {
      stopifnot(inherits(x, "studpermutest")) 
      env.user <- parent.frame()
      curvlists <- split(x$curves, x$curves$groups)
      ngroups <- length(curvlists)
      gnames <- names(curvlists)
      #' check if theoretical functions shall be plottet
      plottheo <- !(is.null(lty.theo) & is.null(col.theo) & is.null(lwd.theo))
      #' prepare plot parameters for groups
      if (is.null(lty)) lty <- 1:ngroups
      if (is.null(col)) col <- 1:ngroups
      if (is.null(lwd)) lwd <- par("lwd")
      if (is.null(col.mean)) col.mean <- col
      if (is.null(lty.mean)) lty.mean <- lty
      lty <- rep(lty, length.out = ngroups)
      col <- rep(col, length.out = ngroups)
      lwd <- rep(lwd, length.out = ngroups)
      col.mean <- rep(col.mean, length.out = ngroups)
      lty.mean <- rep(lty.mean, length.out = ngroups)
      if (plottheo){
        if (is.null(lty.theo)) lty.theo <- ngroups + 1#par("lty")
        if (is.null(col.theo)) col.theo <- ngroups + 1 #par("col")
        if (is.null(lwd.theo)) lwd.theo <- par("lwd")    
      }
      #' transporting the formula in ... unfortunately does not work
      #' for the axis labels, because the fvs contain only one variable.
      #' Have to knit them self
  
      if (is.null(ylab)) {
        if (!missing(fmla)) {
          #' puha. det bliver noget lappevaerk.
          fmla <- as.formula(fmla, env=env.user)
          map <- fvlabelmap(x$curvtheo) 
          lhs <- lhs.of.formula(as.formula(fmla))
          ylab <- eval(substitute(substitute(le, mp),
                                  list(le = lhs, mp = map)))   
        } else ylab <- attr(x$curvtheo, "yexp")
      } 
      if (missing(fmla)) fmla <- attr(x$curvtheo, "fmla")
      if(!is.null(lwd.mean)) lwd.Mean <- lwd.mean*lwd
      if(separately) {
        for (i in seq_along(gnames)) {
          if(!meanonly) 
            plot.fvlist(curvlists[[i]]$fvs, fmla, ..., 
                        col = col[i], lwd = lwd[i], lty= lty[i],
                        xlim = xlim, ylim = ylim, ylab = ylab,
                        main = gnames[i])
          if (!is.null(lwd.mean))
            plot(x$groupmeans[[i]], fmla, ...,
                 col = col.mean[i], lwd = lwd.Mean[i], lty = lty.mean[i], 
                 main = gnames[i], add = !meanonly, ylim = ylim)
          if (plottheo)
            plot(x$curvtheo, fmla, ..., add = TRUE, 
                 col = col.theo, lwd = lwd.theo, lty = lty.theo)
        }
      } else {  
        #' ---- TODO SIMPLIFY! they should all have the same x-range,
        #' just check y-range ----
        lims <- if (meanonly) {
          plot.fvlist(x$groupmeans, fmla,..., limitsonly = TRUE)
        } else {
          as.data.frame(apply(sapply(curvlists, 
            function(C) plot.fvlist(C$fvs, fmla,..., limitsonly = TRUE)),
                              1, range))
        }
        if(is.null(xlim)) xlim <- lims$xlim
        if(is.null(ylim)) ylim <- lims$ylim
        iadd <- add
        for (i in seq_along(gnames)) {
          if(!meanonly)
            plot.fvlist(curvlists[[i]]$fvs, fmla, ..., 
                        col = col[i], lwd = lwd[i], lty= lty[i],
                        xlim = xlim, ylim = ylim, ylab= ylab,
                        main = main,
                        add = iadd)
          iadd <- iadd | !meanonly
          if (!is.null(lwd.mean))
            plot(x$groupmeans[[i]], fmla, ...,
                 col = col.mean[i], lwd = lwd.Mean[i], lty = lty.mean[i],
                 add = iadd,
                 xlim = xlim, ylim = ylim, ylab= ylab, main=main)  
        if (plottheo)
          plot(x$curvtheo, fmla, ..., add = TRUE, 
               col = col.theo, lwd = lwd.theo, lty = lty.theo,
               xlim = xlim, ylim = ylim, ylab= ylab, main=main)  
          iadd <- TRUE
        } 
        if(legend) {
          if(meanonly) {
            lwd <- lwd.Mean
            col <- col.mean
            lty <- lty.mean
          }
          if(plottheo){
            gnames <- c(gnames, "Poisson mean")
            col <- c(col, col.theo)
            lty <- c(lty, lty.theo)
            lwd <- c(lwd, lwd.theo)
          }
          legend(legendpos, gnames, col = col, lty = lty, lwd = lwd, 
                 bty=ifelse(lbox, "o", "n"))
        }  
      }
      return(invisible(NULL))
    }

  #' ------------------ Helper function----------------
  #' flist: list of fv, with plot method

  plot.fvlist <- function(x, fmla, ..., xlim=NULL, ylim=NULL, 
                          add = FALSE, limitsonly = FALSE, main=NULL){
    #' no safety precautions
    if (missing(fmla)) fmla <- attr(x[[1]], "fmla")
    if (!add | limitsonly) {
      lims <- sapply(x, plot, fmla, ..., limitsonly = TRUE)
      if(is.null(xlim)) xlim = range(unlist(lims[1,]))
      if(is.null(ylim)) ylim = range(unlist(lims[2,]))
      lims=list(xlim=xlim, ylim=ylim)
      if(limitsonly) return(lims) 
      plot(x[[1]], fmla, ..., xlim = xlim, ylim = ylim, main = main)
    }
    else plot(x[[1]], fmla,..., add=T)
    for (foo in x[-1]) plot(foo, fmla, ..., add=T)
  }
 
  plot.studpermutest
})

