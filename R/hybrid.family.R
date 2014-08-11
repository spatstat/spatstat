#
#   hybrid.family.R
#
#    $Revision: 1.5 $	$Date: 2011/11/18 04:24:40 $
#
#    Hybrid interactions
#
#    hybrid.family:      object of class 'isf' defining pairwise interaction
#	
# -------------------------------------------------------------------
#	

hybrid.family <-
  list(
       name  = "hybrid",
       print = function(self) {
         cat("Hybrid interaction family\n")
       },
       plot = function(fint, ..., d=NULL, plotit=TRUE, separate=FALSE) {
         # plot hybrid interaction if possible
         verifyclass(fint, "fii")
         inter <- fint$interaction
         if(is.null(inter) || is.null(inter$family)
            || inter$family$name != "hybrid")
           stop("Tried to plot the wrong kind of interaction")
         if(is.null(d)) {
           # compute reach and determine max distance for plots
           dmax <- 1.25 * reach(inter)
           if(!is.finite(dmax)) {
             # interaction has infinite reach
             # Are plot limits specified?
             xlim <- resolve.defaults(list(...), list(xlim=c(0, Inf)))
             if(all(is.finite(xlim))) dmax <- max(xlim) else 
             stop("Interaction has infinite reach; need to specify xlim or d")
           }
           d <- seq(0, dmax, length=256)
         }
         # get fitted coefficients of interaction terms
         # and set coefficients of offset terms to 1         
         Vnames <- fint$Vnames
         IsOffset <- fint$IsOffset
         coeff <- rep(1, length(Vnames))
         names(coeff) <- Vnames
         coeff[!IsOffset] <- fint$coefs[Vnames[!IsOffset]]         
         # extract the component interactions 
         interlist <- inter$par
         # check that they are all pairwise interactions
         families <- unlist(lapply(interlist, function(x) { x$family$name }))
         if(!separate && !all(families == "pairwise")) {
           warning(paste("Cannot compute the resultant function;",
                         "not all components are pairwise interactions;",
                         "plotting each component separately"))
           separate <- TRUE
         }
         # deal with each interaction
         ninter <- length(interlist)
         results <- list()
         for(i in 1:ninter) {
           interI <- interlist[[i]]
           nameI  <- names(interlist)[[i]]
           nameI. <- paste(nameI, ".", sep="")
           # find coefficients with prefix that exactly matches nameI.
           prefixlength <- nchar(nameI.)
           Vprefix <- substr(Vnames, 1, prefixlength)
           relevant <- (Vprefix == nameI.)
           # construct fii object for this component
           fitinI <- fii(interI,
                         coeff[relevant], Vnames[relevant], IsOffset[relevant])
           # convert to fv object
           a <- plot(fitinI, ..., d=d, plotit=FALSE)
           aa <- list(a)
           names(aa) <- nameI
           results <- append(results, aa)
         }
         # computation of resultant is only implemented for fv objects
         if(!separate && !all(unlist(lapply(results, is.fv)))) {
           warning(paste("Cannot compute the resultant function;",
                         "not all interaction components yielded an fv object;",
                         "plotting separate results for each component"))
           separate <- TRUE
         }
         # return separate 'fv' or 'fasp' objects if required
         results <- as.listof(results)
         if(separate) {
           if(plotit) {
             main0 <- "Pairwise interaction components"
             do.call("plot", resolve.defaults(list(results),
                                              list(...),
                                              list(main=main0)))
           }
           return(invisible(results))
         }
         # multiply together to obtain resultant pairwise interaction
         ans <- results[[1]]
         if(ninter >= 2) {
           for(i in 2:ninter) {
             Fi <- results[[i]]
             ans <- eval.fv(ans * Fi)
           }
           copyover <- c("ylab", "yexp", "labl", "desc", "fname")
           attributes(ans)[copyover] <- attributes(results[[1]])[copyover]
         }
         main0 <- "Resultant pairwise interaction"
         if(plotit)
           do.call("plot", resolve.defaults(list(ans),
                                            list(...),
                                            list(main=main0)))
         return(invisible(ans))
       },
       eval  = function(X,U,EqualPairs,pot,pars,correction, ...) {
         # `pot' is ignored; `pars' is a list of interactions
         nU <- length(U$x)
         V <- matrix(, nU, 0)
         IsOffset <- logical(0)
         for(i in 1:length(pars)) {
           # extract i-th component interaction
           interI <- pars[[i]]
           nameI  <- names(pars)[[i]]
           # compute potential for i-th component
           VI <- evalInteraction(X, U, EqualPairs, interI, correction, ...)
           if(ncol(VI) > 0) {
             if(ncol(VI) > 1 && is.null(colnames(VI))) # make up names
               colnames(VI) <- paste("Interaction", seq(ncol(VI)), sep=".")
             # prefix label with name of i-th component 
             colnames(VI) <- paste(nameI, dimnames(VI)[[2]], sep=".")
             # handle IsOffset
             offI <- attr(VI, "IsOffset")
             if(is.null(offI))
               offI <- rep(FALSE, ncol(VI))
             # tack on
             IsOffset <- c(IsOffset, offI)
             # append to matrix V
             V <- cbind(V, VI)
           }
         }
         if(any(IsOffset))
           attr(V, "IsOffset") <- IsOffset
         return(V)
       },
       suffstat = NULL
)

class(hybrid.family) <- "isf"


