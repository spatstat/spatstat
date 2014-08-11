#
# rmhsnoop.R
#
#   visual debug mechanism for rmh
#
#   $Revision: 1.24 $  $Date: 2013/04/25 06:37:43 $
#
#   When rmh is called in visual debug mode (snooping = TRUE),
#   it calls e <- rmhSnoopEnv(...) to create an R environment 'e'
#   containing variables that will represent the current state
#   of the M-H algorithm with initial state X and model reach R.
#
#   The environment 'e' is passed to the C routine xmethas.
#   This makes it possible for data to be exchanged between
#   the C and R code.
#
#   When xmethas reaches the debugger's stopping time,
#   the current state of the simulation and the proposal
#   are copied from C into the R environment 'e'.
#
#   Then to execute the visual display, the C code calls
#   'eval' to execute the R function rmhsnoop().
#
#   The function rmhsnoop uses the 'simplepanel' class
#   to generate a plot showing the state of the simulation
#   and the proposal, and then wait for point-and-click input using
#   locator(). 
#  
#   When rmhsnoop() exits, it returns an integer giving the
#   (user-specified) next stoppping time. This is read back into
#   the C code. Then xmethas resumes simulations.
#
#   I said it was simple! %^]

rmhSnoopEnv <- function(Xinit, Wclip, R) {
  stopifnot(is.ppp(Xinit))
  # Create an environment that will be accessible to R and C code
  e <- new.env()
  # initial state (point pattern)
  X <- Xinit
  assign("Wsim",     as.owin(X),      envir=e)
  assign("xcoords",  coords(X)[,1],   envir=e)
  assign("ycoords",  coords(X)[,2],   envir=e)
  if(is.multitype(X)) {
    mcodes <- as.integer(marks(X)) - 1
    mlevels <- levels(marks(X))
    assign("mcodes",  mcodes,  envir=e)
    assign("mlevels", mlevels, envir=e)
  } else {
    assign("mcodes",  NULL, envir=e)
    assign("mlevels", NULL, envir=e)
  }
  # clipping window
  assign("Wclip",    Wclip,           envir=e)
  # reach of model (could be infinite)
  assign("R",        R,               envir=e)
  # current iteration number 
  assign("irep", 0L,             envir=e)
  # next iteration to be inspected
  assign("inxt",  1L,             envir=e)
  # next transition to be inspected
  assign("tnxt",  1L,             envir=e)
  # proposal type
  assign("proptype", NULL,     envir=e)
  # outcome of proposal
  assign("itype",    NULL,     envir=e)
  # proposal location
  assign("proplocn", NULL,     envir=e)
  # proposal mark 
  assign("propmark", NULL,     envir=e)
  # index of proposal point in existing pattern
  assign("propindx", NULL,     envir=e)
  # Hastings ratio
  assign("numerator",   NULL,  envir=e)
  assign("denominator", NULL,  envir=e)
  # Expression actually evaluated to execute visual debug
  # Expression is evaluated in the environment 'e'
  snoopexpr <-
    expression({
      rslt <- rmhsnoop(Wsim=Wsim, Wclip=Wclip, R=R,
                       xcoords=xcoords,
                       ycoords=ycoords,
                       mlevels=mlevels,
                       mcodes=mcodes,
                       irep=irep,
                       itype=itype,
                       proptype=proptype,
                       proplocn=proplocn,
                       propmark=propmark,
                       propindx=propindx,
                       numerator=numerator,
                       denominator=denominator)
      inxt <- rslt$inxt
      tnxt <- rslt$tnxt
      itype <- if(rslt$accepted) rslt$itype else 0
      storage.mode(tnxt) <-
        storage.mode(inxt) <- storage.mode(itype) <- "integer"
})
  assign("snoopexpr", snoopexpr,  envir=e)
  # callback expression
  assign("callbackexpr", quote(eval(snoopexpr)), envir=e)
  return(e)
}

# visual debug display using base graphics

rmhsnoop <- local({

  rmhsnoop <- function(..., Wsim, Wclip, R,
                       xcoords, ycoords,
                       mlevels, mcodes,
                       irep, itype, 
                       proptype, proplocn, propmark, propindx,
                       numerator, denominator) {
    trap.extra.arguments(..., .Context="In rmhsnoop")
    X <- ppp(xcoords, ycoords, window=Wsim)
    if(!missing(mlevels) && length(mlevels) > 0)
      marks(X) <- factor(mlevels[mcodes+1], levels=mlevels)
    Wclip.orig <- Wclip
    # determine plot arguments
    if(is.mask(Wclip)) {
      parg.Wclip <- list(invert=TRUE, col="grey")
    } else {
      Wclip <- as.psp(Wclip) 
      parg.Wclip <- list(lty=3, lwd=2, col="grey")
    }
    parg.birth <- list(pch=16, cols="green")
    parg.death <- list(pch=4, cols="red", lwd=2)
    parg.birthcircle <- list(col="green", lty=3)
    parg.deathcircle <- list(col="red", lty=3)

    # assemble a layered object representing the state and the proposal
    if(is.null(proptype)) {
      # initial state
      L <- layered(Wsim,
                   Wclip,
                   X)
      layerplotargs(L)$Wclip <- parg.Wclip
      accepted <- TRUE
    } else {
      accepted <- (itype == proptype)
      # add proposal info
      switch(decode.proptype(proptype),
             Reject=
             {
               propname <- "rejected"
               L <- layered(Wsim=Wsim,
                            Wclip=Wclip,
                            X=X)
               layerplotargs(L)$Wclip <- parg.Wclip
             },
             Birth = 
             {
               propname <- "birth proposal"
               U <- ppp(proplocn[1], proplocn[2], window=Wsim)
               D <- if(is.finite(R) && R > 0) {
                 as.psp(disc(R, proplocn))[Wsim]
               } else NULL
               L <- layered(Wsim=Wsim,
                            Wclip=Wclip,
                            PrevState=X,
                            Reach=D,
                            NewPoint=U)
               layerplotargs(L)$Wclip <- parg.Wclip
               layerplotargs(L)$NewPoint <- parg.birth
             },
             Death = 
             {
               propname <- "death proposal"
               # convert from C to R indexing
               propindx <- propindx + 1
               XminI <- X[-propindx]
               XI <- X[propindx]
               D <- if(is.finite(R) && R > 0) {
                 as.psp(disc(R, c(XI$x, XI$y)))[Wsim]
               } else NULL
               L <- layered(Wsim=Wsim,
                            Wclip=Wclip,
                            RetainedPoints=XminI,
                            Reach=D,
                            Deletion=XI)
               layerplotargs(L)$Wclip    <- parg.Wclip
               layerplotargs(L)$Reach    <-  parg.deathcircle
               layerplotargs(L)$Deletion <- parg.death
             },
             Shift = 
             {
               propname <- "shift proposal"
               U <- ppp(proplocn[1], proplocn[2], window=Wsim)
               if(is.finite(R) && R > 0) {
                 DU <- as.psp(disc(R, proplocn))[Wsim]
                 DXI <- as.psp(disc(R, c(XI$x, XI$y)))[Wsim]
               } else { DU <- DXI <- NULL }
               # convert from C to R indexing
               propindx <- propindx + 1
               XminI <- X[-propindx]
               XI <- X[propindx]
               # make layers
               L <- layered(Wsim=Wsim,
                            Wclip=Wclip,
                            OtherPoints=XminI,
                            ReachAfter=DU,
                            AfterShift=U,
                            ReachBefore=DXI,
                            BeforeShift=XI)
               layerplotargs(L)$Wclip       <- parg.Wclip
               layerplotargs(L)$ReachAfter  <- parg.birthcircle
               layerplotargs(L)$AfterShift  <- parg.birth
               layerplotargs(L)$ReachBefore <- parg.deathcircle
               layerplotargs(L)$BeforeShift <- parg.death
             },
             stop("Unrecognised proposal type")
             )
    }
    header <- c(paste("Iteration", irep),
                propname,
                paste("Hastings ratio =",
                      signif(numerator, 4), "/", signif(denominator, 4)))
    info <- list(irep=irep,
                 Wsim=Wsim,
                 Wclip=Wclip.orig,
                 X=X,
                 proptype=proptype,
                 proplocn=proplocn,
                 propindx=propindx,
                 propmark=propmark,
                 accepted=accepted,
                 numerator=numerator,
                 denominator=denominator)
    inspectProposal(L, info, title=header)
  }

  decode.proptype <- function(n) {
    if(n < 0 || n > 3) stop(paste("Unrecognised proposal type:", n))
    switch(n+1, "Reject", "Birth", "Death", "Shift")
  }
  encode.proptype <- function(s) {
    switch(s, Reject=0, Birth=1, Death=2, Shift=3)
  }
  
  inspectProposal <- function(X, info, ..., title) {
    if(missing(title)) title <- short.deparse(substitute(X))
    if(!inherits(X, "layered"))
      X <- layered(X)
    lnames <- names(X)
    if(sum(nzchar(lnames)) != length(X))
      lnames <- paste("Layer", seq_len(length(X)))
    # Find window and bounding box (validates X)
    W <- as.owin(X)
    BX <- as.rectangle(W)
    # Initialise environment for state variables etc
    # This environment is accessible to the panel button functions
    en <- new.env()
    assign("X", X, envir=en)
    assign("W", W, envir=en)
    assign("BX", BX, envir=en)
    assign("zoomfactor", 1L, envir=en)
    midX <- unlist(centroid.owin(BX))
    assign("midX", midX, envir=en)
    assign("zoomcentre", midX, envir=en)
    assign("irep", info$irep, envir=en)
    assign("inxt", info$irep+1, envir=en) 
    assign("tnxt", -1, envir=en)
    assign("accepted", info$accepted, envir=en)
    assign("proplocn", info$proplocn, envir=en)
    assign("info", info, envir=en)
    # Build interactive panel
    # Start with data panel
    P <- simplepanel(title,
                     BX,
                     list(Data=BX),
                     list(Data=dataclickfun),
                     list(Data=dataredrawfun),
                     snoopexit,
                     en)
    # Add pan buttons
    margin <- max(sidelengths(BX))/4
    panelwidth <- sidelengths(BX)[1]/2
    P <- grow.simplepanel(P, "top", margin, navfuns["Up"], aspect=1)
    P <- grow.simplepanel(P, "bottom", margin, navfuns["Down"], aspect=1)
    P <- grow.simplepanel(P, "left", margin, navfuns["Left"], aspect=1)
    P <- grow.simplepanel(P, "right", margin, navfuns["Right"], aspect=1)
    # Zoom/Pan buttons at right
    P <- grow.simplepanel(P, "right", panelwidth, zoomfuns)
    # Accept/reject buttons at top
    P <- grow.simplepanel(P, "top", margin, accept.clicks, accept.redraws)
    # Dump/print buttons at bottom 
    P <- grow.simplepanel(P, "bottom", margin, dumpfuns)
    # Jump controls at left
    maxchars <- max(4, nchar(names(jump.clicks)))
    P <- grow.simplepanel(P, "left", panelwidth * maxchars/6, jump.clicks)
    # go
    rslt <- run.simplepanel(P)
    clear.simplepanel(P)
    rm(en)
    return(rslt)
  }


# button control functions
  zoomfuns <- 
    rev(list(
             "Zoom In"=function(env, xy) {
               z <- get("zoomfactor", envir=env)
               assign("zoomfactor", z * 2, envir=env)
               return(TRUE)
             },
             "Zoom Out"=function(env, xy) {
               z <- get("zoomfactor", envir=env)
               assign("zoomfactor", z / 2, envir=env)
               return(TRUE)
             },
             "At Proposal"=function(env, xy) {
               proplocn <- get("proplocn", envir=env)
               assign("zoomcentre", proplocn, envir=env)
               return(TRUE)
             },
             Reset=function(env, xy) {
               assign("zoomfactor", 1L, envir=env)
               midX <- get("midX", envir=env)
               assign("zoomcentre", midX, envir=env)
               return(TRUE)
             }))
                           
  navfuns <-
    list(
         Left = function(env, xy) {
           zoom <- get("zoomfactor", envir=env)
           ce <- get("zoomcentre", envir=env)
           BX <- get("BX", envir=env)
           width <- sidelengths(BX)[1]
           stepsize <- (width/4)/zoom
           ce <- ce - c(stepsize, 0)
           assign("zoomcentre", ce, envir=env)
           return(TRUE)
         },
         Right = function(env, xy) {
           zoom <- get("zoomfactor", envir=env)
           ce <- get("zoomcentre", envir=env)
           BX <- get("BX", envir=env)
           width <- sidelengths(BX)[1]
           stepsize <- (width/4)/zoom
           ce <- ce + c(stepsize, 0)
           assign("zoomcentre", ce, envir=env)
           return(TRUE)
         },
         Up = function(env, xy) {
           zoom <- get("zoomfactor", envir=env)
           ce <- get("zoomcentre", envir=env)
           BX <- get("BX", envir=env)
           height <- sidelengths(BX)[2]
           stepsize <- (height/4)/zoom
           ce <- ce + c(0, stepsize)
           assign("zoomcentre", ce, envir=env)
           return(TRUE)
         },
         Down = function(env, xy) {
           zoom <- get("zoomfactor", envir=env)
           ce <- get("zoomcentre", envir=env)
           BX <- get("BX", envir=env)
           height <- sidelengths(BX)[2]
           stepsize <- (height/4)/zoom
           ce <- ce - c(0, stepsize)
           assign("zoomcentre", ce, envir=env)
           return(TRUE)
         })

  accept.clicks <-
    rev(list(
             Accept=function(env, xy) {
               assign("accepted", TRUE, envir=env)
               return(TRUE)
             },
             Reject=function(env, xy) {
               assign("accepted", FALSE, envir=env)
               return(TRUE)
             }))

  accept.redraws <-
    rev(list(
             Accept=function(button, name, env) {
               accepted <- get("accepted", envir=env)
               if(accepted) {
                 plot(button, add=TRUE, col="green")
               } else {
                 plot(button, add=TRUE)
               }
               text(centroid.owin(button), labels=name)
             },
             Reject=function(button, name, env) {
               accepted <- get("accepted", envir=env)
               if(accepted) {
                 plot(button, add=TRUE)
               } else {
                 plot(button, add=TRUE, col="pink")
               }
               text(centroid.owin(button), labels=name)
             }))
             
  jump.clicks <-
    rev(list(
             "Next Iteration"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+1, envir=env)
               return(FALSE)
             },
             "Skip 10"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+10, envir=env)
               return(FALSE)
             },
             "Skip 100"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+100, envir=env)
               return(FALSE)
             },
             "Skip 1000"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+1000, envir=env)
               return(FALSE)
             },
             "Skip 10,000"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+10000, envir=env)
               return(FALSE)
             },
             "Skip 100,000"=function(env, xy) {
               irep <- get("irep", envir=env)
               assign("inxt", irep+100000, envir=env)
               return(FALSE)
             },
             "Next Birth"=function(env, xy) {
               assign("inxt", -1, envir=env)
               assign("tnxt", encode.proptype("Birth"), envir=env)
               return(FALSE)
             },
             "Next Death"=function(env, xy) {
               assign("inxt", -1, envir=env)
               assign("tnxt", encode.proptype("Death"), envir=env)
               return(FALSE)
             },
             "Next Shift"=function(env, xy) {
               assign("inxt", -1, envir=env)
               assign("tnxt", encode.proptype("Shift"), envir=env)
               return(FALSE)
             },
             "Exit Debugger"=function(env, xy) {
               assign("inxt", -1L, envir=env)
               return(FALSE)
             }))

  dataclickfun <- function(env, xy) {
    # function for handling clicks in the data window
    z <- get("zoomfactor", envir=env)
    ce <- get("zoomcentre", envir=env)
    midX <- get("midX", envir=env)
    ce <- ce + (unlist(xy) - midX)/z
    assign("zoomcentre", ce, envir=env)
    return(TRUE)
  }

  dataredrawfun <- function(button, name, env) {                             
    # redraw data window
    X <- get("X", envir=env)
    BX <- get("BX", envir=env)
    W <- get("W", envir=env)
    midX <- get("midX", envir=env)
    z <- get("zoomfactor", envir=env)
    ce <- get("zoomcentre", envir=env)
    scaleX <- shift(affine(shift(X, -ce), diag(c(z,z))), unlist(midX))
    scaleW <- shift(affine(shift(W, -ce), diag(c(z,z))), unlist(midX))
    scaleX <- scaleX[, BX]
    scaleW <- intersect.owin(scaleW, BX, fatal=FALSE)
    # redraw data in 'BX' 
    if(!is.null(scaleW)) {
      if(z == 1 && is.rectangle(scaleW)) {
        plot(scaleW, add=TRUE, lwd=2)
      } else {
        plot(BX, add=TRUE, lty=3, border="red")
        if(!identical(BX, scaleW))
          plot(scaleW, add=TRUE, invert=TRUE)
      }
    }
    if(!is.null(scaleX))
      plot(scaleX, add=TRUE)
    invisible(NULL)
  }

# functions to dump the current state, etc
  dumpfuns <- list(
                   "Dump to file"=function(env, xy) {
                     irep <- get("irep", envir=env)
                     X <- get("X", envir=env)
                     xname <- paste("dump", irep, sep="")
                     assign(xname, X)
                     fname <- paste(xname, ".rda", sep="")
                     eval(substitute(save(x, file=y, compress=TRUE), 
                                     list(x=xname, y=fname)))
                     cat(paste("Saved to", sQuote(fname), "\n"))
                     return(TRUE)
                   },
                   "Print Info"=function(env, xy) {
                     info <- get("info", envir=env)
                     will.accept <- get("accepted", envir=env)
                     with(info, {
                       cat(paste("Iteration", irep, "\n"))
                       cat("Simulation window:\n")
                       print(Wsim)
                       cat("Clipping window:\n")
                       print(Wclip)
                       cat("Current state:\n")
                       print(X)
                       propname <- decode.proptype(proptype)
                       cat(paste("Proposal type:", propname, "\n"))
                       prxy <- function(z) paren(paste(z, collapse=", "))
                       switch(propname,
                              Reject = { },
                              Birth = {
                                cat(paste("Birth of new point at location",
                                          prxy(proplocn), "\n"))
                              },
                              Death = {
                                Xi <- X[propindx]
                                cat(paste("Death of data point", propindx,
                                          "located at",  
                                          prxy(as.numeric(coords(Xi))),
                                          "\n"))
                              },
                              Shift = {
                                Xi <- X[propindx]
                                cat(paste("Shift data point",
                                          propindx,
                                          "from current location",
                                          prxy(as.numeric(coords(Xi))),
                                          "to new location",
                                          prxy(proplocn),
                                          "\n"))
                              })
                       cat(paste("Hastings ratio = ",
                                 numerator, "/", denominator,
                                 "=", numerator/denominator, "\n"))
                       cat(paste("Fate of proposal:",
                                 if(will.accept) "Accepted" else "Rejected",
                                 "\n"))
                       return(TRUE)
                     })
                   })
  
# function to determine return value
                             
  snoopexit <- function(env) {
    ans <- eval(quote(list(inxt=inxt,
                           tnxt=tnxt,
                           accepted=accepted)),
                envir=env)
    return(ans)
  }
                             
  testit <- function() {
    rmhsnoop(Wsim=owin(), Wclip=square(0.7), R=0.1,
             xcoords=runif(40),
             ycoords=runif(40),
             mlevels=NULL, mcodes=NULL,
             irep=3, itype=1,
             proptype=1, proplocn=c(0.5, 0.5), propmark=0, propindx=0,
             numerator=42, denominator=24)
  }
                             
  rmhsnoop
})


