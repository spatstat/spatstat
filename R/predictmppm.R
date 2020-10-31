#
#    predictmppm.R
#
#	$Revision: 1.16 $	$Date: 2020/10/31 13:50:03 $
#
#
# -------------------------------------------------------------------

predict.mppm <- local({

  predict.mppm <- function(object, ..., newdata=NULL, type=c("trend", "cif"),
                           ngrid=40, locations=NULL, verbose=FALSE) {
    ##
    ##	'object' is the output of mppm()
    ##
    model <- object
    verifyclass(model, "mppm")
    ## 
    isMulti <- is.multitype(model)
    depends.on.row <- summary(model)$depends.on.row
    npat.old <- model$npat
    ##
    ## ......................................................................
    if(verbose)
      cat("Inspecting arguments...")
    ## ......................................................................
    ##
    ##   hidden arguments
    selfcheck <- resolve.defaults(list(...), list(selfcheck=FALSE))$selfcheck
    ##
    ##  Argument 'type'
    ##             
    type <- pickoption("type", type, c(trend="trend",
                                       lambda="cif",
                                       cif="cif"), multi=TRUE)
    want.trend <- "trend" %in% type
    want.cif   <- "cif"   %in% type
    ##
    ##  Argument 'newdata'
    ##
    use.olddata <- is.null(newdata)
    if(use.olddata) {
      newdata <- model$data
      newdataname <- "Original data"
    } else {
      stopifnot(is.data.frame(newdata) || is.hyperframe(newdata))
      newdataname <- sQuote("newdata")
      if(depends.on.row && nrow(newdata) != npat.old)
        stop(paste("'newdata' must have the same number of rows",
                   "as the original 'data' argument",
                   paren(paste("namely", npat.old)),
                   "because the model depends on the row index"),
             call.=FALSE)
    }
    ##
    ##   Argument 'locations'
    ##
    if(is.hyperframe(locations)) 
      locations <- locations[,1,drop=TRUE]
    if(is.list(locations))
      cls <- unique(sapply(locations, class))

    loctype <-
      if(is.null(locations)) "null" else
      if(is.data.frame(locations))  "data.frame" else
      if(is.list(locations)) {
        if(any(c("ppp", "quad") %in% cls)) "points"
        else if("owin" %in% cls) {
          if(all(sapply(locations, is.mask)))
            "mask"
          else
            "window"
        } else "unknown"
      } else "unknown"

    ## ......................................................................
    if(verbose)
      cat("done.\nDeciding type of locations for prediction...")
    ## ......................................................................
    
    need.grid <- switch(loctype,
                        null      =TRUE,
                        data.frame=FALSE,
                        points    =FALSE,
                        mask      =FALSE,
                        window    =TRUE,
                        unknown   =stop("Unrecognised format for locations"))
    
    make.image <- need.grid || (loctype == "mask")
    ##  
    locationvars <- c("x", "y", "id", if(isMulti) "marks" else NULL)
    ##  
    ##
    if(need.grid) {
      ## prediction on a grid is required
      if(is.data.frame(newdata))
        stop(paste("Cannot predict model on a grid;", newdataname,
                   "are a data frame"))
    } else {
      ## prediction at  `locations' is required
      if(is.hyperframe(newdata)) {
        ## check consistency between locations and newdata
        nloc <- length(locations)
        nnew <- summary(newdata)$ncases
        if(nloc != nnew)
          stop(paste("Length of argument", sQuote("locations"), paren(nloc),
                     "does not match number of rows in",
                     newdataname, paren(nnew)))
      } else {
        ## newdata is a data frame
        if(!is.data.frame(locations)) 
          stop(paste(newdataname,
                     "is a data frame; locations must be a data frame"))
        else {
          stopifnot(nrow(locations) == nrow(newdata))
          dup <- names(newdata) %in% names(locations)
          if(any(dup))
            for(nam in names(newdata)[dup])
              if(!isTRUE(all.equal(newdata[,nam], locations[,nam])))
                stop(paste("The data frames newdata and locations",
                           "both have a column called", sQuote(nam),
                           "but the entries differ"))
          nbg <- !(locationvars %in% c(names(newdata),names(locations)))
          if(any(nbg))
            stop(paste(ngettext(sum(nbg), "Variable", "Variables"),
                       commasep(locationvars[nbg]),
                       "not provided"))
          ## merge the two data frames
          newdata <- cbind(newdata[,!dup], locations)
          locations <- NULL
        }
      }
    }

    ## ......................................................................
    if(verbose)
      cat("done.\nExtracting details of point process model...")
    ## ......................................................................

    ## extract fitted glm/gam/glmm object
    FIT <- model$Fit$FIT
    MOADF <- model$Fit$moadf
    
    ## extract names of interaction variables
    Vnamelist <- model$Fit$Vnamelist
    vnames <- unlist(Vnamelist)
    
    ## determine which interaction is applicable on each row
    interactions <- model$Inter$interaction
    ninter <- if(is.hyperframe(interactions)) nrow(interactions) else 1
    nnew <- nrow(newdata)
    if(ninter != nnew && ninter != 1) {
      if(!all(model$Inter$constant))
        stop(paste("Number of rows of newdata", paren(nnew),
                   "does not match number of interactions in model",
                   paren(ninter)))
      interactions <- interactions[1, ]
    }

    ## extract possible types, if model is multitype
    if(isMulti) {
      levlist <- unique(lapply(data.mppm(model), levelsofmarks))
      if(length(levlist) > 1)
        stop("Internal error: the different point patterns have inconsistent marks",
             call.=FALSE)
      marklevels <- levlist[[1L]]
     } else marklevels <- list(NULL) # sic
    
    ## ......................................................................
    if(verbose) {
      cat("done.\n")
      if(use.olddata) splat("Using original hyperframe of data") else 
      splat("newdata is a", if(is.data.frame(newdata)) "data frame" else "hyperframe")
    }
    ## ......................................................................
    ##
    if(is.data.frame(newdata)) {
      ##  
      ##              newdata is a DATA FRAME
      ##
      if(need.grid)
        stop("Cannot predict model on a grid; newdata is a data frame")
      if(verbose)
        cat("Computing prediction..")
      ## use newdata as covariates
      nbg <- !(locationvars %in% names(newdata))
      if(any(nbg))
        stop(paste(ngettext(sum(nbg), "variable", "variables"),
                   commasep(locationvars[nbg]),
                   "not provided"))
      ## create output data frame
      answer <- as.data.frame(matrix(, nrow=nrow(newdata), ncol=0),
                              row.names=row.names(newdata))
      if(want.trend) {
        ## add interaction components, set to zero (if any)
        if(length(vnames) > 0)
          newdata[, vnames] <- 0
        ## compute fitted values
        answer$trend <- Predict(FIT, newdata=newdata, type="response")
      }
      if(want.cif) {
        if(is.poisson(object)) {
          ## cif = trend
          answer$cif <- if(want.trend) answer$trend else
                        Predict(FIT, newdata=newdata, type="response")
        } else {
          warning("Computation of the cif is not yet implemented when newdata is a data frame")
          ## split data frame by 'id'
          ## compute interaction components using existing point patterns
          ## compute fitted values
        }
      }
      if(verbose) cat("done.\n")
      return(answer)
    }
  
    ## ......................................................................
    ##        newdata is a HYPERFRAME
    ##
    if(verbose)
      cat("Building data for prediction...")
    sumry.new <- summary(newdata)
    npat.new <- sumry.new$ncases
    ## name of response point pattern in model
    Yname <- model$Info$Yname
    ##
    ## Determine response point patterns if known.
    ## Extract from newdata if available
    ## Otherwise from the original data if appropriate
    if(verbose)
      cat("(responses)...")
    Y <- if(Yname %in% sumry.new$col.names) 
      newdata[, Yname, drop=TRUE, strip=FALSE]
    else if(npat.new == npat.old)
      data[, Yname, drop=TRUE, strip=FALSE]
    else NULL
    ##
    if(want.cif && is.null(Y))
      stop(paste("Cannot compute cif:",
                 "newdata does not contain column", dQuote(Yname),
                 "of response point patterns"))
    ##
    ## Determine windows for prediction 
    if(verbose)
      cat("(windows)...")
    Wins <- if(!need.grid)
      lapply(locations, as.owin, fatal=FALSE)
    else if(!is.null(Y))
      lapply(Y, as.owin, fatal=FALSE)
    else NULL
    if(is.null(Wins) || any(sapply(Wins, is.null)))
      stop("Cannot determine windows where predictions should be made")
    ##
    ##
    if(is.null(Y)) {
      ## only want trend; empty patterns will do
      Y <- lapply(Wins, emptypattern)
    }
    
    ## ensure Y contains data points only 
    if(is.quad(Y[[1]]))
      Y <- lapply(Y, getElement, name="data")

    ## Determine locations for prediction
    if(need.grid) {
      ## Generate grids of dummy locations 
      if(verbose)
        cat("(grids)...")
      Gridded <- lapply(Wins, gridsample, ngrid=ngrid)
      Dummies   <- lapply(Gridded, getElement, name="D")
      Templates <- lapply(Gridded, getElement, name="I")
    } else {
      ## locations are given somehow
      if(verbose)
        cat("(locations)...")
      switch(loctype,
             points = {
               Dummies <- locations
             },
             mask = {
               Dummies <- lapply(locations, punctify)
               Templates <- lapply(locations, as.im)
             },
             stop("Internal error: illegal loctype"))
    }

    ## ..........................................
    ## ............... PREDICTION ...............
    ## ..........................................

    ## initialise hyperframe of predicted values
    Answer <- newdata[,integer(0),drop=FALSE]
    if(depends.on.row) Answer$id <- factor(levels(MOADF$id))

    ## Loop over possible types, or execute once:
    ## ///////////////////////////////////////////
  
    for(lev in marklevels) {

      ## Pack prediction locations into quadschemes
      if(verbose) {
        cat("Building quadschemes")
        if(isMulti) cat(paste("with mark", lev))
        cat("...")
      }
      
      if(isMulti) {
        ## assign current mark level to all dummy points
        flev <- factor(lev, levels=marklevels)
        Dummies <- lapply(Dummies, "marks<-", value=flev)
      }
      Quads <- mapply(quad, data=Y, dummy=Dummies,
                      SIMPLIFY=FALSE, USE.NAMES=FALSE)
      
      ## Insert quadschemes into newdata
      newdata[, Yname] <- Quads
    
      ## compute the Berman-Turner frame
      if(verbose)
        cat("done.\nStarting prediction...(Berman-Turner frame)...")
      moadf <- mppm(formula     = model$formula,
                    data        = newdata,
                    interaction = interactions,
                    iformula    = model$iformula,
                    random      = model$random,
                    use.gam     = model$Fit$use.gam,
                    correction  = model$Info$correction,
                    rbord       = model$Info$rbord,
                    backdoor    = TRUE)
      ## compute fitted values
      if(verbose)
        cat("(glm prediction)...")
      values <- moadf[, locationvars]
      if(want.cif)
        values$cif <- Predict(FIT, newdata=moadf, type="response")
      if(want.trend) {
        if(length(vnames) == 0) {
          ## Poisson model: trend = cif 
          values$trend <-
            if(want.cif) values$cif else
            Predict(FIT, newdata=moadf, type="response")
        } else {
          ## zero the interaction components
          moadf[, vnames] <- 0
          ## compute fitted values
          values$trend <- Predict(FIT, newdata=moadf, type="response")
        }
      }
      if(verbose)
        cat("done.\nReshaping results...")
      ##
      ## Reshape results
      ## separate answers for each image
      values <- split(values, values$id)
      ## 
      Trends <- list()
      Lambdas <- list()
      if(!make.image) {
        if(verbose)
          cat("(marked point patterns)...")
        ## values become marks attached to locations
        for(i in seq(npat.new)) {
          Val <- values[[i]]
          Loc <- Dummies[[i]]
          isdum <- !is.data(Quads[[i]])
          if(selfcheck)
            if(length(isdum) != length(Val$trend))
              stop("Internal error: mismatch between data frame and locations")
          if(want.trend)
            Trends[[i]] <- Loc %mark% (Val$trend[isdum])
          if(want.cif)
            Lambdas[[i]] <- Loc %mark% (Val$cif[isdum])
        }
      } else {
        if(verbose)
          cat("(pixel images)...")
        ## assign values to pixel images
        for(i in seq(npat.new)) {
          values.i <- values[[i]]
          Q.i <- Quads[[i]]
          values.i <- values.i[!is.data(Q.i), ]
          Template.i <- Templates[[i]]
          ok.i <- !is.na(Template.i$v)
          if(sum(ok.i) != nrow(values.i))
            stop("Internal error: mismatch between data frame and image")
          if(selfcheck) {
            dx <- rasterx.im(Template.i)[ok.i] - values.i$x
            dy <- rastery.im(Template.i)[ok.i] - values.i$y
            cat(paste("i=", i, "range(dx) =", paste(range(dx), collapse=", "),
                      "range(dy) =", paste(range(dy), collapse=", "), "\n"))
          }
          if(want.trend) {
            Trend.i <- Template.i
            Trend.i$v[ok.i] <- values.i$trend
            Trends[[i]] <- Trend.i
          }
          if(want.cif) {
            Lambda.i <- Template.i
            Lambda.i$v[ok.i] <- values.i$cif
            Lambdas[[i]] <- Lambda.i
          }
        }
      }
      if(verbose)
        cat("done reshaping.\n")
      
      if(want.trend) {
        trendname <- paste0("trend", lev)
        Answer[,trendname] <- Trends
      }
      if(want.cif) {
        cifname <- paste0("cif", lev)
        Answer[,cifname] <- Lambdas
      }
    }   ## ///////////  end loop over possible types //////////////////
    
    return(Answer)
  }

  ## helper functions
  emptypattern <- function(w) { ppp(numeric(0), numeric(0), window=w) }

  levelsofmarks <- function(X) { levels(marks(X)) }
      
  gridsample <- function(W, ngrid) {
    masque <- as.mask(W, dimyx=ngrid)
    xx <- raster.x(masque)
    yy <- raster.y(masque)
    xpredict <- xx[masque$m]
    ypredict <- yy[masque$m]
    Dummy <- ppp(xpredict, ypredict, window=W)
    Image <- as.im(masque)
    return(list(D=Dummy, I=Image))
  }

  punctify <- function(M) { 
    xx <- raster.x(M)
    yy <- raster.y(M)
    xpredict <- xx[M$m]
    ypredict <- yy[M$m]
    return(ppp(xpredict, ypredict, window=M))
  }

  Predict <- function(object, newdata, type=c("link", "response")) {
    type <- match.arg(type)
    if(inherits(object, "glmmPQL")) {
      class(object) <- class(object)[-1L]
      pred <- predict(object, newdata=newdata)
      if(type == "response") pred <- object$family$linkinv(pred)
    } else {
      pred <- predict(object, newdata=newdata, type=type)
    }
    return(as.numeric(pred))
  }
  
  predict.mppm
})
