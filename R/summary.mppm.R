#
# summary.mppm.R
#
# $Revision: 1.11 $  $Date: 2015/02/06 10:26:46 $
#


summary.mppm <- function(object, ..., brief=FALSE) {
  # y will be the summary 
  y <- object[c("Call", "Info", "Inter", "trend", "iformula",
                "npat", "maxlogpl")]
  y$brief <- brief

  Info  <- object$Info
  Inter <- object$Inter
  FIT   <- object$Fit$FIT
  moadf <- object$Fit$moadf

  y$Fit <- object$Fit[c("fitter", "use.gam", "fmla", "Vnamelist")]
  y$Fit$FIT <- summary(FIT)
  y$Fit$moadf <- list(nrow=nrow(moadf), colnames=colnames(moadf))
  
  ninteract    <- Inter$ninteract
  interaction  <- Inter$interaction
  iused        <- Inter$iused
  itags        <- Inter$itags
  processnames <- Inter$processes
  constant     <- Inter$constant
  trivial      <- Inter$trivial

  npat      <- y$npat
  iformula  <- y$iformula
  Vnamelist <- y$Fit$Vnamelist
  allVnames <- unlist(Vnamelist)
  poistags  <- itags[trivial]

#  rownames  <- y$Info$rownames
  
  switch(y$Fit$fitter,
         gam=,
         glm={
           y$coef <- co <- coef(FIT)
           systematic <- !(names(co) %in% c(allVnames, poistags))
           y$coef.syst <- co[systematic]
         })

  # model depends on covariates
  y$depends.covar <- Info$has.covar && (length(Info$used.cov.names) > 0)


  ### Interactions 
  # model is Poisson 
  y$poisson <- all(trivial[iused])
  # Determine how complicated the interactions are:
  # (1) is the interaction formula of the form ~ tag + tag + ... + tag
  isimple  <- identical(sort(variablesinformula(iformula)),
                        sort(termsinformula(iformula)))
  # (2) is it of the form ~tag 
  trivialformula <- (isimple && ninteract == 1)
  # (3) is it of the form ~tag where the interaction is the same in each row
  fixedinteraction <- trivialformula && constant
  
  ### Determine printing of interactions, accordingly ###
  iprint <- list()
  if(fixedinteraction) {    
    # exactly the same interaction for all patterns
    interaction <- interaction[1,1,drop=TRUE]
    fi.all <- fii(interaction, co, Vnamelist[[1]]) 
    iprint <- list("Interaction for all patterns"=fi.all)
    printeachrow <- FALSE
    toohard      <- FALSE
  } else if(trivialformula) {
    # same type of process for all patterns
    pname <-  unlist(processnames)[iused]
    iprint <- list("Interaction for each pattern" = pname)
    printeachrow <- TRUE
    toohard      <- FALSE
  } else if(isimple && all(constant)) {
    # several interactions involved, each of which is the same for all patterns
    iprint <- list("Interaction formula"=iformula,
                   "Interactions defined for each pattern"=NULL)
    for(j in (1:ninteract)[iused]) {
      name.j <- paste("Interaction", sQuote(itags[j]))
      int.j <- Inter$interaction[1,j,drop=TRUE]
      Vnames.j <- Vnamelist[[j]]
      fii.j <- fii(int.j, co, Vnames.j)
      extra.j <- list(fii.j)
      names(extra.j) <- name.j
      iprint <- append(iprint, extra.j)
    }
    printeachrow <- FALSE
    toohard      <- FALSE
  } else {
    # general case
    # determine which interaction(s) are active on each row
    active <- active.interactions(object)
    if(ninteract > 1 || !all(active)) 
      iprint <- list("Active interactions"=active)
    printeachrow <- TRUE
    toohard <- any(rowSums(active) > 1)
  }

  y$ikind <- list(
                  isimple           =isimple,
                  trivialformula    =trivialformula,
                  fixedinteraction  =fixedinteraction,
                  toohard           =toohard,
                  printeachrow      =printeachrow)

  if(toohard)
    iprint <- append(iprint,
                     list("(Sorry, cannot interpret fitted interactions)"))
  else if(printeachrow) {
    subs <- subfits(object, what="interactions")
    names(subs) <- paste("Interaction", 1:npat)
    iprint <- append(iprint, subs)
  }

  y$iprint <- iprint

  class(y) <- c("summary.mppm", class(list))
  return(y)
}


print.summary.mppm <- function(x, ..., brief=x$brief) {
  # NB: x is an object of class "summary.mppm"
  npat <- x$npat
  Inter <- x$Inter
#  ninteract   <- Inter$ninteract
#  interaction   <- Inter$interaction
#  iused     <- Inter$iused
#  constant <- Inter$constant
#  iformula <- x$iformula
#  processnames   <- Inter$processes
  itags   <- Inter$itags
  trivial  <- Inter$trivial

  FIT <- x$Fit$FIT
  Vnamelist <- x$Fit$Vnamelist
  
#  allVnames <- unlist(Vnamelist)
#  poistags <- itags[trivial]

  terselevel <- spatstat.options("terse")
#  rownames <- x$Info$rownames

  splat("Point process model fitted to", npat, "point patterns")
  if(waxlyrical('gory', terselevel))
    splat("Call:", x$Call$callstring)
  splat("Log trend formula:", pasteFormula(x$trend))
  switch(x$Fit$fitter,
         gam=,
         glm={
           cat("Fitted trend coefficients:\n")
           print(x$coef.syst)
           co <- coef(FIT)
         })

  if(!brief && waxlyrical('extras', terselevel)) {
    cat("All fitted coefficients:\n")
    print(co)
  }
    
  parbreak(terselevel)


  ### Print interaction information ###
  if(waxlyrical('extras', terselevel)) {
    iprint <- x$iprint 
    nama <- names(iprint) %orifnull% rep("", length(iprint))
    for(i in seq_along(iprint)) {
      nami <- nama[i]
      vali <- iprint[[i]]
      if(brief && is.matrix(vali))
        vali <- paren(paste(nrow(vali), "x", ncol(vali), "matrix"))
      if(nami != "") {
        inline <- inherits(vali, "formula") ||
                  is.character(vali) ||
                  (brief && inherits(vali, "fii"))
        if(inline) cat(paste0(nami, ":\t")) else splat(paste0(nami, ":"))
      }
      if(!is.null(vali)) {
        if(inherits(vali, "fii")) {
          print(vali, tiny=brief)
        } else if(is.character(vali)) {
          splat(vali)
        } else {
          print(vali)
        } 
      }
      parbreak(terselevel)
    }
  }

  if(!brief && waxlyrical('gory', terselevel)) {
    splat("--- Gory details: ---")
    splat("Combined data frame has", x$Fit$moadf$nrow, "rows")
    print(FIT)
  }
  invisible(NULL)
}

