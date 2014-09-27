#
#	tstat.R		Estimation of T function
#
#	$Revision: 1.7 $	$Date: 2013/04/25 06:37:43 $
#

Tstat <- local({
  
  # helper functions
  edgetri.Trans <- function(X, triid, trim=spatstat.options("maxedgewt")) {
    triid <- as.matrix(triid)
    ntri <- nrow(triid)
    if(ntri == 0) return(numeric(0))
    W <- rescue.rectangle(as.owin(X))
    if(W$type != "rectangle")
      stop("Translation correction is only implemented for rectangular windows")
    x <- matrix(X$x[triid], nrow=ntri)
    y <- matrix(X$y[triid], nrow=ntri)
    dx <- apply(x, 1, function(z) diff(range(z)))
    dy <- apply(y, 1, function(z) diff(range(z)))
    wide <- diff(W$xrange)
    high <- diff(W$yrange)
    weight <- wide * high/((wide - dx) * (high - dy))
    weight <- pmin.int(trim, weight)
    return(weight)
  }
  # helper function
  implemented.for.T <- function(correction, windowtype, explicit) {
    rect <- (windowtype == "rectangle")
    if(any(correction == "best")) {
      # select best available correction
      correction <- if(rect) "translate" else "border"
    } else {
      # available selection of edge corrections depends on window
      if(!rect) {
        tra <- (correction == "translate") 
        if(any(tra)) {
          whinge <- "Translation correction is only implemented for rectangular windows"
          if(explicit) {
            if(all(tra)) stop(whinge) else warning(whinge)
          }
          correction <- correction[!tra]
        }
      }
    }
    return(correction)
  }
  # .......... main function ....................
  Tstat <- function(X, ..., r=NULL, rmax=NULL,
                    correction=c("border", "translate"),
                    ratio=FALSE,
                    verbose=TRUE)
    {
      verifyclass(X, "ppp")
      rfixed <- !is.null(r) 
      npts <- npoints(X)
      W <- Window(X)
      areaW <- area(W)
      lambda <- npts/areaW
      lambda2 <- (npts * (npts - 1))/(areaW^2)
      lambda3 <- (npts * (npts - 1) * (npts - 2))/(areaW^3)

      rmaxdefault <- if(!is.null(rmax)) rmax else rmax.rule("K", W, lambda)
      breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
      r <- breaks$r
      rmax <- breaks$max

      # choose correction(s)
      correction.given <- !missing(correction) && !is.null(correction)
      if(!correction.given)
        correction <- c("border", "bord.modif", "translate")
      correction <- pickoption("correction", correction,
                               c(none="none",
                                 border="border",
                                 "bord.modif"="bord.modif",
                                 trans="translate",
                                 translate="translate",
                                 translation="translate",
                                 best="best"),
                               multi=TRUE)
      correction <- implemented.for.T(correction, W$type, correction.given)
  
      # recommended range of r values
      alim <- c(0, min(rmax, rmaxdefault))

      # this will be the output data frame
      TT <- data.frame(r=r, theo= (pi/2) * (pi - 3 * sqrt(3)/4) * r^4)
      desc <- c("distance argument r", "theoretical Poisson %s")
      TT <- fv(TT, "r", quote(T(r)),
               "theo", , alim, c("r","%s[pois](r)"), desc, fname="T")

      # save numerator and denominator?
      if(ratio) {
        denom <- lambda2 * areaW
        numT <- eval.fv(denom * TT)
        denT <- eval.fv(denom + TT * 0)
        attributes(numT) <- attributes(denT) <- attributes(T)
        attr(numT, "desc")[2] <- "numerator for theoretical Poisson %s"
        attr(denT, "desc")[2] <- "denominator for theoretical Poisson %s"
      }
  
      # identify all close pairs
      rmax <- max(r)
      close <- closepairs(X, rmax, ordered=FALSE)
      I <- close$i
      J <- close$j
      DIJ <- close$d

      nI <- length(I)
  
      # estimate computation time
      if(verbose) {
        nTmax <- nI * (nI-1) /2
        esttime <- exp(1.25 * log(nTmax) - 21.5)
        message(paste("Searching", nTmax, "potential triangles;",
                      "estimated time", codetime(esttime)))
      }

      # find triangles with their diameters
      tri <- trianglediameters(I, J, DIJ, nvert=npts)
      stopifnot(identical(colnames(tri), c("i", "j", "k", "diam")))
      # reassemble so each triangle appears 3 times, once for each vertex
      II <- with(tri, c(i, j, k))
      DD <- with(tri, rep.int(diam, 3))
  
      if(any(correction == "none")) {
        # uncorrected! For demonstration purposes only!
        wh <- whist(DD, breaks$val)  # no weights
        numTun <- cumsum(wh)
        denTun <- lambda3 * areaW
        # uncorrected estimate of T
        Tun <- numTun/denTun
        TT <- bind.fv(TT, data.frame(un=Tun), "hat(%s)[un](r)",
                      "uncorrected estimate of %s",
                      "un")
        if(ratio) {
          # save numerator and denominator
          numT <- bind.fv(numT, data.frame(un=numTun), "hat(%s)[un](r)",
                          "numerator of uncorrected estimate of %s",
                          "un")
          denT <- bind.fv(denT, data.frame(un=denTun), "hat(%s)[un](r)",
                          "denominator of uncorrected estimate of %s",
                          "un")
        }
      }
  
      if(any(correction == "border" | correction == "bord.modif")) {
      # border method
      # Compute distances to boundary
        b <- bdist.points(X)
        bI <- b[II]
      # apply reduced sample algorithm
        RS <- Kount(DD, bI, b, breaks)
        if(any(correction == "bord.modif")) {
          # modified border correction
          denom.area <- eroded.areas(W, r)
          numTbm <- RS$numerator
          denTbm <- lambda3 * denom.area
          Tbm <- numTbm/denTbm
          TT <- bind.fv(TT, data.frame(bord.modif=Tbm), "hat(%s)[bordm](r)",
                        "modified border-corrected estimate of %s",
                        "bord.modif")
          if(ratio) {
            # save numerator and denominator
            numT <- bind.fv(numT, data.frame(bord.modif=numTbm),
                            "hat(%s)[bordm](r)",
                      "numerator of modified border-corrected estimate of %s",
                            "bord.modif")
            denT <- bind.fv(denT, data.frame(bord.modif=denTbm),
                            "hat(%s)[bordm](r)",
                      "denominator of modified border-corrected estimate of %s",
                            "bord.modif")
          }
        }
        if(any(correction == "border")) {
          numTb <- RS$numerator
          denTb <- lambda2 * RS$denom.count
          Tb <- numTb/denTb
          TT <- bind.fv(TT, data.frame(border=Tb), "hat(%s)[bord](r)",
                        "border-corrected estimate of %s",
                        "border")
          if(ratio) {
            numT <- bind.fv(numT, data.frame(border=numTb), "hat(%s)[bord](r)",
                            "numerator of border-corrected estimate of %s",
                            "border")
            denT <- bind.fv(denT, data.frame(border=denTb), "hat(%s)[bord](r)",
                            "denominator of border-corrected estimate of %s",
                            "border")
          }
        }
      }

      if(any(correction == "translate")) {
        # translation correction
        # apply to triangle list
        edgewt <- edgetri.Trans(X, tri[, 1:3])
        wh <- whist(tri$diam, breaks$val, edgewt)
        numTtrans <- 3 * cumsum(wh)
        denTtrans <- lambda3 * areaW
        Ttrans <- numTtrans/denTtrans
        h <- diameter(W)/2
        Ttrans[r >= h] <- NA
        TT <- bind.fv(TT, data.frame(trans=Ttrans), "hat(%s)[trans](r)",
                      "translation-corrected estimate of %s",
                      "trans")
        if(ratio) {
          numT <- bind.fv(numT, data.frame(trans=numTtrans),
                          "hat(%s)[trans](r)",
                          "numerator of translation-corrected estimate of %s",
                          "trans")
          denT <- bind.fv(denT, data.frame(trans=denTtrans),
                          "hat(%s)[trans](r)",
                          "denominator of translation-corrected estimate of %s",
                          "trans")
        }
      }
      # default plot will display all edge corrections
      formula(TT) <- . ~ r
      unitname(TT) <- unitname(X)
      #
      if(ratio) {
        # finish up numerator & denominator
        formula(numT) <- formula(denT) <- . ~ r
        unitname(numT) <- unitname(denT) <- unitname(TT)
        # tack on to result
        TT <- rat(TT, numT, denT, check=FALSE)
      }
      return(TT)
    }
  
  Tstat
})



