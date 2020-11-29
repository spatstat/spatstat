#'
#' metriccontact.R
#'
#' Metric contact distribution
#' (corresponding distance transforms are defined in 'metricPdt.R')
#'
#'     $Revision: 1.1 $  $Date: 2020/11/29 07:41:37 $

rectcontact <- function(X, ..., asp=1.0, npasses=4,
                        eps=NULL, r=NULL, breaks=NULL,
                        correction=c("rs", "km")) {
  verifyclass(X, "im")
  rorbgiven <- !is.null(r) || !is.null(breaks)
  checkspacing <- !isFALSE(list(...)$checkspacing)
  testme       <- isTRUE(list(...)$testme)
  
  check.1.real(asp)
  stopifnot(asp > 0)
  
  if(X$type != "logical")
    stop("X should be a logical-valued image")

  if(!missing(eps))
    X <- as.im(X, eps=eps)
  
  W <- as.mask(X)      # the region that is defined
  Y <- solutionset(X)  # the region that is TRUE
  fullframe  <- all(W$m)
  emptyframe <- !any(W$m)
  
  ## histogram breakpoints
  rmaxdefault <- rmax.rule("F", W)
  breaks <- handle.r.b.args(r, breaks, W, eps, rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max

  if(testme || (rorbgiven && checkspacing))
    check.finespacing(rvals,
                      if(is.null(eps)) NULL else eps/4,
                      W,
                      rmaxdefault=if(rorbgiven) NULL else rmaxdefault,
                      action="fatal",
                      rname="r", 
                      context="in rectcontact(X, r)")
                                
  correction <- pickoption("correction", correction,
                           c(border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             best="km"),
                           multi=TRUE)
  
  ##  compute distances and censoring distances
  if(!emptyframe) {
    dist <- rectdistmap(Y, asp, npasses=npasses)
    if(fullframe) {
      bdry <- attr(dist, "bdist")
    } else {
      bdry <- rectdistmap(complement.owin(W), asp, npasses=npasses)
    }
    #' extract corresponding values
    dist <- dist[W, drop=TRUE, rescue=FALSE]
    bdry <- bdry[W, drop=TRUE, rescue=FALSE]
    ## censoring indicators
    d <- (dist <= bdry)
    ##  observed distances
    o <- pmin.int(dist, bdry)
  }

  ## calculate Kaplan-Meier and/or border corrected (Reduced Sample) estimators
  want.rs <- "rs" %in% correction
  want.km <- "km" %in% correction
  selection <- c(want.rs, want.km)
  tags <- c("rs", "km")[selection]
  labels <- c("hat(%s)[bord](r)", "hat(%s)[km](r)")[selection]
  descr <- c("border corrected estimate of %s",
             "Kaplan-Meier estimate of %s")[selection]
  if(emptyframe) {
    df <- as.data.frame(matrix(0, length(rvals), length(tags)))
    names(df) <- tags
  } else {
    df  <- km.rs.opt(o, bdry, d, breaks, KM=want.km, RS=want.rs)
    df <- as.data.frame(df[tags])
  }
  ## create fv object
  df <- cbind(data.frame(r=rvals), df)
  Z <- fv(df, "r", quote(H(r)),
          if(want.km) "km" else "rs",
          . ~ r,
          c(0,rmax),
          c("r", labels),
          c("distance argument r", descr),
          fname="H")

  fvnames(Z, ".") <- rev(fvnames(Z, "."))
  attr(Z, "alim") <- with(Z, range(.x[is.finite(.y) & .y <= 0.95]))
  attr(Z, "conserve") <- list(checkspacing=FALSE)
  return(Z)
}

