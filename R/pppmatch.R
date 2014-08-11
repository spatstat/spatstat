#
# pppmatch.R
#
# $Revision: 1.11 $  $Date: 2007/10/26 16:03:17 $
#
# Code by Dominic Schuhmacher
#
#
# -----------------------------------------------------------------
# The standard functions for the new class pppmatching
#
# Objects of class pppmatching consist of two point patterns pp1 and pp2,
# and either an adjacency matrix ((i,j)-th entry 1 if i-th point of pp1 and j-th
# point of pp2 are matched, 0 otherwise) for "full point matchings" or
# a "generalized adjacency matrix" (or flow matrix; positive values are
# no longer limited to 1, (i,j)-th entry gives the "flow" between
# the i-th point of pp1 and the j-th point of pp2) for "fractional matchings".
# Optional elements are the type
# of the matching, the cutoff value for distances in R^2, the order
# of averages taken, and the resulting distance for the matching.
# Currently recognized types are "spa" (subpattern assignment,
# where dummy points at maximal dist are introduced if cardinalities differ), 
# "ace" (assignment if cardinalities equal, where dist is maximal if cards differ),
# and "mat" (mass transfer, fractional matching that belongs to the
# Wasserstein distance obtained if point patterns are normalized to probability measures).
# -----------------------------------------------------------------

pppmatching <- function(X, Y, am, type = NULL, cutoff = NULL,
   q = NULL, mdist = NULL) {
   verifyclass(X, "ppp")
   verifyclass(Y, "ppp")
   n1 <- X$n
   n2 <- Y$n
   am <- as.matrix(am)
   am <- apply(am, c(1,2), as.numeric)
   if (length(am) == 0) {
      if (min(n1,n2) == 0) 
         am <- matrix(am, nrow=n1, ncol=n2)
      else
         stop("Adjacency matrix does not have the right dimensions")
   }
   if (dim(am)[1] != n1 || dim(am)[2] != n2)
      stop("Adjacency matrix does not have the right dimensions")
   res <- list("pp1" = X, "pp2" = Y, "matrix" = am, "type" = type, "cutoff" = cutoff, 
      "q" = q, "distance" = mdist)
   class(res) <- "pppmatching"
   res
}

# currently, for fractional matchings all the flows are plotted the same way
# irrespective of their weights
plot.pppmatching <- function(x, addmatch = NULL, main = NULL, ...) {
   if (is.null(main))
      main <- short.deparse(substitute(x))
   pp1 <- x$pp1
   pp2 <- x$pp2
   plot.owin(pp1$window, main = main, ...)
   here <- which((x$matrix > 0), arr.ind = TRUE)
   if (!is.null(addmatch)) {
      addhere <- which((addmatch > 0), arr.ind = TRUE)
      seg <- as.psp(from=pp1[addhere[,1]], to=pp2[addhere[,2]])
      plot(seg, add=TRUE, lty = 2, col="gray70")
   }
   if (length(here) > 0) {
     seg <- as.psp(from=pp1[here[,1]], to=pp2[here[,2]])
     plot(seg, add=TRUE, ...)
   }
   points(x$pp1, pch=20, col=2, ...)
   points(x$pp2, pch=20, col=4, ...)
   return(invisible(NULL))
}

print.pppmatching <- function(x, ...) {
   n1 <- x$pp1$n
   n2 <- x$pp2$n
   if (is.null(x$type) || is.null(x$q) || is.null(x$cutoff))
     cat("Generic matching of two planar point patterns \n")
   else
     cat(x$type, "-", x$q, " matching of two planar point patterns (cutoff = ",
       x$cutoff, ") \n", sep = "")
   cat("pp1:", n1, ngettext(n1, "point", "points"), "\n")
   cat("pp2:", n2, ngettext(n2, "point", "points"), "\n")
   print.owin(x$pp1$window)
   npair <- sum(x$matrix > 0)
   if (npair == 0)
     cat("matching is empty \n") 
   else {
     if (any(x$matrix != trunc(x$matrix)))
       cat("fractional matching,", npair, ngettext(npair, "flow", "flows"), "\n")
     else
       cat("point matching,", npair, ngettext(npair, "line", "lines"), "\n")
   }
   if (!is.null(x$distance))
     cat("distance:", x$distance, "\n") 
   return(invisible(NULL))
}

summary.pppmatching <- function(object, ...) {
   X <- object$pp1
   Y <- object$pp2
   n1 <- X$n
   n2 <- Y$n
   if (is.null(object$type) || is.null(object$q) || is.null(object$cutoff))
     cat("Generic matching of two planar point patterns \n")
   else
     cat(object$type, "-", object$q, " matching of two planar point patterns (cutoff = ",
       object$cutoff, ") \n", sep = "")
   cat("pp1:", n1, ngettext(n1, "point", "points"), "\n")
   cat("pp2:", n2, ngettext(n2, "point", "points"), "\n")
   print.owin(X$window)
   npair <- sum(object$matrix > 0)
   if (npair == 0)
     cat("matching is empty \n") 
   else {
     if (any(object$matrix != trunc(object$matrix))) {
       cat("fractional matching,", npair, ngettext(npair, "flow", "flows"), "\n")
     }
     else {
       cat("point matching,", npair, ngettext(npair, "line", "lines"), "\n")
       rowsum <- apply(object$matrix, 1, "sum")
       colsum <- apply(object$matrix, 2, "sum")
       lt <- ifelse(min(rowsum) >= 1, TRUE, FALSE)
       ru <- ifelse(max(rowsum) <= 1, TRUE, FALSE)
       rt <- ifelse(min(colsum) >= 1, TRUE, FALSE)
       lu <- ifelse(max(colsum) <= 1, TRUE, FALSE)
       if (lt && ru && rt && lu)
         cat("matching is 1-1 \n")
       else if (any(lt, ru, rt, lu)) {
         cat("matching is",
                   ifelse(lt, " left-total", ""),
                   ifelse(lu, " left-unique", ""),
                   ifelse(rt, " right-total", ""),
                   ifelse(ru, " right-unique", ""),
                   "\n", sep="")
         }
     }
   }
   if (!is.null(object$distance))
     cat("distance:", object$distance, "\n") 
   return(invisible(NULL))
}


# -----------------------------------------------------------------
# matchingdist computes the distance associated with a certain kind of matching.
# Any of the arguments type, cutoff and order (if supplied) override the 
# the corresponding arguments in the matching.
# This function is useful for verifying the distance element of an
# object of class pppmatching as well as for comparing different
# (typically non-optimal) matchings.
# -----------------------------------------------------------------

matchingdist <- function(matching, type = NULL, cutoff = NULL, q = NULL) {
  verifyclass(matching, "pppmatching")
  if (is.null(type))
    if (is.null(matching$type))
      stop("Type of matching unknown. Distance cannot be computed")
    else
      type <- matching$type
  if (is.null(cutoff))
    if (is.null(matching$cutoff))
      stop("Cutoff value unknown. Distance cannot be computed")
    else
      cutoff <- matching$cutoff
  if (is.null(q))
    if (is.null(matching$q))
      stop("Order unknown. Distance cannot be computed")
    else
      q <- matching$q

  X <- matching$pp1
  Y <- matching$pp2
  n1 <- X$n
  n2 <- Y$n
  Lpexpect <- function(x, w, p) {
    f <- max(x)
      return(ifelse(f==0, 0, f * sum((x/f)^p * w)^(1/p)))
  }

  if (type == "spa") {
    n <- max(n1,n2) # divisor for Lpexpect
    if (n == 0)
      return(0)
    else if (min(n1,n2) == 0)
      return(cutoff)
    shortdim <- which.min(c(n1,n2))
    shortsum <- apply(matching$matrix, shortdim, sum)
    if (any(shortsum != 1))
      warning("matching does not attribute mass 1 to each point of point pattern with smaller cardinality")
    dfix <- apply(crossdist(X,Y), c(1,2), function(x) { min(x,cutoff) })
    if (is.finite(q))
      resdist <- (Lpexpect(dfix, matching$matrix/n, q)^q + abs(n2-n1)/n * cutoff^q)^(1/q)
    else
      resdist <- ifelse(n1==n2, max(dfix[matching$matrix > 0]), cutoff)
  }
  else if (type == "ace") {
    n <- n1 # divisor for Lpexpect
    if (n1 != n2)
      return(cutoff)
    if (n == 0)
      return(0)
    rowsum <- apply(matching$matrix, 1, sum)
    colsum <- apply(matching$matrix, 2, sum)
    if (any(c(rowsum, colsum) != 1))
      warning("matching is not 1-1")
    dfix <- apply(crossdist(X,Y), c(1,2), function(x) { min(x,cutoff) })
    if (is.finite(q))
      resdist <- Lpexpect(dfix, matching$matrix/n, q)
    else
      resdist <- max(dfix[matching$matrix > 0])
  }
  else if (type == "mat") {
    n <- min(n1,n2) # divisor for Lpexpect
    if (min(n1,n2) == 0)
      return(NaN)
    shortdim <- which.min(c(n1,n2))
    shortsum <- apply(matching$matrix, shortdim, sum)
    if (any(shortsum != 1))
      warning("matching does not attribute mass 1 to each point of point pattern with smaller cardinality")
    dfix <- apply(crossdist(X,Y), c(1,2), function(x) { min(x,cutoff) })
    if (is.finite(q))
      resdist <- Lpexpect(dfix, matching$matrix/n, q)
    else
      resdist <- max(dfix[matching$matrix > 0])
  }
  else 
    stop(paste("Unrecognised type", sQuote(type)))
  return(resdist)
}


# -----------------------------------------------------------------
# The main function for computation of distances and finding optimal
# matchings between point patterns: pppdist
# -----------------------------------------------------------------
#
# pppdist uses several helper functions not normally called by the user 
#
# The arguments of pppdist are 
#
# x and y of class ppp (the two point patterns for which we want to compute
#   a distance)
# The type of distance to be computed; any one of "spa" (default), "ace", "mat".
#   For details of this and the following two arguments see above (description
#   for class "pppmatching")
# cutoff and order q of the distance
# Set matching to TRUE if the full point matching (including distance)
#   should be returned; otherwise only the distance is returned
# If ccode is FALSE R code is used where available. This may be useful if q
#   is high (say above 10) and severe warning messages pop up. R can
#   (on most machines) deal with a higher number of significant digits per
#   number than C (at least with the code used below)
# precision should only be entered by advanced users. Empirically reasonable defaults
#   are used otherwise. As a rule of thumb, if ccode is TRUE, precision should
#   be the highest value that does not give an error (typically 9); if ccode
#   is FALSE, precision should be balanced (typically between 10 and 100) in
#   such a way that the sum of the  number of zeroes and pseudo-zeroes given in the
#   warning messages is minimal
# approximation: if q = Inf, by the distance of which order should 
#   the true distance be approximated. If approximation is Inf, brute force
#   computation is used, which is only practicable for point patterns with
#   very few points (see also the remarks just before the pppdist.prohorov
#   function below).  
# show.rprimal=TRUE shows at each stage of the algorithm what the current restricted
#   primal problem and its solution are (algorithm jumps between restricted primal
#   and dual problem until the solution to the restricted primal (a partial
#   matching of the point patterns) is a full matching)
# timelag gives the number of seconds of pause added each time a solution to
#   the current restricted primal is found (has only an effect if show.primal=TRUE) 
# -----------------------------------------------------------------

pppdist <- function(X, Y, type = "spa", cutoff = 1, q = 1, matching = TRUE,
  ccode = TRUE, precision = NULL, approximation = 10, show.rprimal = FALSE, timelag = 0) {

  verifyclass(X, "ppp")
  verifyclass(Y, "ppp")
  if (!ccode && type == "mat") {
    warning("R code is not available for type = ", dQuote("mat"), ". C code is
    used instead")
    ccode <- TRUE
  }
  if (!ccode && is.infinite(q) && is.infinite(approximation)) {
    warning("R code is not available for q = Inf and approximation = Inf. C code is
    used instead")
    ccode <- TRUE
  }
  if (ccode && is.infinite(q) && is.infinite(approximation) && type == "spa" && X$n != Y$n) {
    warning("approximation = Inf not available for type = ",
        dQuote("spa"), " and point patterns with differing cardinalities")
    approximation <- 10
  }
  if (is.infinite(q) && is.infinite(approximation) && type == "mat") {
    warning("approximation = Inf not available for type = ",
        dQuote("mat"))
    approximation <- 10
  }
  if (show.rprimal) {
    ccode <- FALSE
      if (type != "ace"){
        warning("show.rprimal = TRUE not available for type = ",
        dQuote(type), ". Type is changed to ", dQuote("ace"))
        type <- "ace"
    }
  }

  if (is.null(precision)) {
    if (ccode)
      precision <- trunc(log10(.Machine$integer.max))
    else {
      db <- .Machine$double.base
      minprec <- trunc(log10(.Machine$double.base^.Machine$double.digits))
      if (is.finite(q))
        precision <- min(max(minprec,2*q),(.Machine$double.max.exp-1)*log(db)/log(10))
      else
        precision <- min(max(minprec,2*approximation),(.Machine$double.max.exp-1)*log(db)/log(10))
      }
  }

  if (type == "spa") {
    if (X$n == 0 && Y$n == 0) {
      if (!matching)
        return(0)
      else {
        return(pppmatching(X, Y, matrix(0, nrow=0,ncol=0), type, cutoff, q, 0))
      }
    }
    n1 <- X$n
    n2 <- Y$n
    n <- max(n1,n2)
    dfix <- matrix(cutoff,n,n)
    if (min(n1,n2) > 0)
      dfix[1:n1,1:n2] <- crossdist(X,Y)
    d <- dfix <- apply(dfix, c(1,2), function(x) { min(x,cutoff) })
    if (is.infinite(q)) {
      if (n1 == n2 || matching)
        return(pppdist.prohorov(X, Y, n, d, type, cutoff, matching, ccode,
        precision, approximation))
      else
        return(cutoff)
      # in the case n1 != n2 the distance is clear, and in a sense any
      # matching would be correct. We go here the extra mile and call
      # pppdist.prohorov in order to find (approximate) the matching
      # that is intuitively most interesting (i.e. the one that
      # pairs the points of the
      # smaller cardinality point pattern with the points of the larger
      # cardinality point pattern in such a way that the maximal pairing distance
      # is minimal (for q < Inf the q-th order pairing distance before the introduction
      # of dummy points is automatically minimal if it is minimal after the
      # introduction of dummy points)
      # which would be the case for the obtained pairing if q < Inf
    }
  }
  else if (type == "ace") {
    if (X$n != Y$n) {
      if (!matching)
        return(cutoff)
      else {
        return(pppmatching(X, Y, matrix(0, nrow=X$n, ncol=Y$n), type, cutoff, q, cutoff))
      }
    }
    if (X$n == 0) {
      if (!matching)
        return(0)
      else {
        return(pppmatching(X, Y, matrix(0, nrow=0,ncol=0), type, cutoff, q, 0))
      }
    }
    n <- n1 <- n2 <- X$n
    dfix <- crossdist(X,Y)
    d <- dfix <- apply(dfix, c(1,2), function(x) { min(x,cutoff) })
    if (is.infinite(q))
      return(pppdist.prohorov(X, Y, n, d, type, cutoff, matching, ccode,
      precision, approximation))
  }
  else if (type == "mat") {
    if (!ccode)
      warning("R code is not available for type = ", dQuote("mat"), ". C code is used instead")
    return(pppdist.mat(X, Y, cutoff, q, matching, precision, approximation))
  }
  else stop(paste("Unrecognised type", sQuote(type)))

  d <- d/max(d)
  d <- round((d^q)*(10^precision))
  nzeroes <- sum(d == 0 & dfix > 0)
  if(nzeroes > 0)
    warning(paste(nzeroes, ngettext(nzeroes, "zero", "zeroes"), "introduced, while rounding the q-th powers of distances"))
  if(ccode & any(d > .Machine$integer.max))
    stop("integer overflow, while rounding the q-th powers of distances")
  if(!ccode) {
    if (any(is.infinite(d)))
      stop("Inf obtained, while taking the q-th powers of distances")
    maxd <- max(d)
    npszeroes <- sum(maxd/d[d>0] >= .Machine$double.base^.Machine$double.digits)
    if (npszeroes > 0)
      warning(paste(npszeroes, ngettext(npszeroes, "pseudo-zero", "pseudo-zeroes"), "introduced, while taking the q-th powers of distances"))
      # a pseudo-zero is a value that is positive but contributes nothing to the
      # q-th order average because it is too small compared to the other values
  }

  Lpmean <- function(x, p) {
    f <- max(x)
    return(ifelse(f==0, 0, f * mean((x/f)^p)^(1/p)))
  }
    
  if (show.rprimal && type == "ace") {
    assig <- acedist.show(X, Y, n, d, timelag)
    am <- matrix(0, n, n)
    am[cbind(1:n, assig[1:n])] <- 1
  }
  else if (ccode) {
    res <- .C("dwpure",
             as.integer(d),
             as.integer(rep(1,n)),
             as.integer(rep(1,n)),
             as.integer(n),
             as.integer(n),
             flowmatrix = as.integer(rep(0,n^2)),
             PACKAGE="spatstat")
    am <- matrix(res$flowmatrix, n, n)
  }
  else {
    assig <- acedist.noshow(X, Y, n, d)
    am <- matrix(0, n, n)
    am[cbind(1:n, assig[1:n])] <- 1
  }
  resdist <- Lpmean(dfix[am == 1], q)
  if (!matching)
    return(resdist)
  else {
    amsmall <- suppressWarnings(matrix(am[1:n1,1:n2], nrow=n1, ncol=n2))
    # previous line solves various problems associated with min(n1,n2) = 0 or = 1
    return(pppmatching(X, Y, amsmall, type, cutoff, q, resdist))
  }
}   

#
#
# ===========================================================
# ===========================================================
#                   Anything below:
#    Internal functions usually not to be called by user
# ===========================================================
# ===========================================================
#

#
#   Called if show.rprimal is true
#

acedist.show <- function(X, Y, n, d, timelag = 0) {
      plot(pppmatching(X, Y, matrix(0, n, n)))
      # initialization of dual variables
      u <- apply(d, 1, min)
      d <- d - u
      v <- apply(d, 2, min)
      d <- d - rep(v, each=n)
      # the main loop
      feasible <- FALSE
      while (!feasible) {
         rpsol <- maxflow(d)  # rpsol = restricted primal, solution
         am <- matrix(0, n, n)
         for (i in 1:n) {
            if (rpsol$assignment[i] > -1) am[i, rpsol$assignment[i]] <- TRUE
         }
         Sys.sleep(timelag)
         channelmat <- (d == 0 & !am)
         plot(pppmatching(X, Y, am), addmatch = channelmat)
         # if the solution of the restricted primal is not feasible for  
         # the original primal, update dual variables
         if (min(rpsol$assignment) == -1) {
            w1 <- which(rpsol$fi_rowlab > -1)
            w2 <- which(rpsol$fi_collab == -1)
            subtractor <- min(d[w1, w2])
            d[w1,] <- d[w1,] - subtractor
            d[,-w2] <- d[,-w2] + subtractor 
         }
         # otherwise break the loop
         else {
            feasible <- TRUE
         }   
      }
      return(rpsol$assignment)
}

#
#   R-version of hungarian algo without the pictures
#   useful if q is large
#

acedist.noshow <- function(X, Y, n, d) {
      # initialization of dual variables
      u <- apply(d, 1, min)
      d <- d - u
      v <- apply(d, 2, min)
      d <- d - rep(v, each=n)
      # the main loop
      feasible <- FALSE
      while (!feasible) {
         rpsol <- maxflow(d)  # rpsol = restricted primal, solution
         am <- matrix(0, n, n)
         for (i in 1:n) {
            if (rpsol$assignment[i] > -1) am[i, rpsol$assignment[i]] <- TRUE
         }
         channelmat <- (d == 0 & !am)
         # if the solution of the restricted primal is not feasible for  
         # the original primal, update dual variables
         if (min(rpsol$assignment) == -1) {
            w1 <- which(rpsol$fi_rowlab > -1)
            w2 <- which(rpsol$fi_collab == -1)
            subtractor <- min(d[w1, w2])
            d[w1,] <- d[w1,] - subtractor
            d[,-w2] <- d[,-w2] + subtractor 
         }
         # otherwise break the loop
         else {
            feasible <- TRUE
         }   
      }
      return(rpsol$assignment)
}

#  
# Solution of restricted primal
# 

maxflow <- function(costm) {
  stopifnot(is.matrix(costm))
  stopifnot(nrow(costm) == ncol(costm))
  if(!all(apply(costm == 0, 1, any)))
    stop("Each row of the cost matrix must contain a zero")
  
  m <- dim(costm)[1]   # cost matrix is square m * m
  assignment <- rep(-1, m)   # -1 means no pp2-point assigned to i-th pp1-point
   # initial assignment or rowlabel <- source label (= 0) where not possible
   for (i in 1:m) {
      j <- match(0, costm[i,])
      if (!(j %in% assignment))
         assignment[i] <- j
   }
   newlabelfound <- TRUE
   while (newlabelfound) {
     rowlab <- rep(-1, m)   # -1 means no label given, 0 stands for source label
     collab <- rep(-1, m)
     rowlab <- ifelse(assignment == -1, 0, rowlab)
     # column and row labeling procedure until either breakthrough occurs
     # (which means that there is a better point assignment, i.e. one that
     # creates more point pairs than the current one (flow can be increased))
     # or no more labeling is possible
     breakthrough <- -1
     while (newlabelfound && breakthrough == -1) { 
         newlabelfound <- FALSE
         for (i in 1:m) {
            if (rowlab[i] != -1) {
               for (j in 1:m) {
                  if (costm[i,j] == 0 && collab[j] == -1) {
                     collab[j] <- i
                     newlabelfound <- TRUE
                     if (!(j %in% assignment) && breakthrough == -1)
                        breakthrough <- j
                  }
               }
            }
         }
         for (j in 1:m) {
            if (collab[j] != -1) {
               for (i in 1:m) {
                  if (assignment[i] == j && rowlab[i] == -1) {
                     rowlab[i] <- j
                     newlabelfound <- TRUE
                  }
               }
            }
         }
      }
      # if the while-loop was left due to breakthrough,
      # reassign points (i.e. redirect flow) and restart labeling procedure
      if (breakthrough != -1) {
         l <- breakthrough
         while (l != 0) {
            k <- collab[l]
            assignment[k] <- l
            l <- rowlab[k] 
         }
      }
   }
   # the outermost while-loop is left, no more labels can be given; hence
   # the maximal number of points are paired given the current restriction
   # (flow is maximal given the current graph)
   return(list("assignment"=assignment, "fi_rowlab"=rowlab, "fi_collab"=collab))  
}

# 
# Prohorov distance computation/approximation (called if q = Inf in pppdist
#   and type = "spa" or "ace")
# Exact brute force computation of distance if approximation = Inf,
#   scales very badly, should not be used for cardinality n larger than 10-12
# Approximation by order q distance gives often (if the warning messages 
#   are not too extreme) the right matching and therefore the exact Prohorov distance,
#   but in very rare cases the result can be very wrong. However, it is always
#   an exact upper bound of the Prohorov distance (since based on *a* pairing
#   as opposed to optimal pairing.
#

pppdist.prohorov <- function(X, Y, n, dfix, type, cutoff = 1, matching = TRUE,
  ccode = TRUE, precision = 9, approximation = 10) {
  n1 <- X$n
  n2 <- Y$n
  d <- dfix/max(dfix)
  if (is.finite(approximation)) {
      warning(paste("distance with parameter q = Inf is approximated by distance with parameter q =", approximation))
    d <- round((d^approximation)*(10^precision)) 
    nzeroes <- sum(d == 0 & dfix > 0)
    if (nzeroes > 0)
      warning(paste(nzeroes, ngettext(nzeroes, "zero", "zeroes"), "introduced, while rounding distances"))
    if (ccode) {
      if (any(d > .Machine$integer.max))
        stop("integer overflow, while rounding the q-th powers of distances")
      res <- .C("dwpure",
               as.integer(d),
               as.integer(rep(1,n)),
               as.integer(rep(1,n)),
               as.integer(n),
               as.integer(n),
               flowmatrix = as.integer(rep(0,n^2)),
               PACKAGE="spatstat")
      am <- matrix(res$flowmatrix, n, n)
    }
    else {
      if (any(is.infinite(d)))
        stop("Inf obtained, while taking the q-th powers of distances")
      maxd <- max(d)
      npszeroes <- sum(maxd/d[d>0] >= .Machine$double.base^.Machine$double.digits)
      if (npszeroes > 0)
        warning(paste(npszeroes, ngettext(npszeroes, "pseudo-zero", "pseudo-zeroes"), "introduced, while taking the q-th powers of distances"))
      assig <- acedist.noshow(X, Y, n, d)
      am <- matrix(0, n, n)
      am[cbind(1:n, assig[1:n])] <- 1
    }
  }
  else {
    d <- round(d*(10^precision))
    nzeroes <- sum(d == 0 & dfix > 0)
    if (nzeroes > 0)
      warning(paste(nzeroes, ngettext(nzeroes, "zero", "zeroes"), "introduced, while rounding distances"))
    if (any(d > .Machine$integer.max))
      stop("integer overflow, while rounding the q-th powers of distances")
    res <- .C("dinfty_R",
             as.integer(d),
             as.integer(n),
             assignment = as.integer(rep(-1,n)),
             PACKAGE="spatstat")
    assig <- res$assignment
    am <- matrix(0, n, n)
    am[cbind(1:n, assig[1:n])] <- 1
  }
  if (n1 == n2)
    resdist <- max(dfix[am == 1])
  else
    resdist <- cutoff
  if (!matching)
    return(resdist)
  else {
    amsmall <- suppressWarnings(matrix(am[1:n1,1:n2], nrow=n1, ncol=n2))
    # previous line solves various problems associated with min(n1,n2) = 0 or = 1
    return(pppmatching(X, Y, amsmall, type, cutoff, Inf, resdist))
  }
}   

# 
# Computation of "pure Wasserstein distance" for any q (called if type="mat"
#   in pppdist, no matter if q finite or not).
# If q = Inf, approximation using ccode is enforced
# (approximation == Inf is not allowed here).
#

pppdist.mat <- function(X, Y, cutoff = 1, q = 1, matching = TRUE, precision = 9,
  approximation = 10) {
  n1 <- X$n
  n2 <- Y$n
  n <- min(n1,n2)
  if (n == 0) {
    if (!matching)
      return(NaN)
    else
      return(pppmatching(X, Y, matrix(0, nrow=0,ncol=0), "mat", cutoff, q, NaN))
  }

  dfix <- crossdist(X,Y)
  d <- dfix <- apply(dfix, c(1,2), function(x) { min(x,cutoff) })
  d <- d/max(d)
  if (is.infinite(q)) {
    if (is.infinite(approximation))
      stop("approximation = Inf")
    warning(paste("distance with parameter q = Inf is approximated by distance with parameter q =", approximation))
    d <- round((d^approximation)*(10^precision)) 
    nzeroes <- sum(d == 0 & dfix > 0)
    if (nzeroes > 0)
      warning(paste(nzeroes, "zeroes introduced, while rounding distances"))
    if (any(d > .Machine$integer.max))
      stop("integer overflow, while rounding the q-th powers of distances")
    gcd <- greatest.common.divisor(n1,n2)
    mass1 <- n2/gcd
    mass2 <- n1/gcd

    res <- .C("dwpure",
             as.integer(d),
             as.integer(rep(mass1,n1)),
             as.integer(rep(mass2,n2)),
             as.integer(n1),
             as.integer(n2),
             flowmatrix = as.integer(rep(0,n1*n2)),
             PACKAGE="spatstat")
    am <- matrix(res$flowmatrix/(max(n1,n2)/gcd), n1, n2)
    resdist <- max(dfix[am > 0])
  }
  else {
    d <- round((d^q)*(10^precision))
    nzeroes <- sum(d == 0 & dfix > 0)
    if(nzeroes > 0)
      warning(paste(nzeroes, ngettext(nzeroes, "zero", "zeroes"), "introduced, while rounding the q-th powers of distances"))
    if(any(d > .Machine$integer.max))
      stop("integer overflow, while rounding the q-th powers of distances")
    gcd <- greatest.common.divisor(n1,n2)
    mass1 <- n2/gcd
    mass2 <- n1/gcd

    Lpexpect <- function(x, w, p) {
      f <- max(x)
      return(ifelse(f==0, 0, f * sum((x/f)^p * w)^(1/p)))
    }

    res <- .C("dwpure",
             as.integer(d),
             as.integer(rep(mass1,n1)),
             as.integer(rep(mass2,n2)),
             as.integer(n1),
             as.integer(n2),
             flowmatrix = as.integer(rep(0,n1*n2)),
             PACKAGE="spatstat")
    am <- matrix(res$flowmatrix/(max(n1,n2)/gcd), n1, n2)
    # our "adjacency matrix" in this case is standardized to have
    # rowsum 1 if n1 <= n2 and colsum 1 if n1 >= n2
    resdist <- Lpexpect(dfix, am/n, q)
  }
  if (!matching)
    return(resdist)
  else {
   amsmall <- suppressWarnings(matrix(am[1:n1,1:n2], nrow=n1, ncol=n2))
   # previous line solves various problems associated with min(n1,n2) = 0 or = 1
   return(pppmatching(X, Y, amsmall, "mat", cutoff, q, resdist))
  }
}

