#
# connected.R
#
# connected component transform
#
#    $Revision: 1.16 $  $Date: 2013/08/23 07:34:24 $
#
# Interpreted code for pixel images by Julian Burgos <jmburgos@u.washington.edu>
# Rewritten in C by Adrian Baddeley
#
# Code for point patterns by Adrian Baddeley

connected <- function(X, ...) {
  UseMethod("connected")
}

connected.im <- function(X, ..., background=NA, method="C") {
  method <- pickoption("algorithm choice", method,
                       c(C="C", interpreted="interpreted"))
  if(!is.na(background))
    X <- solutionset(X != background)
  else
    X <- as.mask(as.owin(X))
  connected.owin(X, method=method)
}

connected.owin <- function(X, ..., method="C") {
  method <- pickoption("algorithm choice", method,
                       c(C="C", interpreted="interpreted"))
  # convert X to binary mask
  X <- as.mask(X)
  #     
  Y <- X$m
  nr <- X$dim[1]
  nc <- X$dim[2]

  if(method == "C") {
################ COMPILED CODE #########################
# Pad border with FALSE
    M <- rbind(FALSE, Y, FALSE)
    M <- cbind(FALSE, M, FALSE)
    # assign unique label to each foreground pixel 
    L <- M
    L[M] <- seq_len(sum(M))
    L[!M] <- 0
    # resolve labels
    #        (warning: 'mat' is overwritten, so DUP=TRUE)
    z <- .C("cocoImage",
            mat=as.integer(t(L)),
            nr=as.integer(nr),
            nc=as.integer(nc))
#            PACKAGE="spatstat")
    # unpack
    Z <- matrix(z$mat, nr+2, nc+2, byrow=TRUE)
  } else {
################ INTERPRETED CODE #########################
# by Julian Burgos
#  
# Pad border with zeros
    padY <- rbind(0, Y, 0)
    padY <- cbind(0, padY, 0)
    # Initialise 
    Z <- matrix(0, nrow(padY), ncol(padY))
    currentlab <- 1
    todo <- as.vector(t(Y))
    equiv <- NULL
  
    # ........ main loop ..........................
    while(any(todo)){
      # pick first unresolved pixel
      one <- which(todo)[1]
      onerow <- ceiling(one/nc)
      onecol <- one -((onerow-1)*nc)
      parow=onerow+1 # Equivalent rows & column in padded matrix
      pacol=onecol+1
      #Examine four previously scanned neighbors
      # (use padded matrix to avoid edge issues)
      nbrs <- rbind(c(parow-1,pacol-1),
                    c(parow-1,pacol),
                    c(parow,  pacol-1),
                    c(parow-1,pacol+1))
      px <- sum(padY[nbrs])
      if (px==0){
        # no neighbours: new component
        Z[parow,pacol] <- currentlab
        currentlab <- currentlab+1
        todo[one] <- FALSE
      } else if(px==1) {
        # one neighbour: assign existing label
        labs <- unique(Z[nbrs], na.rm=TRUE)
        labs <- labs[labs != 0]
        Z[parow,pacol] <- labs[1]
        currentlab <- max(Z)+1
        todo[one] <- FALSE
      } else {
        # more than one neighbour: possible merger of labels
        labs <- unique(Z[nbrs], na.rm=TRUE)
        labs <- labs[labs != 0]
        labs <- sort(labs)
        equiv <- rbind(equiv,c(labs,rep.int(0,times=4-length(labs))))
        Z[parow,pacol] <- labs[1]
        currentlab <- max(Z)+1
        todo[one] <- FALSE
      }
    }
    # ........... end of loop ............
    # Resolve equivalences ................

    if(length(equiv)>1){
      merges <- (equiv[,2] > 1)
      nmerge <- sum(merges)
      if(nmerge==1)
        equiv <- equiv[which(merges), , drop=FALSE]
      else if(nmerge > 1) {
        relevant <- (equiv[,2] > 0)
        equiv <- equiv[relevant, , drop=FALSE]
        equiv <- equiv[fave.order(equiv[,1]),]
      }
      for (i in 1:nrow(equiv)){
        current <- equiv[i, 1]
        for (j in 2:4){
          twin <- equiv[i,j]
          if (twin>0){
            # Change labels matrix
            Z[which(Z==twin)] <- current
            # Update equivalence table
            equiv[which(equiv==twin)] <- current
          }
        }
      }
    }
  }

  ########### COMMON CODE ############################
    
  # Renumber labels sequentially
  mapped <- (Z != 0)
  usedlabs <- sort(unique(as.vector(Z[mapped])))
  nlabs <- length(usedlabs)
  labtable <- cumsum(seq_len(max(usedlabs)) %in% usedlabs)
  Z[mapped] <- labtable[Z[mapped]]

  # banish zeroes
  Z[!mapped] <- NA
  
  # strip borders
  Z <- Z[2:(nrow(Z)-1),2:(ncol(Z)-1)]
  # dress up 
  Z <- im(factor(Z, levels=1:nlabs),
          xcol=X$xcol, yrow=X$yrow, unitname=unitname(X))
  return(Z)
}

connected.ppp <- function(X, R, ...) {
  stopifnot(is.ppp(X))
  check.1.real(R, "In connected.ppp")
  stopifnot(R >= 0)
  internal <- resolve.1.default("internal", list(...), list(internal=FALSE))
  nv <- npoints(X)
  cl <- closepairs(X, R, what="indices")
  ie <- cl$i - 1L
  je <- cl$j - 1L
  ne <- length(ie)
  DUP <- spatstat.options("dupC")
  zz <- .C("cocoGraph",
           nv=as.integer(nv),
           ne=as.integer(ne),
           ie=as.integer(ie),
           je=as.integer(je),
           label=as.integer(integer(nv)),
           status=as.integer(integer(1)),
           DUP=DUP)
#           PACKAGE="spatstat")
  if(zz$status != 0)
    stop("Internal error: connected.ppp did not converge")
  if(internal)
    return(zz$label)
  lab <- zz$label + 1L
  # Renumber labels sequentially 
  lab <- as.integer(factor(lab))
  # Convert labels to factor
  lab <- factor(lab)
  # Apply to points
  Y <- X %mark% lab
  return(Y)
}

