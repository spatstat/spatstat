#'  perspective plot of 3D 
#'
#'  $Revision: 1.6 $ $Date: 2018/04/30 09:31:33 $
#'


project3Dhom <- local({

  check3dvector <- function(x) {
    xname <- deparse(substitute(x))
    if(!(is.numeric(x) && length(x) == 3))
      stop(paste(xname, "should be a numeric vector of length 3"),
           call.=FALSE)
    return(NULL)
  }

  normalise <- function(x) {
    len <- sqrt(sum(x^2))
    if(len == 0) stop("Attempted to normalise a vector of length 0")
    return(x/len)
  }

  innerprod <- function(a, b) sum(a*b)

  crossprod <- function(u, v) {
    c(u[2] * v[3] - u[3] * v[2],
      -(u[1] * v[3] - u[3] * v[1]),
      u[1] * v[2] - u[2] * v[1])
  }

  project3Dhom <- function(xyz, eye=c(0,-3,1), org=c(0,0,0), vert=c(0,0,1)) {
    ## xyz: data to be projected (matrix n * 3)
    stopifnot(is.matrix(xyz) && ncol(xyz) == 3)
    ## eye: eye position (x,y,z)
    check3dvector(eye)
    ## org: origin (x,y,z) becomes middle of projection plane
    check3dvector(org)
    ## vert: unit vector in direction to become the 'vertical'
    if(!missing(vert)) {
      check3dvector(vert)
      vert <- normalise(vert)
    }
    ## vector pointing into screen
    vin <- normalise(org - eye)
    ## projection of vertical onto screen
    vup <- normalise(vert - innerprod(vert, vin) * vin)
    ## horizontal axis in screen
    vhoriz <- crossprod(vin, vup)
    ##
    dbg <- FALSE
    if(dbg) {
      cat("vin=")
      print(vin)
      cat("vup=")
      print(vup)
      cat("vhoriz=")
      print(vhoriz)
    }
    ## homogeneous coordinates
    hom <- t(t(xyz) - eye) %*% cbind(vhoriz, vup, vin)
    colnames(hom) <- c("x", "y", "d")
    return(hom)
  }

  project3Dhom
})

plot3Dpoints <- local({

  plot3Dpoints <- function(xyz, eye=c(2,-3,2), org=c(0,0,0),
                           ...,
                           type=c("p", "n", "h"),
                           xlim=c(0,1), ylim=c(0,1), zlim=c(0,1),
                           add=FALSE, box=TRUE, 
                           main, cex=par('cex'), 
                           box.back=list(col="pink"),
                           box.front=list(col="blue", lwd=2)
                           ) {
    if(missing(main)) main <- short.deparse(substitute(xyz))
    type <- match.arg(type)
    #'
    if(is.null(box.back) || (is.logical(box.back) && box.back))
      box.back <- list(col="pink")
    if(is.null(box.front) || (is.logical(box.front) && box.front))
      box.front <- list(col="blue", lwd=2)
    stopifnot(is.list(box.back) || is.logical(box.back))
    stopifnot(is.list(box.front) || is.logical(box.front))
    #'
    stopifnot(is.matrix(xyz) && ncol(xyz) == 3)
    if(nrow(xyz) > 0) {
      if(missing(xlim)) xlim <- range(pretty(xyz[,1]))
      if(missing(ylim)) ylim <- range(pretty(xyz[,2]))
      if(missing(zlim)) zlim <- range(pretty(xyz[,3]))
      if(missing(org)) org <- c(mean(xlim), mean(ylim), mean(zlim))
    }
    if(!add) {
      #' initialise plot
      bb <- plot3Dbox(xlim, ylim, zlim, eye=eye, org=org, do.plot=FALSE)
      plot(bb$xlim, bb$ylim, axes=FALSE, asp=1, type="n",
           xlab="", ylab="", main=main)
    }
    if(is.list(box.back)) {
      #' plot rear of box
      do.call(plot3DboxPart,
              resolve.defaults(list(xlim=xlim,
                                    ylim=ylim,
                                    zlim=zlim,
                                    eye=eye, org=org,
                                    part="back"),
                               box.back,
                               list(...)))
    }
    if(type != "n") {
      #' plot points
      uv <- project3Dhom(xyz, eye=eye, org=org)
      uv <- as.data.frame(uv)
      dord <- order(uv$d, decreasing=TRUE)
      uv <- uv[dord, , drop=FALSE]
      if(type == "h") {
        xy0 <- cbind(xyz[,1:2,drop=FALSE], zlim[1])
        uv0 <- as.data.frame(project3Dhom(xy0, eye=eye, org=org))
        uv0 <- uv0[dord, , drop=FALSE]
        do.call.matched(segments,
                        list(x0=with(uv0, x/d),
                             y0=with(uv0, y/d),
                             x1=with(uv,  x/d),
                             y1=with(uv,  y/d),
                             ...))
      }
      with(uv, points(x/d, y/d, cex=cex * min(d)/d, ...))
    }
    if(is.list(box.front)) 
      do.call(plot3DboxPart,
              resolve.defaults(list(xlim=xlim,
                                    ylim=ylim,
                                    zlim=zlim,
                                    eye=eye, org=org,
                                    part="front"),
                               box.front,
                               list(...)))
    return(invisible(NULL))
  }

  vertexind <- data.frame(i=rep(1:2,4),
                          j=rep(rep(1:2,each=2),2),
                          k=rep(1:2, each=4))

  edgepairs <- data.frame(from=c(1, 1, 2, 3, 1, 2, 5, 3, 5, 4, 6, 7),
                          to = c(2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8))
  
  vertexfrom <- vertexind[edgepairs$from,]
  vertexto   <- vertexind[edgepairs$to,]

  hamming <- function(a, b) sum(abs(a-b))

  ## determine projected positions of box vertices
  ## and optionally plot the box
  plot3Dbox <- function(xlim=c(0,1), ylim=xlim, zlim=ylim,
                        eye=c(0,-3,1), org=c(0,0,0),
                        do.plot=TRUE) {
    fromxyz <- with(vertexfrom, cbind(xlim[i], ylim[j], zlim[k]))
    toxyz   <- with(vertexto,   cbind(xlim[i], ylim[j], zlim[k]))
    fromuv <-  project3Dhom(fromxyz, eye=eye, org=org)
    touv <-  project3Dhom(toxyz, eye=eye, org=org)
    xfrom <- fromuv[,1]/fromuv[,3]
    xto   <- touv[,1]/touv[,3]
    yfrom <- fromuv[,2]/fromuv[,3]
    yto   <- touv[,2]/touv[,3]
    if(do.plot) 
      segments(xfrom, yfrom, xto, yto)
    return(invisible(list(xlim=range(xfrom, xto), ylim=range(yfrom, yto))))
  }

  ## plot either back or front of box
  plot3DboxPart <- function(xlim=c(0,1), ylim=xlim, zlim=ylim,
                            eye=c(0,-3,1), org=c(0,0,0),
                            part=c("front", "back"), ...) {
    part <- match.arg(part)
    boxvert <- with(vertexind, cbind(xlim[i], ylim[j], zlim[k]))
    pvert <- project3Dhom(boxvert, eye=eye, org=org)
    xyvert <- pvert[,c("x","y")]/pvert[,"d"]
    ## find vertex which is furthest away
    nback <- which.max(pvert[,"d"])
    nearback <- with(edgepairs, (from==nback) | (to==nback))
    ind <- if(part == "back") nearback else !nearback
    ## draw lines
    with(edgepairs[ind,],
         do.call.matched(segments,
                         list(x0=xyvert[from, 1],
                              y0=xyvert[from, 2],
                              x1=xyvert[to,   1],
                              y1=xyvert[to,   2],
                              ...)))
  }

  plot3Dpoints
})

