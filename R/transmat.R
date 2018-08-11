## transmat.R
##
## transform matrices between different spatial indexing conventions
##
##  $Revision: 1.1 $  $Date: 2015/03/04 07:13:10 $

transmat <- local({

  euro <- matrix(c(0,-1,1,0), 2, 2)
  spat <- matrix(c(0,1,1,0), 2, 2)
  cart <- diag(c(1,1))
  dimnames(euro) <- dimnames(spat) <- dimnames(cart) <- 
    list(c("x","y"), c("i","j"))

  known <- list(spatstat=spat,
                cartesian=cart,
                Cartesian=cart,
                european=euro,
                European=euro)

  cmap <- list(x=c(1,0),
               y=c(0,1),
               i=c(1,0),
               j=c(0,1))
  
  maptocoef <- function(s) { 
    e <- parse(text=s)[[1]]
    eval(eval(substitute(substitute(f, cmap), list(f=e)))) 
  }

  
  as.convention <- function(x) {
    if(is.character(x) && length(x) == 1) {
      k <- pmatch(x, names(known))
      if(is.na(k)) 
        stop(paste("Unrecognised convention", sQuote(x)), call.=FALSE)
      return(known[[k]])
    }
    if(is.list(x) && is.character(unlist(x))) {
      xx <- lapply(x, maptocoef)
      if(all(c("x", "y") %in% names(xx))) z <- rbind(xx$x, xx$y) else
      if(all(c("i", "j") %in% names(xx))) z <- cbind(xx$x, xx$y) else 
      stop("entries should be named i,j or x,y", call.=FALSE)
      dimnames(z) <- list(c("x","y"), c("i","j"))
      if(!(all(z == 0 | z == 1 | z == -1) && 
           all(rowSums(abs(z)) == 1) && 
           all(colSums(abs(z)) == 1)))
        stop("Illegal convention", call.=FALSE)
      return(z)
    }
    stop("Unrecognised format for spatial convention", call.=FALSE)
  }  

  transmat <- function(m, from, to) {
    m <- as.matrix(m)
    from <- as.convention(from)
    to <- as.convention(to)
    conv <- solve(from) %*% to
    flip <- apply(conv == -1, 2, any)
    if(flip[["i"]]) m <- m[nrow(m):1, , drop=FALSE]
    if(flip[["j"]]) m <- m[         , ncol(m):1, drop=FALSE]
    if(all(diag(conv) == 0))
       m <- t(m)
    return(m)
  }
  
  transmat
})
