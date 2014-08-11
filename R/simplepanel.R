#
# simplepanel.R
#
#  A simple, robust point & click interface
#     used in rmh visual debugger.
#
#  $Revision: 1.11 $  $Date: 2013/12/08 00:08:23 $
#

simplepanel <- function(title, B, boxes, clicks, redraws=NULL, exit=NULL, env) {
  stopifnot(is.rectangle(B))
  stopifnot(is.list(boxes))
  if(!all(unlist(lapply(boxes, is.rectangle))))
    stop("some of the boxes are not rectangles")
  if(!all(unlist(lapply(boxes, is.subset.owin, B=B))))
    stop("Some boxes do not lie inside the bounding box B")
  stopifnot(is.list(clicks) && length(clicks) == length(boxes))
  if(!all(unlist(lapply(clicks, is.function))))
    stop("clicks must be a list of functions")
  if(is.null(redraws)) {
    redraws <- rep.int(list(dflt.redraw), length(boxes))
  } else {
    stopifnot(is.list(redraws) && length(redraws) == length(boxes))
    if(any(isnul <- unlist(lapply(redraws, is.null))))
      redraws[isnul] <- rep.int(list(dflt.redraw), sum(isnul))
    if(!all(unlist(lapply(redraws, is.function))))
      stop("redraws must be a list of functions")
  }
  if(is.null(exit)) {
    exit <- function(...) { NULL}
  } else stopifnot(is.function(exit))
  stopifnot(is.environment(env))
  n <- length(boxes)
  got.boxnames <- (sum(nzchar(names(boxes))) == n)
  got.clicknames <- (sum(nzchar(names(clicks))) == n)
  nama <- if(got.boxnames && !got.clicknames) names(boxes) else
          if(got.clicknames && !got.boxnames) names(clicks) else
          paste("Button", seq_len(n))
  out <- list(title=title, B=B,
              nama=nama, boxes=boxes, clicks=clicks, redraws=redraws,
              exit=exit, env=env)
  class(out) <- c("simplepanel", class(out))
  return(out)
}

grow.simplepanel <- function(P, side=c("right","left","top","bottom"),
                             len=NULL,
                             new.clicks, new.redraws=NULL, ..., aspect) {
  verifyclass(P, "simplepanel")
  side <- match.arg(side)
  stopifnot(is.list(new.clicks))
  if(!all(unlist(lapply(new.clicks, is.function))))
    stop("new.clicks must be a list of functions")
  if(is.null(new.redraws)) {
    new.redraws <- rep.int(list(dflt.redraw), length(new.clicks))
  } else {
    stopifnot(is.list(new.redraws) && length(new.redraws) == length(new.clicks))
    if(!all(unlist(lapply(new.redraws, is.function))))
      stop("new.redraws must be a list of functions")
  }
  if(missing(aspect) || is.null(aspect)) {
    # determine aspect ratio from length of longest text string
    n <- length(new.clicks)
    nama <- names(new.clicks)
    if(sum(nzchar(nama)) != n)
      nama <- names(new.redraws)
    if(sum(nzchar(nama)) != n)
      nama <- paste("Box", seq_len(n))
    aspect <- 3/max(4, nchar(nama))
  }
  B <- P$B
  n <- length(new.clicks)
  switch(side,
         right={
           new.width <- if(!is.null(len)) len else sidelengths(B)[1]/2
           extraspace <- owin(B$xrange[2] + c(0, new.width), B$yrange)
           new.boxes <- layout.boxes(extraspace, n, ..., aspect=aspect)
         },
         left={
           new.width <- if(!is.null(len)) len else sidelengths(B)[1]/2
           extraspace <- owin(B$xrange[1] - c(new.width, 0), B$yrange)
           new.boxes <- layout.boxes(extraspace, n, ..., aspect=aspect)
         },
         top={
           new.height <- if(!is.null(len)) len else sidelengths(B)[2]/2
           extraspace <- owin(B$xrange, B$yrange[2] + c(0, new.height))
           new.boxes <- layout.boxes(extraspace, n, ..., aspect=aspect,
                                     horizontal=TRUE)
         },
         bottom={
           new.height <- if(!is.null(len)) len else sidelengths(B)[2]/2
           extraspace <- owin(B$xrange, B$yrange[1] - c(new.height, 0))
           new.boxes <- layout.boxes(extraspace, n, ..., aspect=aspect,
                                     horizontal=TRUE)
         })
  with(P, simplepanel(title,
                      bounding.box(B, extraspace),
                      append(boxes, new.boxes),
                      append(clicks, new.clicks),
                      append(redraws, new.redraws),
                      exit, env))
}

                             
redraw.simplepanel <- function(P, verbose=FALSE) {
  verifyclass(P, "simplepanel")
  if(verbose)
    cat("Redrawing entire panel\n")
  with(P, {
    ntitle <- sum(nzchar(title))
    plot(B, type="n", main=title)
    for(j in seq_along(nama)) 
      (redraws[[j]])(boxes[[j]], nama[j], env)
  })
  invisible(NULL)
}

clear.simplepanel <- function(P) {
  verifyclass(P, "simplepanel")
  plot(P$B, main="")
  invisible(NULL)
}
                             
run.simplepanel <- function(P, popup=TRUE, verbose=FALSE) {
  verifyclass(P, "simplepanel")
  if(popup) dev.new()
  ntitle <- sum(nzchar(P$title))
  opa <- par(mar=c(0,0,ntitle+0.2,0),ask=FALSE)
  with(P, {
    # interaction loop
    more <- TRUE
    while(more) {
      redraw.simplepanel(P, verbose=verbose)
      xy <- locator(1)
      if(is.null(xy)) {
        if(verbose) cat("No (x,y) coordinates\n")
        break
      }
      found <- FALSE
      for(j in seq_along(boxes)) {
        if(inside.owin(xy$x, xy$y, boxes[[j]])) {
          found <- TRUE
          if(verbose) cat(paste("Caught click on", sQuote(nama[j]), "\n"))
          more <- (clicks[[j]])(env, xy)
          if(!is.logical(more) || length(more) != 1) {
            warning(paste("Click function for",
                          sQuote(nama[j]),
                          "did not return TRUE/FALSE"))
            more <- FALSE
          }
          if(verbose) cat(if(more) "Continuing\n" else "Terminating\n")
          break
        }
      }
      if(verbose && !found)
        cat(paste("Coordinates", paren(paste(xy, collapse=",")),
                  "not matched to any box\n"))
    }
  })
  if(verbose)
    cat("Calling exit function\n")

  rslt <- with(P, exit(env))
  
  # revert to original graphics parameters
  par(opa)
  # close popup window?
  if(popup) dev.off()
  
  # return value of 'exit' function
  return(rslt)
}

layout.boxes <- function(B, n, horizontal=FALSE, aspect=0.5, usefrac=0.9){
  # make n boxes in B
  stopifnot(is.rectangle(B))
  stopifnot(n > 0)
  width <- sidelengths(B)[1]
  height <- sidelengths(B)[2]
  if(!horizontal) {
    heightshare <- height/n
    useheight <- min(width * aspect, heightshare * usefrac)
    usewidth <-  min(useheight /aspect, width * usefrac)
    lostwidth <- width - usewidth
    lostheightshare <- heightshare - useheight
    template <- owin(c(0, usewidth), c(0, useheight))
    boxes <- list()
    boxes[[1]] <- shift(template,
                        c(B$xrange[1]+lostwidth/2,
                          B$yrange[1] + lostheightshare/2))
    if(n > 1) 
      for(j in 2:n) 
        boxes[[j]] <- shift(boxes[[j-1]], c(0, heightshare))
  } else {
    boxes <- layout.boxes(flipxy(B), n,
                            horizontal=FALSE, aspect=1/aspect, usefrac=usefrac)
    boxes <-  lapply(boxes, flipxy)
  }
  return(boxes)
}

# default redraw function for control buttons

dflt.redraw <- function(button, name, env) {
  plot(button, add=TRUE, border="pink")
  text(centroid.owin(button), labels=name)
}

print.simplepanel <- function(x, ...) {
  nama <- x$nama
  cat("simplepanel object\n")
  cat(paste("\tTitle:", sQuote(x$title), "\n"))
  cat("\tPanel names:")
  for(i in seq_along(nama)) {
    if(i %% 6 == 1) cat("\n\t")
    cat(paste0(sQuote(nama[i]), "  "))
  }
  cat("\n")
  return(invisible(NULL))
}
