#'
#' $Revision: 1.1 $ $Date: 2017/06/05 10:31:58 $
#'

clicklpp <- local({

  clicklpp <- function(L, n=NULL, types=NULL, ...,
                       add=FALSE, main=NULL, hook=NULL) {
    if(!inherits(L, "linnet"))
      stop("L should be a linear network", call.=FALSE)
    instructions <-
      if(!is.null(n)) paste("click", n, "times in window") else
      paste("add points: click left mouse button in window\n",
            "exit: press ESC or another mouse button")
    if(is.null(main))
      main <- instructions
    W <- Window(L)
  
    ####  single type #########################
    if(is.null(types)) {
      plot(L, add=add, main=main)
      if(!is.null(hook))
        plot(hook, add=TRUE)
      xy <- if(!is.null(n)) spatstatLocator(n=n, ...) else spatstatLocator(...)
      ok <- inside.owin(xy, w=W)
      if((nbad <- sum(!ok)) > 0) 
        warning(paste("Ignored",
	              nbad,
	              ngettext(nbad, "point", "points"),
		      "outside window"),
	        call.=FALSE)
      X <- as.lpp(xy$x[ok], xy$y[ok], L=L)
      return(X)
    }
  
    ##### multitype #######################
    
    ftypes <- factor(types, levels=types)
    #' input points of type 1 
    X <- getem(ftypes[1L], instructions, n=n, L=L, add=add, ..., pch=1)
    X <- X %mark% ftypes[1L]
    #' input points of types 2, 3, ... in turn
    for(i in 2:length(types)) {
      Xi <- getem(ftypes[i], instructions, n=n, L=L, add=add,
                  ..., hook=X, pch=i)
      Xi <- Xi %mark% ftypes[i]
      X <- superimpose(X, Xi, L=L)
    }
    if(!add) 
      plot(X, main="Final pattern")
    return(X)
  }

  getem <- function(i, instr, ...) {
    main <- paste("Points of type", sQuote(i), "\n", instr)
    do.call(clicklpp, resolve.defaults(list(...), list(main=main)))
  }

  clicklpp
})


