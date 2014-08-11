#
#	exactPdt.R
#	R function exactPdt() for exact distance transform of pixel image
#
#	$Revision: 4.13 $	$Date: 2012/04/06 09:49:56 $
#

"exactPdt"<-
  function(w)
{
  verifyclass(w, "owin")
  if(w$type != "mask")
    stop(paste("Input must be a window of type", sQuote("mask")))
#	
  nr <- w$dim[1]
  nc <- w$dim[2]
# input image will be padded out with a margin of width 2 on all sides
  mr <- mc <- 2
  # full dimensions of padded image
  Nnr <- nr + 2 * mr
  Nnc <- nc + 2 * mc
  N <- Nnr * Nnc
  # output image (subset): rows & columns (R indexing)
  rmin <- mr + 1
  rmax <- Nnr - mr
  cmin <- mc + 1
  cmax <- Nnc - mc
  # do padding
  x <- matrix(FALSE, nrow=Nnr, ncol=Nnc)
  x[rmin:rmax, cmin:cmax] <- w$m
  #
  DUP <- spatstat.options("dupC")
  res <- .C("ps_exact_dt_R",
            as.double(w$xrange[1]),
            as.double(w$yrange[1]),
            as.double(w$xrange[2]),
            as.double(w$yrange[2]),
            nr = as.integer(nr),
            nc = as.integer(nc),
            mr = as.integer(mr),
            mc = as.integer(mc),
            inp = as.integer(t(x)),
            distances = as.double (double(N)),
            rows      = as.integer(integer(N)),
            cols      = as.integer(integer(N)),
            boundary  = as.double (double(N)),
            DUP=DUP,
            PACKAGE="spatstat"
            )
  dist <- matrix(res$distances,
                 ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
  rows <- matrix(res$rows,
                 ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
  cols <- matrix(res$cols,
                 ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
  bdist<- matrix(res$boundary,
                 ncol=Nnc, nrow=Nnr, byrow = TRUE)[rmin:rmax, cmin:cmax]
  # convert from C to R indexing
  rows <- rows + 1L
  cols <- cols + 1L
  return(list(d=dist,row=rows,col=cols,b=bdist, w=w))
}
