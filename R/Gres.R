#
#	Gres.R
#
#	Residual G 
#
#	$Revision: 1.3 $	$Date: 2013/04/25 06:37:43 $
#
#############################################################################
#

Gres <- function(object, ...) {
  if(!is.fv(object)) {
    # usual case where 'object' is a ppm, ppp or quad
    G <- Gcom(object, ...)
  } else {
    # case where 'object' is the output of 'Gcom'
    a <- attr(object, "maker")
    if(is.null(a) || a != "Gcom")
      stop("fv object was not created by Gcom")
    G <- object
    if(length(list(...)) > 0)
      warning("Extra arguments ignored")
  }
  # initialise fv object
  df <- data.frame(r=G$r, theo=numeric(length(G$r)))
  desc <- c("distance argument r", "value 0 corresponding to perfect fit")
  ans <- fv(df, "r", substitute(bold(R)~hat(G)(r), NULL),
            "theo", . ~ r,
            attr(G, "alim"), c("r","bold(R)~%s[theo](r)"), desc, fname="G")
  # add residual estimates
  nam <- names(G)
  if(all(c("border","bcom") %in% nam))
    ans <- bind.fv(ans,
                    data.frame(bres=with(G, border-bcom)),
                    "bold(R)~hat(%s)[bord](r)",
                    "border corrected residual of %s",
                    "bres")
  if(all(c("han","hcom") %in% nam))
    ans <- bind.fv(ans,
                    data.frame(hres=with(G, han-hcom)),
                    "bold(R)~hat(%s)[han](r)",
                    "Hanisch corrected residual of %s",
                    "hres")
  if("hvar" %in% nam) {
    savedotnames <- fvnames(ans, ".")
    hsd <- with(G, sqrt(hvar))
    ans <- bind.fv(ans,
                    data.frame(hvar=with(G, hvar),
                               hsd = hsd,
                               hi =  2*hsd,
                               lo = -2*hsd),
                    c("bold(C)^2~hat(%s)[han](r)",
                      "sqrt(bold(C)^2~hat(%s)[han](r))",
                      "bold(R)~hat(%s)[Hi](r)",
                      "bold(R)~hat(%s)[Lo](r)"),
                    c("pseudovariance of Hanisch corrected residual %s",
                      "pseudo-SD of Hanisch corrected residual %s",
                      "upper critical band for Hanisch corrected residual %s",
                      "lower critical band for Hanisch corrected residual %s"),
                    "hres")
    ans <- bind.fv(ans,
                   data.frame(hstdres=with(ans, hres/hsd)),
                   "bold(T)~hat(%s)[han](r)",
                   "standardised Hanisch-corrected residual %s",
                   "hres")
    fvnames(ans, ".") <- c(savedotnames, c("hi", "lo"))
  }
  unitname(ans) <- unitname(G)
  return(ans)
}
