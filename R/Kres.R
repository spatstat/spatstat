#
#	Kres.R
#
#	Residual K
#
#	$Revision: 1.2 $	$Date: 2013/02/26 04:32:27 $
#
#############################################################################
#

Kres <- function(object, ...) {
  if(!is.fv(object)) {
    # usual case where 'object' is a ppm, ppp or quad
    K <- Kcom(object, ...)
  } else {
    # case where 'object' is the output of 'Kcom'
    a <- attr(object, "maker")
    if(is.null(a) || a != "Kcom")
      stop("fv object was not created by Kcom")
    K <- object
    if(length(list(...)) > 0)
      warning("Extra arguments ignored")
  }
  # initialise fv object
  df <- data.frame(r=K$r, theo=rep(0, length(K$r)))
  desc <- c("distance argument r", "value 0 corresponding to perfect fit")
  ans <- fv(df, "r", substitute(bold(R)~hat(K)(r), NULL),
            "theo", . ~ r ,
            attr(K, "alim"), c("r","bold(R)~%s[theo](r)"), desc, fname="K")
  # add residual functions
  nam <- names(K)
  if("border" %in% nam)
    ans <- bind.fv(ans,
                    data.frame(bres=with(K, border-bcom)),
                    "bold(R)~hat(%s)[bord](r)",
                    "residual function %s based on border correction",
                    "bres")
  if(all(c("trans","tcom") %in% nam))
    ans <- bind.fv(ans,
                    data.frame(tres=with(K, trans-tcom)),
                    "bold(R)~hat(%s)[trans](r)",
                    "residual function %s based on translation correction",
                    "tres")
  if(all(c("iso","icom") %in% nam))
    ans <- bind.fv(ans,
                    data.frame(ires=with(K, iso-icom)),
                    "bold(R)~hat(%s)[iso](r)",
                    "residual function %s based on isotropic correction",
                    "ires")
  if("ivar" %in% nam) {
    savedotnames <- fvnames(ans, ".")
    ans <- bind.fv(ans,
                   as.data.frame(K)[, c("ivar", "isd", "ihi", "ilo")],
                    c("bold(C)^2~hat(%s)[iso](r)",
                      "sqrt(bold(C)^2~hat(%s)[iso](r))",
                      "bold(R)~hat(%s)[Hi](r)",
                      "bold(R)~hat(%s)[Lo](r)"),
                    c("pseudovariance of isotropic-corrected residual %s",
                      "pseudo-SD of isotropic-corrected residual %s",
                      "upper critical band for isotropic-corrected residual %s",
                      "lower critical band for isotropic-corrected residual %s"),
                    "ires")
    ans <- bind.fv(ans,
                   data.frame(istdres=with(ans, ires/isd)),
                   "bold(T)~hat(%s)[iso](r)",
                   "standardised isotropic-corrected residual %s",
                   "ires")
    fvnames(ans, ".") <- c(savedotnames, c("ihi", "ilo"))
  }
  unitname(ans) <- unitname(K)
  return(ans)
}
