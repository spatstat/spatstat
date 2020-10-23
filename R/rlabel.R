#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.13 $   $Date: 2020/10/23 15:18:00 $
#
#
rlabel <- local({

  resample <- function(x, replace=FALSE) {
    x[sample(length(x), replace=replace)]
  }

  rlabel <- function(X, labels=marks(X), permute=TRUE,
                     group=NULL, ..., nsim=1, drop=TRUE) {
    stopifnot(is.ppp(X) || is.lpp(X) || is.pp3(X) || is.ppx(X) || is.psp(X))
    if(is.null(labels))
      stop("labels not given and marks not present")
    singlecolumn <- (length(dim(labels)) < 2)
    nthings <- nobjects(X)
    things <- if(is.psp(X)) "segments" else "points"
    nlabels <- if(singlecolumn) length(labels) else nrow(labels)
    if((nlabels != nthings) && (permute || !is.null(group))) 
      stop(paste(if(singlecolumn) "Length" else "Number of rows",
                 "of labels does not match the number of",
                 things),
           call.=FALSE)
    ##
    if(is.null(group)) {
      Y <- replicate(nsim, {
        X %mark% marksubset(labels, sample(nlabels, nthings, replace=!permute))
      },
      simplify=FALSE)
    } else {
      group <- marks(cut(X, group, ...))
      seqn <- seq_len(nlabels)
      pieces <- split(seqn, group)
      Y <- replicate(nsim, {
        X %mark% marksubset(labels,
                            unsplit(lapply(pieces, resample, replace=!permute),
                                    group))
      },
      simplify=FALSE)
    }
    ## 
    return(simulationresult(Y, nsim, drop))
  }

  rlabel
})
