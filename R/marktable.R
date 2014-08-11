#
#	marktable.R
#
#	Tabulate mark frequencies in r-neighbourhood of each point 
#	for multitype point patterns
#
#	$Revision: 1.6 $	$Date: 2013/02/22 05:34:02 $
#
#       Requested by Ian Robertson <igr@stanford.edu>


"marktable" <- 
function(X, R, exclude=TRUE) 
{
	verifyclass(X, "ppp")
	if(!is.marked(X, dfok=FALSE))
		stop("point pattern has no marks")
        stopifnot(is.numeric(R) && length(R) == 1 && R > 0)
        stopifnot(is.logical(exclude) && length(exclude) == 1)

        m <- marks(X)
        if(!is.factor(m))
          stop("marks must be a factor")
        
        # identify close pairs
        p <- closepairs(X,R,what="indices")
        pi <- p$i
        pj <- p$j
        if(!exclude) {
          # add identical pairs
          n <- X$n
          pi <- c(pi, 1:n)
          pj <- c(pj, 1:n)
        }

        # tabulate
        i <- factor(pi, levels=seq_len(X$n))
        mj <- m[pj]
        mat <- table(point=i, mark=mj)

        return(mat)
}

