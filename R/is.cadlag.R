is.cadlag <- function (s) 
{
if(!is.stepfun(s)) stop("s is not a step function.\n")
r <- knots(s)
h <- s(r)
n <- length(r)
r1 <- c(r[-1L],r[n]+1)
rm <- (r+r1)/2
hm <- s(rm)
isTRUE(all.equal(h,hm))
}
