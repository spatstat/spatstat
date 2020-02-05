## Tests for shift:
# Check ppp and ppx shift are the same
Xppp <- cells
Xppx <- ppx(coords(cells), domain = boxx(0:1,0:1))
Xppx2 <- shift(Xppx, vec = c(1,1))
all.equal(coords(Xppp2), coords(Xppx2), check.attributes = FALSE)
all.equal(domain(Xppp2), as.owin(domain(Xppx2)), check.attributes = FALSE) 

# Check a single numeric for vec in shift.ppx
identical(Xppx2, shift(Xppx, vec = 1))


## Tests for scale:
dat <- data.frame(x=1:3, y=1:3, m=letters[1:3])
xrange <- yrange <- c(0,4)
cent <- c(2,2)
scal <- c(5,5)
X <- as.ppp(dat, W = owin(xrange, yrange))
Xscaled <- affine(shift(X, vec = -cent), mat = diag(1/scal))

# Check ppx without domain:
Y <- ppx(dat, coord.type = c("spatial", "spatial", "mark"))
Yscaled <- scale(Y, center = cent, scale = scal)
all.equal(coords(Xscaled),coords(Yscaled), check.attributes = FALSE)

# Check ppx with domain:
Y$domain <- boxx(xrange, yrange)
Yscaled <- scale(Y, center = cent, scale = scal)
identical(as.boxx(Window(Xscaled)), domain(Yscaled))
