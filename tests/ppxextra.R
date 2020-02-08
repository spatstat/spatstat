require(spatstat)

## Tests for shift:
# Check ppp and ppx shift are the same
Xppp <- cells
Xppx <- ppx(coords(cells), domain = boxx(0:1,0:1))
Xppx2 <- shift(Xppx, vec = c(1,1))
stopifnot(all.equal(coords(Xppp2), coords(Xppx2), check.attributes = FALSE))
stopifnot(all.equal(domain(Xppp2), as.owin(domain(Xppx2)), check.attributes = FALSE))
# Check a single numeric for vec in shift.ppx
stopifnot(identical(Xppx2, shift(Xppx, vec = 1)))


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
stopifnot(all.equal(coords(Xscaled),coords(Yscaled), check.attributes = FALSE))
# Check ppx with domain:
Y$domain <- boxx(xrange, yrange)
Yscaled <- scale(Y, center = cent, scale = scal)
stopifnot(identical(as.boxx(Window(Xscaled)), domain(Yscaled)))


## Tests for intersect.boxx:
# Should be unit 2D box:
A <- intersect.boxx(boxx(c(-1,1),c(0,2)), boxx(c(0,3),c(0,1)))
stopifnot(identical(A, boxx(c(0,1),c(0,1))))
# Should be empty (NULL)
B <- intersect.boxx(boxx(c(-1,1),c(0,2)), boxx(c(0,3),c(0,1)), boxx(c(1,2), c(-1,1)))
stopifnot(is.null(B))
# Should be unit 3D box:
C <- intersect.boxx(boxx(c(-1,1),c(0,2),c(-1,1)), boxx(c(0,3),c(0,1),c(0,4)))
stopifnot(identical(C, boxx(c(0,1),c(0,1),c(0,1))))
# Should be empty (NULL)
D <- intersect.boxx(boxx(c(-1,1),c(0,2),c(-1,1)), boxx(c(0,3),c(0,1),c(0,4)), NULL)
stopifnot(is.null(D))


## Tests for [.boxx with clip:
# Check ppp and ppx subset with clip are the same
X <- cells
WX <- shift(domain(X), vec = c(.5,.5))
X2 <- X[WX, clip=TRUE]
Y <- ppx(coords(X), domain = boxx(c(0,1),c(0,1)))
WY <- shift(domain(Y), vec = c(.5,.5))
Y2 <- Y[WY, clip=TRUE]
stopifnot(all.equal(coords(X2), coords(Y2), check.attributes = FALSE))
stopifnot(all.equal(domain(X2), as.owin(domain(Y2))))
