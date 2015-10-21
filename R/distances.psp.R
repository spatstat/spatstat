#
#  distances.psp.R
#
#  Hausdorff distance and Euclidean separation for psp objects
#
#  $Revision: 1.11 $ $Date: 2015/10/21 09:06:57 $
#
#

pairdist.psp <- function(X, ..., method="C", type="Hausdorff") {
  verifyclass(X, "psp")
  if(X$n == 0)
    return(matrix(, 0, 0))
  type <- pickoption("type", type,
                     c(Hausdorff="Hausdorff",
                       hausdorff="Hausdorff",
                       separation="separation"))

  D12 <- AsymmDistance.psp(X, X, metric=type, method=method)

  switch(type,
         Hausdorff={
           # maximum is Hausdorff metric
           D <- array(pmax.int(D12, t(D12)), dim=dim(D12))
         },
         separation={
           # Take minimum of endpoint-to-segment distances
           D <- array(pmin.int(D12, t(D12)), dim=dim(D12))
           # Identify any pairs of segments which cross
           cross <- test.selfcrossing.psp(X)
           # Assign separation = 0 to such pairs
           D[cross] <- 0
         })
  return(D)
}

crossdist.psp <- function(X, Y, ..., method="C", type="Hausdorff") {
  verifyclass(X, "psp")
  Y <- as.psp(Y)
  if(X$n * Y$n == 0)
    return(matrix(, X$n, Y$n))

  type <- pickoption("type", type,
                     c(Hausdorff="Hausdorff",
                       hausdorff="Hausdorff",
                       separation="separation"))
  
  DXY <- AsymmDistance.psp(X, Y, metric=type, method=method)
  DYX <- AsymmDistance.psp(Y, X, metric=type, method=method)
  
  switch(type,
         Hausdorff={
           # maximum is Hausdorff metric
           D <- array(pmax.int(DXY, t(DYX)), dim=dim(DXY))
         },
         separation={
           # Take minimum of endpoint-to-segment distances
           D <- array(pmin.int(DXY, t(DYX)), dim=dim(DXY))
           # Identify pairs of segments which cross
           cross <- test.crossing.psp(X, Y)
           # Assign separation = 0 to such pairs
           D[cross] <- 0
         })
  return(D)
}

nndist.psp <- function(X, ..., k=1, method="C") {
  verifyclass(X, "psp")
  if(!(is.vector(k) && all(k %% 1 == 0) && all(k >= 1)))
    stop("k should be a positive integer or integers")
  n <- nobjects(X)
  kmax <- max(k)
  lenk <- length(k)
  result <- if(lenk == 1) numeric(n) else matrix(, nrow=n, ncol=lenk)
  if(n == 0)
    return(result)
  if(kmax >= n) {
    # not enough objects 
    # fill with Infinite values
    result[] <- Inf
    if(any(ok <- (kmax < n))) {
      # compute the lower-order nnd's
      result[, ok] <- nndist.psp(X, ..., k=k[ok], method=method)
    }
    return(result)
  }
  # normal case:
  D <- pairdist.psp(X, ..., method=method)
  diag(D) <- Inf
  if(kmax == 1) 
    NND <- apply(D, 1, min)
  else 
    NND <- t(apply(D, 1, orderstats, k=k))[, , drop=TRUE]
  return(NND)
}

# .....  AsymmDistance.psp .....
#
# If metric="Hausdorff":
#     this function computes, for each pair of segments A = X[i] and B = Y[j],
#     the value max_{a in A} d(a,B) = max_{a in A} min_{b in B} ||a-b||
#     which appears in the definition of the Hausdorff metric.
#     Since the distance function d(a,B) of a segment B is a convex function,
#     the maximum is achieved at an endpoint of A. So the algorithm
#     actually computes h(A,B) = max (d(e_1,B), d(e_2,B)) where e_1, e_2
#     are the endpoints of A. And H(A,B) = max(h(A,B),h(B,A)).
#
# If metric="separation":
#     the function computes, for each pair of segments A = X[i] and B = Y[j],
#     the MINIMUM distance from an endpoint of A to any point of B.
#        t(A,B) = min (d(e_1,B), d(e_2,B))
#     where e_1, e_2 are the endpoints of A.
#     Define the separation distance
#        s(A,B) = min_{a in A} min_{b in B} ||a-b||.
#     The minimum (a*, b*) occurs either when a* is an endpoint of A,
#     or when b* is an endpoint of B, or when a* = b* (so A and B intersect).
#     (If A and B are parallel, the minimum is still achieved at an endpoint)
#     Thus s(A,B) = min(t(A,B), t(B,A)) unless A and B intersect.


AsymmDistance.psp <- function(X, Y, metric="Hausdorff",
                              method=c("C", "Fortran", "interpreted")) {
  method <- match.arg(method)
  # Extract endpoints of X
  EX <- endpoints.psp(X, "both")
  idX <- attr(EX, "id")
  # compute shortest dist from each endpoint of X to each segment of Y
  DPL <- distppll(cbind(EX$x,EX$y), Y$ends, mintype=0, method=method)
  # for each segment in X, maximise or minimise over the two endpoints
  Dist <- as.vector(DPL)
  Point <- as.vector(idX[row(DPL)])
  Segment <- as.vector(col(DPL))
  switch(metric,
         Hausdorff={
           DXY <- tapply(Dist, list(factor(Point), factor(Segment)), max)
         },
         separation={
           DXY <- tapply(Dist, list(factor(Point), factor(Segment)), min)
           })
  return(DXY)
}
