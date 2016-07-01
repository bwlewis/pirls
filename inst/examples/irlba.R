library(irlba)
library(ddR)

### Utility functions

# Compute a %*% x where a is a darray and x is a vector
mvec = function(a, x)
{
  collect(
    dmapply(function(x, y) x %*% y,
            parts(a),
            replicate(totalParts(a), x, FALSE),
            output.type="darray", combine="rbind", nparts=nparts(a)))
}
# Compute x %*% a where a is a darray and x is a vector
vecm = function(x, a)
{
  colSums(
    dmapply(function(x, y) x %*% y,
            split(x, rep(1:nparts(a)[1], psize(a)[, 1])),
            parts(a),
            output.type="darray", combine="rbind", nparts=nparts(a)))
}
setMethod("%*%", signature(x="ParallelObj", y="numeric"), function(x ,y)
{
  stopifnot(ncol(x) == length(y))
  mvec(x, y)
})
setMethod("%*%", signature(x="numeric", y="ParallelObj"), function(x ,y)
{
  stopifnot(length(x) == nrow(y))
  vecm(x, y)
})

x = dmapply(function(x) matrix(runif(4), 25, 25), 1:4, output.type="darray", combine="rbind", nparts=c(4, 1))
print(irlba(x, nv=3)$d)
print(svd(collect(x))$d[1:3])
