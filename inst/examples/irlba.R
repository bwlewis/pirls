library(irlba)
library(ddR)

setMethod("%*%", signature(x="ParallelObj", y="numeric"), function(x ,y)
{
  stopifnot(ncol(x) == length(y))
  collect(
    dmapply(function(a, b) a %*% b,
            parts(x),
            replicate(totalParts(x), y, FALSE),
            output.type="darray", combine="rbind", nparts=nparts(x)))
})
setMethod("%*%", signature(x="numeric", y="ParallelObj"), function(x ,y)
{
  stopifnot(length(x) == nrow(y))
  colSums(
    dmapply(function(x, y) x %*% y,
            split(x, rep(1:nparts(y)[1], psize(y)[, 1])),
            parts(y),
            output.type="darray", combine="rbind", nparts=nparts(y)))
})

x = dmapply(function(x) matrix(runif(4), 25, 25), 1:4, output.type="darray", combine="rbind", nparts=c(4, 1))
print(irlba(x, nv=3)$d)
print(svd(collect(x))$d[1:3])
