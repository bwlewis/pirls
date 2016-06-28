library(ddR)

# Compute a %*% x where x is a vector
# XXX rbind is expensive, consider allocating the result and explicitly filling in
mvec = function(a, x)
{
  Reduce(rbind, Map(function(j) collect(dmapply(function(A, x) A %*% x, parts(a), replicate(totalParts(a), x, FALSE), output.type="darray",combine="rbind", nparts=nparts(a)), j), seq(1,totalParts(a))))
}

# Cross product t(b) %*% a for vector b and darray a, returns a vector
cross = function(a, b)
{
  Reduce(`+`, Map(function(j) collect(dmapply(function(x, y) crossprod(x, y), parts(a), split(b, rep(1:nparts(a)[1], psize(a)[,1])), output.type="darray",combine="rbind", nparts=nparts(a)), j), seq(1,totalParts(a))))
}

# Compute t(a) %*% diag(w) %*% a for darray a and vector w, returning a matrix
wcross = function (a, w)
{
  Reduce(`+`, Map(function(j) collect(dmapply(function(x, y) crossprod(x, y*x), parts(a), split(w, rep(1:nparts(a)[1], psize(a)[,1])), output.type="darray",combine="rbind", nparts=nparts(a)), j), seq(1,totalParts(a))))
}


# decided to use a slightly modified function instead of just overloading
# everything because the weighted crossproduct is cheaper this way.

dirls =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  x = rep(0,ncol(A))
  for(j in 1:maxit)
  {
    eta    = mvec(A, x)
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    xold   = x
    x      = solve(wcross(A, W), cross(A, W*z), tol=2*.Machine$double.eps)
    if(sqrt(crossprod(x - xold)) < tol) break
  }
  list(coefficients=x, iterations=j)
}

x = dmapply(function(x) matrix(runif(4),2,2),1:4,output.type="darray",combine="rbind",nparts=c(4,1))
y = 1:8

print(dirls(x, y, gaussian)$coef)
print(coef(lm.fit(collect(x), y)))
