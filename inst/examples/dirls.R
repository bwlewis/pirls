library(pirls)
library(ddR)

# Compute a %*% x where a is a darray and x is a vector
# XXX rbind is expensive, consider allocating the result and explicitly filling in XXX
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


#' Distributed iteratively reweighted least squares (IRLS) algorithm using ddR
#' @param x design matrix of dimension n * p
#' @param y observation vector of length n
#' @param family a description of the error distribution and link function to
#'          be used in the model specified as a family function name
#'          the third option is supported.  (See ‘family’ for details of
#'          family functions.)
#' @param maxit maximum number of IRLS iterations
#' @param tol convergence tolerance, stop iterations when the
#' Euclidean norem of the update of the
#' model coefficients is less than tol.
#' @return a list with entries \code{coefficients} (model coefficients)
#' and \code{iterations} (number of IRLS iterations performed).
#' @seealso \code{\link{glm.fit}}
#' @examples
#' x = dmapply(function(x) matrix(runif(4), 2, 2), 1:4, output.type="darray", combine="rbind", nparts=c(4, 1))
#' y = 1:8
#' print(coef(dirls(x, y, gaussian)))
#' print(coef(lm.fit(collect(x), y)))
dirls =
function(x, y, family=binomial, maxit=25, tol=1e-08)
{
  b = rep(0, ncol(x))
  for(j in 1:maxit)
  {
    eta    = mvec(x, b)
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (y - g) / gprime
    W      = drop(gprime^2 / family()$variance(g))
    bold   = b
    b      = solve(wcross(x, W), cross(x, W*z), tol=2*.Machine$double.eps)
    if(sqrt(crossprod(b - bold)) < tol) break
  }
  list(coefficients=b, iterations=j)
}




# A trivially basic parallel darray times vector product, used in quadratic_loss
setMethod("%*%", signature(x="ParallelObj", y="numeric"), function(x ,y)
{
  stopifnot(ncol(x) == length(y))
  mvec(x, y)
})
setMethod("%*%", signature(x="ParallelObj", y="matrix"), function(x ,y)
{
  stopifnot(ncol(x) == nrow(y))
  mvec(x, drop(y))
})

# Compute a[, -i] %*% x[-i] where x is a vector and a is a darray
# XXX rbind is expensive, consider allocating the result and explicitly filling in XXX
mvec_holdout = function(a, x, i)
{
  Reduce(rbind, Map(function(j) collect(dmapply(function(A, x) A[, -i, drop=FALSE] %*% drop(x)[-i], parts(a), replicate(totalParts(a), x, FALSE), output.type="darray",combine="rbind", nparts=nparts(a)), j), seq(1,totalParts(a))))
}


# Prototype darray coordinate descent that can be used with pirls
dcoordinate_descent = function(X, W, z, lambda, alpha, beta, maxit)
{

  dupdate_coordinates = function(X, W, z, lambda, alpha, beta)
  {
    beta_old = beta
    for (i in 1:length(beta)) {
#      beta[i] = soft_thresh_r(sum(W*X[1:nrow(X), i]*(z - X[,-i] %*% beta_old[-i])),
      beta[i] = soft_thresh_r(sum(W*X[1:nrow(X), i]*(z - mvec_holdout(X, beta_old, i))),
                                 sum(W)*lambda*alpha)
    }
#    beta / (colSums(W*X^2) + lambda*(1-alpha))
    beta / (colSums(dmapply(function(a, b) b  * a ^ 2, parts(x), split(W, rep(1:nparts(X)[1], psize(X)[,1])), output.type="darray", combine="rbind", nparts=nparts(X))) + lambda * (1 - alpha))
  }

  quad_loss = pirls:::quadratic_loss(X, W, z, lambda, alpha, beta)
  for(i in 1:maxit) {
    beta_old = beta
    quad_loss_old = quad_loss
    beta = dupdate_coordinates(X, W, z, lambda, alpha, beta)
    quad_loss = pirls:::quadratic_loss(X, W, z, lambda, alpha, beta)
    if(quad_loss >= quad_loss_old) {
      beta = beta_old
      break
    }
  }
  if (i == maxit && quad_loss <= quad_loss_old) {
    warning("Coordinate descent did not converge.")
  }
  beta
}



# EXAMPLES
x = dmapply(function(x) matrix(runif(4), 2, 2), 1:4, output.type="darray", combine="rbind", nparts=c(4, 1))
y = 1:8

# IRLS
print(coef(dirls(x, y, gaussian)))
print(coef(lm.fit(collect(x), y)))

# PIRLS
print(pirls(collect(x), y, lambda=0, alpha=1, family=gaussian, beta=matrix(0, nrow=2, ncol=1)))
print(pirls(x, y, lambda=0, alpha=1, family=gaussian, beta=matrix(0, nrow=2, ncol=1), beta_update=dcoordinate_descent))
