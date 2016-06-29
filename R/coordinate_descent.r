update_coordinates = function(X, W, z, lambda, alpha, beta) {
  beta_old = beta
  for (i in 1:length(beta)) {
    beta[i] = soft_thresh_r(sum(W*X[,i]*(z - X[,-i, drop=FALSE] %*% beta_old[-i])),
                               sum(W)*lambda*alpha)
  }
  beta / (colSums(W*X^2) + lambda*(1-alpha))
}

quadratic_loss = function(X, W=1, z, lambda, alpha, beta) {
  1/nrow(X)/2 * sum(W*(z - X %*% beta)^2) - 
      lambda * ((1-alpha) * sum(beta^2)/2 + alpha * sum(abs(beta)))
}

coordinate_descent = function(X, W, z, lambda, alpha, beta, maxit) {
  quad_loss = quadratic_loss(X, W, z, lambda, alpha, beta)
  for(i in 1:maxit) {
    beta_old = beta
    quad_loss_old = quad_loss
    beta = update_coordinates(X, W, z, lambda, alpha, beta)
    quad_loss = quadratic_loss(X, W, z, lambda, alpha, beta)
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
