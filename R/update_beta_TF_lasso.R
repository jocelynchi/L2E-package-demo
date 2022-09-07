#' Beta update in L2E trend filtering regression using Lasso
#'
#' \code{update_beta_TF_lasso} updates beta in L2E trend filtering regression using the Lasso penalty
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param D The fusion matrix
#' @param lambda The tuning parameter
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#'
update_beta_TF_lasso <- function(y,X,beta,tau,D,lambda,max_iter=1e2,tol=1e-4) {

  n <- nrow(X)
  for (i in 1:max_iter) {

    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))

    beta <- sol_TF_lasso(y, X, w, D, lambda)$beta


    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }

  return(list(beta=beta,iter=i))
}







## QP to solve trend filtering lasso problem
#' @importFrom Matrix Diagonal
#' @importFrom Matrix Matrix
#' @importFrom osqp osqpSettings
#' @importFrom osqp osqp
sol_TF_lasso <- function(y, X, w, D, lambda) {

  n <- length(y)
  p <- ncol(X)
  r <- nrow(D)

  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X

  P <- Matrix(0, nrow = p+2*r, ncol = p+2*r)
  P[1:p, 1:p] <- XtWX

  q <- c(as.vector(-XtWy), rep(lambda, 2*r))

  A <- rbind(cbind(D, Diagonal(n=r, x=-1), Diagonal(n=r, x=1)),
             cbind(-D, Diagonal(n=r, x=1), Diagonal(n=r, x=-1)),
             cbind(Matrix(0, nrow=r, ncol = p,sparse = TRUE), Diagonal(n=r, x=1), Diagonal(n=r, x=0)),
             cbind(Matrix(0, nrow=r, ncol = p,sparse = TRUE), Diagonal(n=r, x=0), Diagonal(n=r, x=1)))
  l <- rep(0,  4*r)
  u <- rep(Inf, 4*r)

  settings <- osqpSettings(verbose = FALSE)

  model  <- osqp(P=P, q=q, A=A, l=l, u=u, settings)

  result <- model$Solve()

  beta <- result$x[1:p]

  v <- result$x[(p+1):(p+r)]-result$x[(p+r+1):(p+2*r)]

  return(list(beta = beta, v = v))

}
