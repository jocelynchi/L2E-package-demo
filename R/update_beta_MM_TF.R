#' Beta update in L2E trend filtering regression - MM
#'
#' \code{update_beta_MM_TF} updates beta in L2E trend filtering regression using the distance penalty
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param D The fusion matrix
#' @param k The number of nonzero entries in D*beta
#' @param rho The parameter in the proximal distance algorithm
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#'
update_beta_MM_TF <- function(y,X,beta,tau,D,k,rho,max_iter=1e2,tol=1e-4) {

  n <- nrow(X)

  for (i in 1:max_iter) {
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))


    beta <- gradient_descent(y, X, w, beta_last, D, rho, k)$beta

    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }

  return(list(beta=beta,iter=i))
}










## Project x onto the set C={x| x has at most k zero entries}
#' @importFrom utils tail
Proj_sparse <- function(x, k){

  res <- sort(abs(x), method="quick", index.return=TRUE)
  ind <- tail(res$ix, k)
  xnew <- rep(0, length(x))
  xnew[ind] <- x[ind]
  return(xnew)
}

f <- function(X, y, beta, k, rho){
  pjbeta <- Proj_sparse(beta, k)
  s1 <- rho*norm(beta - pjbeta, "2")^2/2
  Xbeta <- X%*%beta
  s2 <- norm(y-Xbeta, "2")^2/2
  return(s1+s2)
}


#' @importFrom Matrix Diagonal
gradient_descent <- function(y, X, w, beta, D, rho, k, max_iter=1e2, tol=1e-5){
  n <- length(w)
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X

  for (i in 1:max_iter) {
    beta_last <- beta
    a <- XtWX%*%beta_last
    Dbeta <- as.vector(D%*%beta_last)
    Dbeta_proj <- Proj_sparse(Dbeta, k)

    gradient <- as.vector(a - XtWy+rho*t(D)%*%(Dbeta- Dbeta_proj))

    Av <- XtWX%*%gradient
    vAv <- t(gradient)%*%Av
    Dv <- as.vector(D%*%gradient)

    stepsize <- as.numeric(norm(gradient, "2")^2/(vAv + rho*norm(Dv, "2")^2))
    beta <- beta_last - stepsize*gradient

    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }

  return(list(beta=beta))
}
