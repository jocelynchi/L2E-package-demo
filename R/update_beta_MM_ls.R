#' Beta update in L2E multivariate regression - MM
#'
#' \code{update_beta_MM_ls} updates beta for L2E multivariate regression using MM
#'
#' @param y Response
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Precision estimate
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#' @importFrom stats lm
#'
update_beta_MM_ls <- function(y,X,beta,tau,max_iter=1e2,tol=1e-4) {

  n <- length(y)
  for (i in 1:max_iter) {

    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))

    beta <- lm(y~X-1, weights = sqrt(w))$coefficients
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  return(list(beta=beta,iter=i))
}
