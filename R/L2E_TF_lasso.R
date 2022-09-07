#' Solution path of the L2E trend filtering regression with Lasso
#'
#' \code{L2E_TF_lasso} computes the solution path of the robust trend filtering regression under the L2 criterion with Lasso penalty
#'
#' @param y Response vector
#' @param X Design matrix. Default is the identity matrix.
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param D The fusion matrix
#' @param lambdaSeq A decreasing sequence of values for the tuning parameter lambda, can be omitted
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (matrix) and
#' tau (vector) for each value of the tuning parameter lambda,
#' the run time (vector) for each lambda,
#' and the sequence of lambda used in the regression (vector)
#' @importFrom stats mad
#' @export
#' @examples
#' ## Completes in 10 seconds
#'
#' set.seed(12345)
#' n <- 100
#' x <- 1:n
#' f <- matrix(rep(c(-2,5,0,-10), each=n/4), ncol=1)
#' y <- y0 <- f + rnorm(length(f))
#'
#' ## Clean Data
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#'
#' D <- myGetDkn(1, n)
#' lambda <- 10^seq(-1, -2, length.out=20)
#' sol <- L2E_TF_lasso(y=y, D=D, lambdaSeq=lambda)
#'
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' lines(x, sol$Beta[,1], col='blue', lwd=3) ## 1st lambda
#'
#' ## Contaminated Data
#' ix <- sample(1:n, 10)
#' y[ix] <- y0[ix] + 2
#'
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#'
#' sol <- L2E_TF_lasso(y=y, D=D, lambdaSeq=lambda)
#'
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' lines(x, sol$Beta[,1], col='blue', lwd=3) ## 1st lambda
#'
L2E_TF_lasso <- function(y,X,beta0,tau0,D,lambdaSeq,max_iter=1e2,tol=1e-4,Show.Time=TRUE){

  if(missing(X)){
    X <- diag(nrow = length(y))  # initial X
  }

  if(missing(beta0)){
    beta0 <- rep(mean(y), ncol(X))  # initial beta
  }

  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }

  if (tau0 <= 0) stop("Entered non-positive tau0")

  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of lambda
  }


  Nlambda <- length(lambdaSeq)

  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)

  time <- proc.time()
  for (j in 1:Nlambda) {

    res <- l2e_regression_TF_lasso(y, X, beta=beta0, tau=tau0, D=D, lambda=lambdaSeq[j],
                                 max_iter=max_iter, tol=tol, Show.Time = FALSE)
    Beta[, j] <- beta0 <- res$beta
    Tau[j] <- res$tau # warm start tau is not good

  }
  runtime <-  proc.time() - time
  if(Show.Time) print(runtime)

  return(list(Beta=Beta, Tau=Tau, runtime=runtime, lambdaSeq=lambdaSeq))

}
