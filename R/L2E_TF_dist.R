#' Solution path of the L2E trend filtering regression with distance penalization
#'
#' \code{L2E_TF_dist} computes the solution path of the robust trend filtering regression under the L2 criterion with distance penalty
#'
#' @param y Response vector
#' @param X Design matrix. Default is the identity matrix.
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param D The fusion matrix
#' @param kSeq A sequence of tuning parameter k, the number of nonzero entries in D*beta
#' @param rhoSeq An increasing sequence of tuning parameter rho, can be omitted
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (matrix) and
#' tau (vector) for each value of the tuning parameter k,
#' the path of estimates for beta (list of matrices) and tau (matrix) for each value of rho,
#' the run time (vector) for each k,
#' and the sequence of rho and k used in the regression (vectors)
#' @importFrom stats mad
#' @importFrom signal filter
#' @importFrom signal MedianFilter
#' @export
#' @examples
#' ## Completes in 15 seconds
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
#' k <- c(4,3,2)
#' rho <- 10^8
#' # (not run)
#' # sol <- L2E_TF_dist(y=y, D=D, kSeq=k, rhoSeq=rho)
#'
#' # plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' # lines(x, f, lwd=3)
#' # lines(x, sol$Beta[,3], col='blue', lwd=3) ## k=2
#' # lines(x, sol$Beta[,2], col='red', lwd=3) ## k=3
#' # lines(x, sol$Beta[,1], col='dark green', lwd=3) ## k=4
#'
#' ## Contaminated Data
#' ix <- sample(1:n, 10)
#' y[ix] <- y0[ix] + 2
#'
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#'
#' # (not run)
#' # sol <- L2E_TF_dist(y=y, D=D, kSeq=k, rhoSeq=rho)
#'
#' # plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' # lines(x, f, lwd=3)
#' # lines(x, sol$Beta[,3], col='blue', lwd=3) ## k=2
#' # lines(x, sol$Beta[,2], col='red', lwd=3) ## k=3
#' # lines(x, sol$Beta[,1], col='dark green', lwd=3) ## k=4
#'
L2E_TF_dist <- function(y,X,beta0,tau0,D,kSeq,rhoSeq,max_iter=1e2,tol=1e-4,Show.Time=TRUE) {

  if(missing(X)){
    X <- diag(nrow = length(y))  # initial X
  }

  if(missing(beta0)){
    beta0 <- filter(MedianFilter(9), y)  # initial beta, random initial is very bad
  }

  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }

  if (tau0 <= 0) stop("Entered non-positive tau0")

  if(missing(rhoSeq)){
    rhoSeq <- 10^seq(0, 4, length.out = 20)  # set a sequence of rho
  }


  Nk <- length(kSeq)
  Nrho <- length(rhoSeq)

  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nk)
  Beta_path <- list()
  Tau <- double(Nk)
  Tau_path <- matrix(0, Nk, Nrho)
  time <- double(Nk)

  for (i in 1:Nk) {

    k <- kSeq[i]
    beta <- matrix(0, nrow = ncol(X), ncol = Nrho)
    tau <- double(Nrho)

    start_time <- proc.time()
    for (j in 1:Nrho) {

      res <- l2e_regression_TF_dist(y, X, beta = beta0, tau = tau0, D=D, k=k, rho=rhoSeq[j],
                                      max_iter=max_iter, tol=tol, Show.Time = FALSE)
      beta[, j] <- beta0 <- res$beta
      tau[j] <- tau0 <- res$tau # but warm start tau here helps reduce running time

    }
    runtime <-  proc.time() - start_time
    if(Show.Time) print(runtime)

    Beta_path[[i]] <- beta
    Beta[, i] <- beta0
    Tau_path[i, ] <- tau
    Tau[i] <- tau0
    time[i] <- runtime[[3]]
  }


  return(list(Beta=Beta, Beta_path=Beta_path, Tau=Tau, Tau_path=Tau_path, time = time,
              rhoSeq = rhoSeq, kSeq = kSeq))

}
