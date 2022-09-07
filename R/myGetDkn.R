#' Compute kth order differencing matrix
#'
#' \code{myGetDkn} computes the kth order differencing matrix for use in
#' trend filtering regression
#'
#' @param k Order of the differencing matrix
#' @param n Number of time points
#' @return Returns a Matrix object as the kth order differencing matrix
#' @importFrom Matrix Matrix
#' @export
myGetDkn <- function(k, n) {

  D1n<-function(n){         # a function generates 1-st order difference matrix with n columns
    D1n=matrix(0, nrow = (n-1), ncol = n)
    for (i in 1:(n-1)) {
      D1n[i,i]=-1
      D1n[i,i+1]=1
    }
    return(D1n)
  }
  if(k==1){
    return(D1n(n))
  } else{
    Dkn=D1n(n)
    for (i in 1:(k-1)) {
      m=n-i
      Dkn=D1n(m)%*%Dkn
    }
    A=Matrix(Dkn, sparse=TRUE)
    return(A)
  }
}
