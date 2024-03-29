% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_beta_MM_TF.R
\name{update_beta_MM_TF}
\alias{update_beta_MM_TF}
\title{Beta update in L2E trend filtering regression - MM}
\usage{
update_beta_MM_TF(y, X, beta, tau, D, k, rho, max_iter = 100, tol = 1e-04)
}
\arguments{
\item{y}{Response vector}

\item{X}{Design matrix}

\item{beta}{Initial vector of regression coefficients}

\item{tau}{Initial precision estimate}

\item{D}{The fusion matrix}

\item{k}{The number of nonzero entries in D*beta}

\item{rho}{The parameter in the proximal distance algorithm}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Relative tolerance}
}
\value{
Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
}
\description{
\code{update_beta_MM_TF} updates beta in L2E trend filtering regression using the distance penalty
}
