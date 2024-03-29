% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_beta_TF_lasso.R
\name{update_beta_TF_lasso}
\alias{update_beta_TF_lasso}
\title{Beta update in L2E trend filtering regression using Lasso}
\usage{
update_beta_TF_lasso(y, X, beta, tau, D, lambda, max_iter = 100, tol = 1e-04)
}
\arguments{
\item{y}{Response vector}

\item{X}{Design matrix}

\item{beta}{Initial vector of regression coefficients}

\item{tau}{Initial precision estimate}

\item{D}{The fusion matrix}

\item{lambda}{The tuning parameter}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Relative tolerance}
}
\value{
Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
}
\description{
\code{update_beta_TF_lasso} updates beta in L2E trend filtering regression using the Lasso penalty
}
