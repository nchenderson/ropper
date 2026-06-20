\name{PPEP}
\alias{PPEP}
\title{
Population posterior expected percentiles
}
\description{
Computes estimates of cluster-level population posterior expected percentiles given an estimate of 
the regression coefficients and an estimate of the random effects variance.
}
\usage{
PPEP(beta.hat, Y, X, ses, tau.sq)
}
\arguments{
  \item{beta.hat}{
The 1 x p vector of regression coefficient estimates.
}
  \item{Y}{
The K x 1 vector of cluster-specific effect estimates.
}
  \item{X}{
The K x p design matrix containing cluster-specific covariates.
}
  \item{ses}{
The K x 1 vector sampling standard errors squared.
}
  \item{tau.sq}{
Estimate of the random effects variance. This must 
be a scalar greater than 0. 
}
}
\details{
The population posterior expected percentile for cluster k is defined
\eqn{E(\rho_{k} | Y_{k} ) = \Phi[ V_{k}( Y_{k} - X_{k}^{T}\hat{\beta})] }, where
\eqn{\Phi} is the cumulative distribution of a standard normal distribution, 
\eqn{\hat{\beta}} is the vector of regression coefficient estimates, and \eqn{V_{k}}
is defined as \eqn{V_{k} = \tau^{2}/[(\tau^{2} + ses_{k})\sqrt{2ses_{k} + \tau^{2}}] }.
The parameter \eqn{\rho_{k}} is the "population percentile" which is defined
as \eqn{\rho_{k} = \Phi(v_{k}/\tau)}, where \eqn{v_{k}} is the random effect
for the cluster-k effect estimate \eqn{Y_{k}}.
}
\value{
A K x 1 vector of population posterior expected percentile estimates.
}
\references{
Henderson, N. C., & Hartman, N. (2025). Targeted Parameter Estimation for Robust Empirical Bayes Ranking. 
\emph{arXiv preprint}. \url{https://arxiv.org/abs/2511.16530}
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
See also \code{\link{ropper}}, \code{\link{BLUP}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
\keyword{ regression }
\keyword{ robust }
\concept{ empirical Bayes }
\concept{ ranking }
