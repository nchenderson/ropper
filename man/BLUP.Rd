\name{BLUP}
\alias{BLUP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior expected means
}
\description{
Computes estimates of cluster-level posterior expected means given an estimate of 
the regression coefficients and an estimate of the random effects variance. For
a fixed value of the random effects variance, the collection of posterior expected
mean estimates is commonly referred to as the "best linear unbiased predictor",
or BLUP.
}
\usage{
BLUP(beta.hat, Y, X, ses, tau.sq)
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
The BLUP for cluster k is defined as
\eqn{E(v_{k} | Y_{k} ) = B_{k}( Y_{k} - X_{k}^{T}\hat{\beta}) }, where
\eqn{\hat{\beta}} is the vector of regression coefficient estimates and \eqn{V_{k}}
is defined as \eqn{B_{k} = \tau^{2}/(\tau^{2} + ses_{k}) }
}
\value{
A K x 1 vector of posterior mean estimates.
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{PPEP}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{ regression }
\concept{ empirical Bayes }
\concept{ shrinkage }
