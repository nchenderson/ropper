\name{REMLTausqEst}
\alias{REMLTausqEst}
\title{
Restricted maximum likelihood estimation of random effects variance
}
\description{
Computes an estimate of the random effects variance by maximizing
the restricted maximum likelihood objective function. The standard
errors for the cluster-level effect estimates are assumed to be fixed.
}
\usage{
REMLTausqEst(Y, X, ses, trim.quant = 0)
}
\arguments{
  \item{Y}{
The K x 1 vector of cluster-specific effect estimates.
}
  \item{X}{
The K x p design matrix containing cluster-specific covariates.
}
  \item{ses}{
The K x 1 vector sampling standard errors squared.
}
  \item{trim.quant}{
The proportion of extreme estimates to exclude when estimating the random effects
variance. Cluster-specific estimates lower than the \code{trim.quant} quantile will be 
set equal to the \code{trim.quant} quantile, and estimates larger than the
1 - \code{trim.quant} quantile will be set equal to the 1 - \code{trim.quant} quantile.
This must be a number greater than or equal to 0 and less than 0.5.
}
}
\details{
Estimate of the random effects variance \eqn{\tau^{2}} is found by maximizing
the restricted maximum likelihood objective function 
\eqn{-\log\det(\mathbf{W}) - \log\det(\mathbf{X}^{T}\mathbf{W}^{-1}\mathbf{X}) + \mathbf{Y}^{T}\mathbf{P}\mathbf{Y} }
with respect to \eqn{\tau^{2}}. Here, \eqn{\mathbf{Y}} is the K x 1 vector of 
cluster-level effect estimates, \eqn{\mathbf{W}} is the K x K diagonal matrix
with diagonal elements \eqn{\tau^{2} + ses_{k}^{2}}, and \eqn{\mathbf{P}} is
the K x K matrix defined as 
\eqn{\mathbf{P} = \mathbf{W}^{-1} - \mathbf{W}^{-1}\mathbf{X}(\mathbf{X}^{T}\mathbf{W}^{-1}\mathbf{X})^{-1}\mathbf{X}^{T}\mathbf{W}^{-1}}
}
\value{
The estimate of the random effects variance. A vector of length 1.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
See also \code{\link{ropper}},
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
\keyword{ optimize }
\keyword{ regression }
