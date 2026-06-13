\name{KnnTausqEst}
\alias{KnnTausqEst}
\title{
Nearest neighbors estimation of random effects variance
}
\description{
Computes an estimate of the random effects variance using a nonparametric, 
k nearest neighbors approach. In this approach, the standard
errors for the cluster-level effect estimates are assumed to be fixed.
}
\usage{
KnnTausqEst(Y, X, ses, trim.quant = 0, k = 1)
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
\item{k}{
The number of neighbors to include in the k nearest neighbors procedure. This should
be an integer both greater than or equal to 1 and less than or equal to K.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The estimate of the random effects variance. A vector of length 1.
}
\references{
Devroye, L., Györfi, L., Lugosi, G., & Walk, H. (2018). A nearest neighbor 
estimate of the residual variance. Electronic Journal of Statistics, 12(1), 
1752–1778. \doi{10.1214/18-EJS1438}
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\note{
%%  ~~further notes~~
}
\seealso{
See also \code{\link{REMLTausqEst}}, \code{\link{ropper}}
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
\keyword{ regression }
\keyword{ robust }
