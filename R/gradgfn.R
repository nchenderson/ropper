gradgfn <- function(beta, Y, X, ses, tau.sq) {
  resids <- as.numeric(Y - X%*%beta)
  B <- tau.sq/(ses + tau.sq)
  VV <- sqrt(B/(2*ses + tau.sq))
  rr <- VV*resids
  
  num <- 2*dnorm(rr)*(pnorm(rr) - 1/2)
  den <- 1 - (pnorm(rr) - 0.5)^2
  fixed_grad <- num/den
  
  #ans <- crossprod(X, VV*fixed_grad)
  return(fixed_grad)
}