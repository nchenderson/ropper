PEPP <- function(beta.hat, Y, X, ses, tau.sq) {
  ## Function to compute the vector of posterior expected 
  ## population percentiles (PEPP)
  
  ## Compute residuals
  resids <- as.numeric(Y - X%*%beta.hat)
  ## Compute vector of ranking "shrinkage" terms
  B <- tau.sq/(ses + tau.sq)
  VV <- sqrt(B/(2*ses + tau.sq))
  ## Compute vector of population posterior expected ranks.
  post.perc <- pnorm(VV*resids)
  return(post.perc)
}