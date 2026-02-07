BLUP <- function(beta.hat, Y, X, ses, tau.sq) {
  ## Function to compute the posterior expected means (BLUPs)
  
  ## Compute residuals
  resids <- as.numeric(Y - X%*%beta.hat)
  
  ## Compute vector of shrinkage terms and vector of BLUPs
  B <- tau.sq/(ses + tau.sq)
  blup <- B*resids
  return(blup)
}