Q_REML <- function(tau.sq, y, X, ses) {
  ######################################################################
  ##
  ## Input:
  ##    tau.sq - the value of tau.sq (tau.sq = var(vi))
  ##    Y - the vector (of length K) of responses
  ##    X - the K x p design matrix
  ##    ses - var(Y_k|\theta_k) = ses_k
  ##
  ##  Output:
  ##    the value of the objective function used in REML estimation
  ########################################################################
  
  
  mvar <- ses + tau.sq  ## marginal variances
  inv_mvar <- 1/mvar
  XtDX <- crossprod(X, X*inv_mvar)
  XtDy <- crossprod(X, y*inv_mvar)
  Py <- y - X%*%solve(XtDX, XtDy)
  
  t1 <- sum(log(mvar))
  t2 <- determinant(XtDX, logarithm=TRUE)$modulus
  t3 <- sum(y*inv_mvar*Py)
  ans <- t1 + t2 + t3
  return(ans)
}
