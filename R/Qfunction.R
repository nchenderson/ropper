Qfunction <- function(beta.coef, y, X, VV, H) {
  ### The ROPPER objective function to be minimized
  ## Inputs:
  ##.    beta.coef - vector of fixed-effects regression coefficients
  ##     y - length K vector of unit-specific estimates
  ##     X - K x p design matrix
  ##     V - length K vector of ranking shrinkage terms V_k
  ##     H - order of risk function approximation should be H = 1, H = 2, or H = 3
  
  resids <- as.numeric(y - X%*%beta.coef)
  if(H==1) {
    term1 <- sqrt(tau.sq/(2*pi))*mean(VV*dnorm(VV*resids))
  } else if(H==2) {
    der3 <- pnormder(VV*resids, ord=3)
    f1 <- sqrt(1/(2*pi))
    t1 <- sqrt(tau.sq)*f1*mean(VV*dnorm(VV*resids))
    term1 <- t1 - (1/6)*sqrt(tau.sq*tau.sq*tau.sq)*f1*mean((VV^3)*der3) - (1/6)*sqrt(tau.sq)*3*f1*mean(VV*dnorm(VV*resids))
  } else if(H==3) {
    der1 <- pnormder(VV*resids, ord=1)
    der3 <- pnormder(VV*resids, ord=3)
    der5 <- pnormder(VV*resids, ord=5)
    f1 <- sqrt(1/(2*pi))
    t1 <- sqrt(tau.sq)*f1*mean(VV*dnorm(VV*resids))
    t2 <-  t1 - (1/6)*sqrt(tau.sq*tau.sq*tau.sq)*f1*mean((VV^3)*der3) - (1/6)*sqrt(tau.sq)*3*f1*mean(VV*dnorm(VV*resids))
    term1 <- t2 + (1/40)*15*f1*sqrt(tau.sq)*mean(VV*der1) + (1/40)*10*f1*sqrt(tau.sq^3)*mean((VV^3)*der3) + (1/40)*f1*sqrt(tau.sq^5)*mean((VV^5)*der5)
  }
  midrank_discr <- pnorm(VV*resids) - 0.5
  #ans <- sum(midrank_discr*midrank_discr)
  ans <- -2*term1 + mean(midrank_discr*midrank_discr)
  return(ans)
}