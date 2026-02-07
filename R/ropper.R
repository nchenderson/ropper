ropper <- function(y, X, ses, tau.sq, H=1, method="optim") {
  ##############################################
  ## Inputs:
  ##   y - length K vector of unit-specific estimates
  ##   X - K x p design matrix
  ##   ses - length K vector of squared standard errors
  ##   H - order of risk function approximation should be H = 1, H = 2, or H = 3
  
  ## add estimation of tau.sq here.
  
  B <- tau.sq/(ses + tau.sq)
  VV <- sqrt(B/(2*ses + tau.sq))
  
  Qfn <- function(beta, H, VV) {
    resids <- as.numeric(y - X%*%beta)
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
  if(ncol(X)==1 & method=="optim") {
    if(H==1) {
      beta.rank <- optimize(Qfn, lower=-10,upper=10, H=1, VV=VV)$minimum
    } else if(H==2) {
      beta.rank <- optimize(Qfn, lower=-10,upper=10, H=2, VV=VV)$minimum
    } else if(H==3) {
      beta.rank <- optimize(Qfn, lower=-10,upper=10, H=3, VV=VV)$minimum
    }
  } else if(ncol(X) > 1 & method=="optim") {
    if(H==1) {
      beta.rank <- optim(rep(0, ncol(X)),fn=Qfn, H=1, VV=VV)$par
    } else if(H==2) {
      beta.rank <- optim(rep(0, ncol(X)),fn=Qfn, H=2, VV=VV)$par
    } else if(H==3) {
      beta.rank <- optim(rep(0, ncol(X)),fn=Qfn, H=3, VV=VV)$par
    }
  } else if(method=="MM") {
    
    lam <- sqrt(pi)/sqrt(2*tau.sq)
    afrac <- 1/3
    #V <- diag(VV)
    XV <- VV*X
    
    ## Use MLE as initial value of beta.
    ww.mle <- ses + tau.sq
    beta.old <- solve(crossprod(X, ww.mle*X), crossprod(X, ww.mle*y))
    
    niter <- 100
    tol <- 1e-6
    ObjFnVals <- rep(NA, niter + 1)
    BetaVals <- matrix(NA, nrow=niter + 1, ncol=ncol(X))
    ObjFnVals[1] <-  Qfn(beta.old, H=1, VV=VV)
    BetaVals[1,] <- beta.old
    
    for(k in 1:niter) {
      Xbeta.old <- as.numeric(X%*%beta.old)
      resids <- y - Xbeta.old
      rr <- VV*resids
      tmp1 <- dnorm(resids, sd=1/VV)
      tmp2 <- lam*(1 - (pnorm(rr) - 0.5)^2)
      w1 <- tmp1/(sum(tmp1) + sum(tmp2))
      w2 <- tmp2/(sum(tmp1) + sum(tmp2))
      
      wcomb <- w1 + w2*afrac
      XWV <- wcomb*VV*X
      XtVX <- crossprod(XV, XWV)
      dvec <- as.numeric(gradgfn(beta.old, y=y, X=X, ses=ses, tau.sq=tau.sq))
      
      vecterm1 <- crossprod(XV, w1*VV*y + w2*dvec + (w2*VV*Xbeta.old)*afrac)
      beta.new <- as.numeric(solve(XtVX, vecterm1))
      
      ObjFnVals[k+1] <- Qfn(beta=beta.new, H=1, VV=VV)
      BetaVals[k+1,] <- beta.new
      
      ss.parchange <- sqrt(sum((beta.new - beta.old)*(beta.new - beta.old)))
      if(ss.parchange < tol) {
        break
      }
      beta.old <- beta.new
    }
    ObjFnVals <- ObjFnVals[!is.na(ObjFnVals)]
  }
  if(method=="optim") {
    return(list(par=beta.rank, objfn=NULL))
  } else {
    return(list(par=beta.new, objfn=ObjFnVals))
  }
}

