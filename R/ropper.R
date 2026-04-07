ropper <- function(y, X, ses, tau.sq = c("reml", "kNN"), H=1, 
                   method=c("MM", "optim"), control = list()) {
  ##############################################
  ## Inputs:
  ##   y - length K vector of unit-specific estimates
  ##   X - K x p design matrix
  ##   ses - length K vector of squared standard errors
  ##   tau.sq - method for estimating tau.sq or an estimate of tau.sq
  ##   H - order of risk function approximation should be H = 1, H = 2, or H = 3
  
  ## Get values of control parameters maxiter and tol:
  control.default <- list(maxiter = 100, tol = 1e-06)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  maxiter <- control$maxiter
  tol <- control$tol
  
  ## Get optimization method and method for estimating tau.sq:
  method <- match.arg(method)
  tau.sq <- match.arg(tau.sq)
  
  if(method != "MM" | method != "optim") {
       stop("The method argument should equal mm or optim")
  }
  
  if(tau.sq=="reml") {
      
  } else if(tau.sq=="kNN") {
      
  }
  if(class(tau.sq) != "numeric") {
       stop("The tau.sq argument should either be a positive number or equal
            to either reml or kNN")
       if(tau.sq <= 0.0) {
           stop("The tau.sq argument should either be a positive number or equal
                to either reml or kNN")
       }
  }
  
  B <- tau.sq/(ses + tau.sq)
  VV <- sqrt(B/(2*ses + tau.sq))
 
  if(ncol(X)==1 & method=="optim") {
    if(H==1) {
      beta.rank <- optimize(Qfunction, lower=-10, upper=10, y=y, X=X,
                            ses=ses, tau.sq=tau.sq, H=1, VV=VV)$minimum
    } else if(H==2) {
      beta.rank <- optimize(Qfn, lower=-10,upper=10, y=y, X=X,
                            ses=ses, tau.sq=tau.sq, H=2, VV=VV)$minimum
    } else if(H==3) {
      beta.rank <- optimize(Qfn, lower=-10,upper=10, y=y, X=X,
                            ses=ses, tau.sq=tau.sq, H=3, VV=VV)$minimum
    }
  } else if(ncol(X) > 1 & method=="optim") {
    if(H==1) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=y, X=X,
                         ses=ses, tau.sq=tau.sq, H=1, VV=VV)$par
    } else if(H==2) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=y, X=X,
                         ses=ses, tau.sq=tau.sq, H=2, VV=VV)$par
    } else if(H==3) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=y, X=X,
                         ses=ses, tau.sq=tau.sq, H=3, VV=VV)$par
    }
  } else if(method=="MM") {
    
      lam <- sqrt(pi)/sqrt(2*tau.sq)
      afrac <- 1/3
      XV <- VV*X
    
      ## Use MLE as initial value of beta.
      ww.mle <- ses + tau.sq
      beta.old <- solve(crossprod(X, ww.mle*X), crossprod(X, ww.mle*y))
    
      ObjFnVals <- rep(NA, maxiter + 1)
      BetaVals <- matrix(NA, nrow=maxiter + 1, ncol=ncol(X))
      ObjFnVals[1] <-  Qfunction(beta.old, y=y, X=X, ses=ses, tau.sq=tau.sq,
                                 H=1)
      BetaVals[1,] <- beta.old
      Qfunction <- function(beta.coef, y, X, ses, tau.sq, H)
      for(k in 1:maxiter) {
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
      
          ObjFnVals[k+1] <- Qfunction(beta=beta.new, y=y, X=X, ses=ses, 
                                      tau.sq=tau.sq, H=1)
          BetaVals[k+1,] <- beta.new
      
          ## look at sum of squares of parameter changes to determine convergence 
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

