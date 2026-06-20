ropper <- function(Y, X, ses, tau.sq = c("REML", "kNN"), H=1, 
                   opt.method=c("MM", "optim"), trim.quant=0.0,
                   control = list()) {
  ##############################################
  ## Inputs:
  ##   y - length K vector of unit-specific estimates
  ##   X - K x p design matrix
  ##   ses - length K vector of squared standard errors
  ##   tau.sq - method for estimating tau.sq or an estimate of tau.sq
  ##   H - order of risk function approximation should be H = 1, H = 2, or H = 3
  
  ## Get values of control parameters maxiter and tol:
  control.default <- list(maxiter = 500, tol = 1e-06)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  maxiter <- control$maxiter
  tol <- control$tol
  
  ## Get optimization method and method for estimating tau.sq:
  opt.method <- match.arg(opt.method)
  
  if(opt.method != "MM" & opt.method != "optim") {
       stop("The opt.method argument should equal MM or optim")
  }
  if(class(tau.sq) == "numeric") {
    if(tau.sq <= 0.0) {
        stop("The tau.sq argument should either be a positive number or equal
                to either REML or kNN")
    }
  }
  if(class(tau.sq) != "numeric") {
     tau.sq <- match.arg(tau.sq)
     if(tau.sq != "REML" & tau.sq !="kNN") {
         stop("The tau.sq argument should either be a positive number or equal
                to either REML or kNN")
     }
     if(tau.sq=="REML") {
        ## Did we use winsorization for the REML estimate?
        tau.sq <- REMLTausqEst(Y=Y, X=X, ses=ses, trim.quant=trim.quant) 
     } else if(tau.sq=="kNN") {
        tau.sq <- KnnTausqEst(Y=Y, X=X, ses=ses, trim.quant=trim.quant)
     }
  }
  
  B <- tau.sq/(ses + tau.sq)
  VV <- sqrt(B/(2*ses + tau.sq))
 
  if(ncol(X)==1 & opt.method=="optim") {
    if(H==1) {
      beta.rank <- optimize(Qfunction, lower=-10, upper=10, y=Y, X=X,
                            VV=VV, tau.sq=tau.sq, H=1)$minimum
    } else if(H==2) {
      beta.rank <- optimize(Qfunction, lower=-10,upper=10, y=Y, X=X,
                            VV=VV, tau.sq=tau.sq, H=2)$minimum
    } else if(H==3) {
      beta.rank <- optimize(Qfunction, lower=-10,upper=10, y=Y, X=X,
                            VV=VV, tau.sq=tau.sq, H=3)$minimum
    }
  } else if(ncol(X) > 1 & opt.method=="optim") {
    if(H==1) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=Y, X=X,
                         VV=VV, tau.sq=tau.sq, H=1)$par
    } else if(H==2) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=Y, X=X,
                         VV=VV, tau.sq=tau.sq, H=2)$par
    } else if(H==3) {
      beta.rank <- optim(rep(0, ncol(X)), fn=Qfunction, y=Y, X=X,
                         VV=VV, tau.sq=tau.sq, H=3)$par
    }
  } else if(opt.method=="MM") {
    
      lam <- sqrt(pi)/sqrt(2*tau.sq)
      afrac <- 1/3
      XV <- VV*X
    
      ## Use MLE as initial value of beta.
      ww.mle <- ses + tau.sq
      beta.old <- solve(crossprod(X, ww.mle*X), crossprod(X, ww.mle*Y))
    
      ObjFnVals <- rep(NA, maxiter + 1)
      BetaVals <- matrix(NA, nrow=maxiter + 1, ncol=ncol(X))
      ObjFnVals[1] <-  Qfunction(beta.coef=beta.old, Y=Y, X=X, VV=VV, 
                                 tau.sq=tau.sq, H=1)
      BetaVals[1,] <- beta.old
      for(k in 1:maxiter) {
          Xbeta.old <- as.numeric(X%*%beta.old)
          resids <- Y - Xbeta.old
          rr <- VV*resids
          tmp1 <- dnorm(resids, sd=1/VV)
          tmp2 <- lam*(1 - (pnorm(rr) - 0.5)^2)
          w1 <- tmp1/(sum(tmp1) + sum(tmp2))
          w2 <- tmp2/(sum(tmp1) + sum(tmp2))
      
          wcomb <- w1 + w2*afrac
          XWV <- wcomb*VV*X
          XtVX <- crossprod(XV, XWV)
          dvec <- as.numeric(gradgfn(beta.old, Y=Y, X=X, ses=ses, tau.sq=tau.sq))
      
          vecterm1 <- crossprod(XV, w1*VV*Y + w2*dvec + (w2*VV*Xbeta.old)*afrac)
          beta.new <- as.numeric(solve(XtVX, vecterm1))

          ObjFnVals[k+1] <- Qfunction(beta.coef=beta.new, Y=Y, X=X, VV=VV, 
                                      tau.sq=tau.sq, H=1)
          BetaVals[k+1,] <- beta.new
      
          ## look at sum of squares of parameter changes to determine convergence 
          ss.parchange <- sqrt(sum((beta.new - beta.old)*(beta.new - beta.old)))
          if(ss.parchange < tol) {
             break
          }
          beta.old <- beta.new
      }
      beta.rank <- beta.new
      ObjFnVals <- ObjFnVals[!is.na(ObjFnVals)]
  }
  ## Compute final residuals
  resids <- as.numeric(Y - X%*%beta.rank)
  ## Compute vector of population posterior expected ranks.
  post.rank <- pnorm(VV*resids)
  ## Should add optimal percentiles to returned list.
  if(opt.method=="optim") {
    return(list(coefficients=beta.rank, ppep=post.rank, objfn=NULL))
  } else {
    return(list(coefficients=beta.rank, ppep=post.rank, objfn=ObjFnVals))
  }
}

