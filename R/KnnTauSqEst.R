KnnTausqEst <- function(Y, X, ses, trim.quant=0, k = 1) {
  ### Note that this function will not work for the case
  ### when X only has a single intercept column
  
  ### This estimator is based on the idea of
  ## "A nearest neighbor estimate of the
  ## residual variance". Devroye Electronic J of Stat (2018)

  n <- nrow(X)
  if(n%%2 == 0) {
    half_samp <- n/2
  } else {
    half_samp <- (n+1)/2
  }
  y.winsor <- Y
  if(trim.quant > 0) {
    qq <- quantile(Y, probs = c(trim.quant, 1 - trim.quant))
    
    # Replace values below the lower bound
    y.winsor[Y < qq[1]] <- qq[1]
    # Replace values above the upper bound
    y.winsor[Y > qq[2]] <- qq[2]
  } 
  ind <- sample(1:n, size=half_samp)
  ind.train <- setdiff(1:n, ind)
  y.test <- y.winsor[ind]
  y.train <- y.winsor[ind.train]
  X.test <- X[ind,]
  X.train <- X[ind.train,]
  ses.test <- ses[ind]
  n.test <- nrow(X.test)
  n.train <- nrow(X.train)
  muhat <- dists <- rep(NA, n.test)

  for (i in 1:n.test) {
    for(j in 1:n.train) {
      # Compute distances to all training points
      dists[j] <- sum((X.test[i,] - X.train[j,])*(X.test[i,] - X.train[j,]))
      # Find indices of k nearest neighbors
    }
    nn_idx <- order(dists)[1:(k+1)]
    nn_idx <- nn_idx[-1]
    # Average the y values of the neighbors (not including i)
    muhat[i] <- mean(y.train[nn_idx])
  }
  ## With n = 50, probs=c(0.025, 0.975) is a reasonable choice
  
  Sstar <- mean(y.test*muhat)
  Lstar <- mean(y.winsor*y.winsor) - Sstar
  tausq_hat <- max(Lstar - mean(ses.test), 0)
  
  return(tausq_hat)
}
