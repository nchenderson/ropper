REMLTausqEst <- function(Y, X, ses, trim.quant=0.0) {
  if(trim.quant < 0 | trim.quant >= 0.5) {
       stop("trim.quant should be a number greater than or equal to 0 and less than 0.5")
  }
  y.winsor <- Y
  if(trim.quant > 0) {
      qq <- quantile(Y, probs = c(trim.quant, 1 - trim.quant))
      
      # Replace values below the lower bound
      y.winsor[Y < qq[1]] <- qq[1]
      # Replace values above the upper bound
      y.winsor[Y > qq[2]] <- qq[2]
  } 
  VY <- var(y.winsor)
  tausq_hat <- optimize(Q_REML, interval=c(0,4*VY), Y=y.winsor, X=X, ses=ses)$minimum
  return(tausq_hat)
}