pnormder <- function(x, ord) {
  ## R function to compute up to 5th derivative
  ## of the pnorm function
  p1 <- dnorm(x)
  p2 <- (-x)*dnorm(x)
  p3 <- (-x)*p2 - p1
  p4 <- (-x)*p3 - 2*p2
  p5 <- (-x)*p4 - 3*p2
  pvec <- c(p1, p2, p3, p4, p5)
  return(pvec[ord])
}