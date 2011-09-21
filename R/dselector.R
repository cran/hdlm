dselector <- function(x,y,lambda){

  n <- nrow(x)
  p <- ncol(x)

  A <- t(x) %*% x
  R <- rbind(A, -A)
  a <- c(as.matrix(t(x) %*% y))
  r <- c(a-lambda, -a-lambda)
  beta <- quantreg::rq.fit.fnc(diag(p), rep(0,p), R=R, r=r)$coefficients
  return(beta)

}
