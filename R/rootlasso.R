# Dervied from matlab code by Alexendre Belloni
# x <- input matrix, y <- response vector, lambda <- tuning paramter
# MaxIter <- maximum number of interations
# OptTolNorm, OptTolObj <- tolorance of stopping criterion 
rootlasso <- function(x,y,lambda,MaxIter=1e5,OptTolNorm=1e-6,OptTolObj=1e-6) {

  # Derived inputs:
  n <- nrow(x)
  p <- ncol(x)
  gamma <- sqrt(diag(t(x) %*% x)/n)

  ## Start at ridge regression:
  RidgeMatrix <- diag(p)
  for(j in 1:p) {
    RidgeMatrix[j,j] <- lambda * gamma[j]
  }
  beta <- solve( t(x) %*% x + RidgeMatrix ) %*% (t(x) %*% y)

  ### Initialization:
  Iter <- 0                     

  XX <- (t(x) %*% x) / n 
  Xy <- (t(x) %*% y) / n

  ERROR = y - x %*% beta
  Qhat = sum(ERROR^2)/n

  ### Main loop:
  while(Iter < MaxIter) {
    Iter <- Iter + 1
    beta_old <- beta

    # Go over each coordinate, and optimize 
    for(j in 1:p) {
      # Compute the Shoot and Update the variable
      S0 <- XX[j,]%*%beta - XX[j,j]*beta[j] - Xy[j]
      if ( abs(beta[j]) > 0 ) {
        ERROR <- ERROR + x[,j] * beta[j]
        Qhat <- sum( ERROR^2 )/n
      }

      if( n^2 < ( lambda * gamma[j] )^2 / XX[j,j]) {
        beta[j] <- 0
      } else if(S0 > (lambda / n) * gamma[j] * sqrt(Qhat)) {
        ### Optimal beta(j) < 0
        # For lasso: beta(j) = (lambda - S0)/XX2(j,j); 
        # for square-root lasso
        beta[j] <- (  ( lambda * gamma[j] / sqrt( n^2 - (lambda * gamma[j])^2 / XX[j,j]  ) )  *
                    sqrt( max(Qhat - (S0^2/XX[j,j]),0) )  - S0 )   /   XX[j,j]

        ERROR <- ERROR - x[,j] * beta[j]
      } else if(S0 < - (lambda / n) * gamma[j] * sqrt(Qhat)) {
        ### Optimal beta(j) < 0
        # For lasso: beta(j) = (lambda - S0)/XX2(j,j); 
        # for square-root lasso
        beta[j] <- ( - ( lambda * gamma[j] / sqrt( n^2 - (lambda * gamma[j])^2 / XX[j,j]  ) )  *
                    sqrt( max(Qhat - (S0^2/XX[j,j]),0) )  - S0 )   /   XX[j,j]
        ERROR <- ERROR - x[,j] * beta[j]
      } else if(abs(S0) <= (lambda/n) * gamma[j] * sqrt(Qhat) ) {
        beta[j] <- 0
      }

    }

    # Update primal and dual value
    fobj <- sqrt( sum((x %*% beta-y)^2)/n )  +  t(lambda*gamma/n) %*% abs(beta)

    if( mean(ERROR^2) > 1e-10) {
      aaa <- sqrt(n)* ERROR / mean(ERROR^2)
      dual <- t(aaa) %*% y / n - t(abs( lambda * gamma / n - abs(t(x) %*% aaa)/n )) %*% abs(beta)
    } else {
      dual <- lambda * t(gamma) %*% abs(beta) / n
    } 

    # Stopping Crierion:
    if( sum (abs(beta - beta_old)) < OptTolNorm & fobj - dual < OptTolObj) break 

  }
  return(as.numeric(beta))
}
