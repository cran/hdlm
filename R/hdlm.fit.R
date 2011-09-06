hdlm.fit <-
function(x, y, N, C, ...) {

  p <- ncol(x)
  n <- nrow(x)
  # Estimate sigma: (can be much more complicated)
  if(is.null(C)) C <- floor(min(c(sqrt(n), p)))
  if(round(C) != C || C <= 0) {
    sigma_hat <- var(y)
  } else {
    sigma_hat <- min(var(y), summary(lm(y ~ x[,sample(1:p,C)]))$sigma)
  }
  lambda <- sigma_hat * (log(p) + 1) / log(n)

  hdlm.stat <- function(x, INDEX) {
    # Estimate model:
    out <- lars(x[INDEX,],y[INDEX])
    res <- NULL
    if(lambda > out$lambda[1]) {
      res <- rep(0,p)
    }
    if(lambda < out$lambda[length(out$lambda)]) {
      res <- out$beta[nrow(out$beta),]
    }
    if(is.null(res)) {
      i <- max(which(out$lambda > lambda))
      betaHat <- out$beta[i:(i+1),]
      lam <- out$lambda[(i+1):i]
      alpha <- 1-(lambda - lam[1])/(lam[2] - lam[1])
      res <- betaHat[1,]*(1-alpha) + betaHat[2,]*alpha
    }
    return(res)
  }

  # Bootstrap to determine standard deviation:
  B  <- boot(x, hdlm.stat, N)

  # Calculate p-values and standard errors
  hdlm.pval <- function(v) {
    1-mean(v[-1] == v[1])*abs(v[1])
  }

  se <- apply(B$t,2,sd)
  bias <- apply(B$t,2,mean) - B$t0
  pvalue <- apply(sign(rbind(B$t0, B$t)), 2, hdlm.pval)

  # Standard estimators:
  point_estimator <- B$t0
  fitted <- x %*% point_estimator
  resid <- fitted - y


  z <- list(coefficients=point_estimator, residuals=resid, effects=NULL, rank=c(n,p),
            fitted.values=fitted, assign=NULL, standard.error=se, bias=bias, p.value=pvalue,
            sigma.hat = sigma_hat, N=N)

  return(z)

}

