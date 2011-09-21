hdlm.fit <-
function(x, y, method=c('mc+', 'scad', 'lasso', 'root-lasso', 'dantzig'),
         p.value.method = c('one-split', 'multi-split', 'bootstrap', 'none'),
         N, C, sigma.hat = NULL, ...) {

  method.num <- match(method,c('mc+', 'scad', 'lasso', 'root-lasso', 'dantzig'))
  p <- ncol(x)
  n <- nrow(x)
  # Estimate sigma:
  if(is.null(C)) C <- floor(min(c(sqrt(n), p)))
  if(!is.null(sigma.hat) || method == 'root-lasso') {
    sigma_hat <- sigma.hat
  } else if(round(C) != C || C <= 0) {
    sigma_hat <- var(y)
  } else {
    sigma_hat <- min(var(y), summary(lm(y ~ x[,sample(1:p,C)]))$sigma)
  }
  lambda <- sigma_hat * sqrt((log(p) + 1) / n)

  if(method == 'lasso') {
    hdlm.stat <- function(x, INDEX) {
      # Estimate model:
      out <- lars::lars(x[INDEX,],y[INDEX])
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
  }
  else if(method == 'mc+') {
    lambda <- lambda/sqrt(n)
    hdlm.stat <- function(x, INDEX) {
      # Estimate model:
      out <- plus::plus(x[INDEX,],y[INDEX], method='mc+', intercept=FALSE)
      res <- NULL
      if(lambda > out$lam.path[1]) {
        res <- rep(0,p)
      }
      if(lambda < out$lam.path[length(out$lam.path)]) {
        res <- out$beta[nrow(out$beta),]
      }
      if(is.null(res)) {
        i <- max(which(out$lam.path > lambda))
        betaHat <- out$beta[i:(i+1),]
        lam <- out$lam.path[(i+1):i]
        alpha <- 1-(lambda - lam[1])/(lam[2] - lam[1])
        res <- betaHat[1,]*(1-alpha) + betaHat[2,]*alpha
      }
      return(res)
    }
  }
  else if(method == 'scad') {
    lambda <- lambda/sqrt(n)
    hdlm.stat <- function(x, INDEX) {
      # Estimate model:
      out <- plus::plus(x[INDEX,],y[INDEX], method='scad', intercept=FALSE)
      res <- NULL
      if(lambda > out$lam.path[1]) {
        res <- rep(0,p)
      }
      if(lambda < out$lam.path[length(out$lam.path)]) {
        res <- out$beta[nrow(out$beta),]
      }
      if(is.null(res)) {
        i <- max(which(out$lam.path > lambda))
        betaHat <- out$beta[i:(i+1),]
        lam <- out$lam.path[(i+1):i]
        alpha <- 1-(lambda - lam[1])/(lam[2] - lam[1])
        res <- betaHat[1,]*(1-alpha) + betaHat[2,]*alpha
      }
      return(res)
    }
  }
  else if(method == 'root-lasso') {
    c <- 1.1
    alpha <- 0.05
    lambda <- c * sqrt(n) * qnorm(1-alpha/(2*p))
    MaxIter <- 10000
    OptTolNorm <- 1e-6
    OptTolObj <- 1e-8

    hdlm.stat <- function(x, INDEX) {
      # Estimate model:
      out <- rootlasso(x[INDEX,],y[INDEX],lambda,MaxIter,OptTolNorm,OptTolObj)
      return(out)
    }
  }  else if(method == 'dantzig') {

    hdlm.stat <- function(x, INDEX) {
      # Estimate model:
      out <- dselector(x[INDEX,],y[INDEX],lambda)
      return(out)
    }
  }

  # Does user want bootstrap standard errors and bias calculations?
  bootSE <- !is.null(N)
  if(bootSE == FALSE & p.value.method == 'bootstrap') N <- 100
  if(bootSE == FALSE & is.null(N)) N <- 1 	
  if(abs(N - round(N)) > .Machine$double.eps^0.5 || N < 2) {
    N <- 1
    bootSE <- FALSE
  }

  # Bootstrap to determine point estimates and standard deviation
  B  <- boot(x, hdlm.stat, N)

  # Standard estimators:
  point_estimator <- B$t0
  fitted <- x %*% point_estimator
  resid <- fitted - y
  if(method == 'root-lasso') sigma_hat <- mean(resid^2)

  # Calculate p-values and standard errors
  hdlm.pval.boot <- function(v) {
    1-mean(v[-1] == v[1])*abs(v[1])
  }

  if(bootSE == TRUE) {
    se <- apply(B$t,2,sd)
    bias <- apply(B$t,2,mean) - B$t0
  } else {
    se <- NULL
    bias <- NULL
  }

  pvalue <- NULL
  if(p.value.method == 'bootstrap') {
    pvalue <- apply(sign(rbind(B$t0, B$t)), 2, hdlm.pval.boot)
  } else if(p.value.method == 'one-split'){
    I <- sample(1:n, floor(n/2))
    point_estimator2 <- hdlm.stat(x, I) 
    index <- point_estimator2 != 0
    if(sum(index) < min(n,p)) {
      pvalue <- rep(1,p)
      if(sum(index) != 0) pvalue[index] <- as.numeric(summary(lm( y[-I] ~ x[-I,index] - 1))[[4]][,4])
    } else {
      warning('cannot computer p-values by one-split method; point estimate not sufficently sparse\n consider setting p.value.method to bootstrap')
    }
  }


  z <- list(coefficients=point_estimator, residuals=resid, effects=NULL, rank=c(n,p),
            fitted.values=fitted, assign=NULL, standard.error=se, bias=bias, p.value=pvalue,
            sigma.hat = sigma_hat, N=N)

  return(z)

}

