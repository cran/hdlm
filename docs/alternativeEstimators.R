# As discussed in the package vignette, we have here code
# for alternative point estimators. These are not included
# in the actual package so as to decrease the number of 
# dependent packages (which grow quite large) and also to 
# make it easier to adapt these to specific usages.

# The current set is only a start and more estimators will
# be added in the coming months.

# Full documentation of each method can be found in the
# relavent packages.

# - Taylor Arnold


###############################
#####  FUNLM alternatives  ####
###############################

# (1) Quantile regression
library("quantreg")
FUNLM <- function(formula) {
  out <- rq(formula)
  class(out) <- "rq_new"
  return(out)
}
summary.rq_new <- function(out) {
  class(out) <- 'rq'
  val <- summary.rq(out, se='nid')
  return(val)
}

# (2) Robust MM-Estimation
library(MASS)
FUNLM <- function(formula) return(rlm(formula, maxit = 50))


###############################
#### FUNCVFIT alternatives ####
###############################

# (1) Dantzig selector
library(quantreg)
FUNCVFIT <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  lambda <- sqrt(log(ncol(x)+1) * nrow(x)) * sd(as.numeric(y))

  A <- t(x) %*% x
  R <- rbind(A, -A)
  a <- c(as.matrix(t(x) %*% y))
  r <- c(a-lambda, -a-lambda)
  beta <- quantreg::rq.fit.fnc(diag(p), rep(0,p), R=R, r=r)$coefficients
  return(round(beta, 6))
}

# (2) MC+ algorithm
library(plus)
FUNCVFIT <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  lambda <- sqrt(log(ncol(x)+1) * nrow(x)) * sd(as.numeric(y))

  out <- plus(x,y, method='mc+')
  if(lambda > max(out$lam.path) | lambda < min(out$lam.path)) lambda <- mean(out$lam.path) 
  beta <- coef(out, lam = lambda)
  return(beta)
}


# (3) SCAD
library(plus)
FUNCVFIT <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  lambda <- sqrt(log(ncol(x)+1) * nrow(x)) * sd(as.numeric(y))

  out <- plus(x,y, method='scad')
  if(lambda > max(out$lam.path) | lambda < min(out$lam.path)) lambda <- mean(out$lam.path) 
  beta <- coef(out, lam = lambda)
  return(beta)
}

# (4) Marginal regression (NEED TO SET PARAMTER 'M' in the model)
FUNCVFIT <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)

  x <- x - matrix(apply(x,2,mean),ncol=ncol(x),nrow=nrow(x),byrow=TRUE)
  x <- x / sqrt(matrix(apply(x^2,2,mean),ncol=ncol(x),nrow=nrow(x),byrow=TRUE))
  y <- y - mean(y)
  y <- y / sqrt(mean(y^2))

  out <- abs(y %*% x)
  out <- (out - min(out)) / (max(out) - min(out))
  return(out)
}








