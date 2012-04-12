% Fits high dimensional models, code adapted from
% File src/library/stats/man/lm.Rd, from R Core Development Team
\name{hdlm}
\alias{hdlm}
\alias{hdlm.fit}
\alias{HDprintCoefmat}
\alias{print.hdlm}
\alias{print.summary.hdlm}
\alias{plot.hdlm}
\alias{residuals.hdlm}
\alias{qr.hdlm}
\alias{simulate.hdlm}
\alias{deviance.hdlm}
\alias{formula.hdlm}
\alias{family.hdlm}
\alias{model.frame.hdlm}
\alias{variable.names.hdlm}
\alias{case.names.hdlm}
\alias{anova.hdlm}
\alias{predict.hdlm}
\alias{effects.hdlm}
\alias{model.matrix.hdlm}
\alias{plot.hdlm}
\alias{labels.hdlm}
\alias{df.residual.hdlm}
\alias{bayes.hdlm.fit}
\alias{bayesianLatentFit}
\alias{mod.cv.glmnet}
\alias{votingRecord}
\alias{wordDataset}
\title{
Fitting High Dimensional Linear Models
}
\description{
  \code{hdlm} is used to fit high dimensional linear models when
  the model matrix is rank deficent. The default usage is the same
  as the lm function in stats; for instance running the code:
  'summary(hdlm(y ~ x))' will produce a regression table. A myriad
  of options are also avaliable, as described below. For technical
  and theoretical details of the underlying methods see the Details
  section below as well.
}
\usage{
hdlm(formula, data, subset, bootstrap = 10, siglevel = 0.05,
    alpha = 0.5, M = NULL, N = NULL, model = TRUE, x = FALSE,
    y = FALSE, scale=TRUE, pval.method=c('median', 'fdr', 'holm', 'QA'),
    ..., FUNCVFIT = NULL, FUNLM = NULL, bayes=FALSE, bayesIters=NULL,
    bayesTune = c(1,1), refit=FALSE) 
 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{lm} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{bootstrap}{
    number of bootstrap trails to conduct. Default is 10.
  }

  \item{siglevel}{
    significance level to use for confidence bounds. Default is 0.05.
  }

  \item{alpha}{
    elastic net mixing parameter sent to glmnet, can be any value in (0,1].
    When alpha = 1, this is the lasso penalty and when alpha = 0 (not supported)
    this is the ridge penalty. See glmnet help pages for more details.
  }

  \item{M}{
    maximum model size sent to the second stage low-dimensional regression. 
    When more than M variables are choosen in the first stage, the model is
    trimmed by succesively taking larger sized coefficents until only M
    remain. If NULL, M is taken to be 90% of the number of observations used
    in the second stage. If M = 0, the model is fit with all of the data once,
    and the estimated parameters are returned as is.
  }

  \item{N}{
    Numer of observations to include in the first stage regression. Default
    is (# samples / 2), so that the data is split evenly amongst the two 
    stages, which will be set when N=NULL.
  }

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.
  }

  \item{scale}{
   Logical; should the variables in the data matrix be scaled.
  }

  \item{pval.method}{
    one of 'median', 'fdr', 'holm', or 'QA'. Signifies the method used to 
    combine p-values when bootstrap is greater than 1. For details
    and relative strengths of the three methods, see the package
    vignette. 
  }

  \item{\dots}{additional arguments to be passed to the low level
    regression fitting functions (see below).}

  \item{FUNCVFIT}{
    Used to pass alternative model selection function. Must accept data matrix
    as its first element and response vector as second element. Return should
    be a vector of length p (the number of regressors), which indicates which
    variables are included in the final model. Zero terms are considered to be
    out of the model; typically all non-zero terms are treated as in the model,
    though if the model size is too large (see 'M' above), it will be trimmed
    relative to the absolute size of each non-zero term. Therefore, it is advised
    to return the model vector in a relative scale rather than an absolute one.
    The default, used when NULL, is the elastic net function from package glmnet,
    with the mixing parameter alpha from above. See package vignette for additional
    details and examples.
  }

  \item{FUNLM}{
    Used to pass alternative second stage, low-dimensional function. Must accept
    as its first argument a formula object. The return class must have a summary
    method and the summary method in turn must have a coef method. The coef.summary
    should return a matrix where the first column are the coefficients and the second
    column are standard errors. Intercepts should be handled according to the passed
    formula. As an example, stats::lm works by default; stats::lm is additionally the 
    default when FUNLM is set to NULL. See package vignette for additional details
    and examples. 
  }

  \item{bayes}{
    logical. Should Bayesian method be used in place of the two stage method. 
  }

  \item{bayesIters}{
    number of iterations to conduct in the Gibbs sampler when bayes=TRUE. A
    total of (bayesIters * 0.1) burn-in steps are included as well. Default
    is 1000, and can be set by setting bayesIters = NULL.
  }

  \item{bayesTune}{
    numerical vector tuning parameter for the Bayes estimator. Defines a 
    Beta(bayesTune[1], bayesTune[2]) prior on the proportion of variables
    included in the true support. 
  }

  \item{refit}{
    Either a logical or number in (0,1]. When not equal to false, the final
    model will be refit from the entire dataset using FUNLM. When a numeric,
    the model is selected by only including variables with p-values less than
    refit. When set to TRUE, any variable corrisponding to a non-zero p-value
    is included. 
  }
}
\details{
  Models for \code{hdlm} are specified symbolically.  A typical model has
  the form \code{response ~ terms} where \code{response} is the (numeric)
  response vector and \code{terms} is a series of terms which specifies a
  linear predictor for \code{response}.  A terms specification of the form
  \code{first + second} indicates all the terms in \code{first} together
  with all the terms in \code{second} with duplicates removed.  A
  specification of the form \code{first:second} indicates the set of
  terms obtained by taking the interactions of all terms in \code{first}
  with all terms in \code{second}.  The specification \code{first*second}
  indicates the \emph{cross} of \code{first} and \code{second}.  This is
  the same as \code{first + second + first:second}.

  If the formula includes an \code{\link{offset}}, this is evaluated and
  subtracted from the response.

  See \code{\link{model.matrix}} for some further details.  The terms in
  the formula will be re-ordered so that main effects come first,
  followed by the interactions, all second-order, all third-order and so
  on: to avoid this pass a \code{terms} object as the formula (see
  \code{\link{aov}} and \code{demo(glm.vr)} for an example).

  A formula has an implied intercept term.  To remove this use either
  \code{y ~ x - 1} or \code{y ~ 0 + x}.  See \code{\link{formula}} for
  more details of allowed formulae. Note that the intercept term will
  not be penalized along with other terms. If you want a penalized
  intercept, add it to directly to the matrix x.

}
\value{
  \code{hdlm} generally returns an object of \code{\link{class}} \code{"hdlm"},
  unless refit is not set to false. In the latter case the output is dependent
  on the choice of funtion FUNLM.

  The function \code{summary} is used to obtain and print a summary of the
  results. The generic accessor functions \code{coefficients},
  \code{effects}, \code{fitted.values} and \code{residuals} extract
  various useful features of the value returned by \code{hdlm}.
}
\references{

  Bickel, P.J., Y. Ritov, and A.B. Tsybakov (2009) "Simultaneous analysis of Lasso
  and Dantzig selector". \emph{The Annals of Statistics} 37.4, pp. 1705--1732.

  Buhlmann, P. and S. Van De Geer (2011) \emph{Statistics for High-Dimensional Data:
  Methods, Theory and Applications}. Springer-Verlag New York Inc.

  Chambers, J. M. (1992) \emph{Linear models.}
  Chapter 4 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
   (with discussion) \emph{Annals of Statistics}; see also \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.

  Fan, J., Y. Feng, and Y. Wu (2009) "Network exploration via the adaptive LASSO and SCAD penalties".
  \emph{Annals of Applied Statistics} 3.2, pp. 521--541.

  Hastie, Tibshirani and Friedman (2002) Elements of Statistical
  Learning, Springer, NY.

  Hans, C. (2009). Brief Technical Report to Accompany the R Package blasso Bayesian Lasso
  Regression. URL http://www.stat.osu.edu/~hans/software/blasso/.

  Wasserman, L., and Roeder, K. (2009), "High Dimensional Variable Selection,"
  The Annals of Statistics, 37, 2178--2201.

}
\author{
  Created by Taylor B. Arnold for point estimation and confidence
  intervals in high-dimensional regression. 

  The Bayesian option for package hdlm is as implemented with Gibbs
  sampling with C code from Chris hans, avaliable as packged with
  package 'blasso' from: \url{www.stat.osu.edu/~hans/software/blasso/}

  The design of the function was inspired by the S/R function
  \code{lm} described in Chambers (1992). 
}
\examples{
  set.seed(1)
  x <- matrix(rnorm(100*40),ncol=100)
  y <- x[,1] + x[,2] * 0.5 + rnorm(40, sd=0.1)
  out <- hdlm(y ~ x)
  summary(out)
}
\note{
  This package focuses on methods which produce sparse estimates. Users who do not
  require sparse estimates are directed to other methods such as ridge regression.
}
