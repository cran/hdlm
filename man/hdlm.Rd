% Fits high dimensional models, code adapted from
% File src/library/stats/man/lm.Rd, from R Core Development Team
\name{hdlm}
\alias{hdlm}
\alias{hdlm.fit}
\alias{HDprintCoefmat}
\alias{print.hdlm}
\alias{print.summary.hdlm}
\alias{rootlasso}
\alias{dselector}
\alias{summary.hdlm}
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
\title{
Fitting High Dimensional Linear Models
%%  ~~function to do ... ~~
}
\description{
\code{hdlm} is used to fit high dimensional linear models when
  the number of predictors is greater than (or close to the same
  order of magnitude as) the number of observations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hdlm(formula, data, subset,
      method=c('mc+',  'root-lasso',  'scad', 'lasso','dantzig'),
      p.value.method = c('one-split', 'bootstrap', 'none'),
      model = TRUE, x = FALSE, y = FALSE, N = NULL, C = NULL,
      sigma.hat = NULL, ...)
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

  \item{method}{method for performing the optimization. See below for details.}

  \item{p.value.method}{method for performing the calculating p-value. The default
    splits the data into two equal sets, uses 'method' for model selection on one
    half and the other half conducts ordinary least squares if reduced model is
    small enough (<= p) for it to run. Otherwise, users can choose to bootstrap
    p-values or not calculate them.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.
  }

  \item{N}{
    Number of bootstrap replicates to run in order to estimate standard error
    and bias. If set to NULL (default) bootstrap estimates will not be calculated.
    Notice that if p.value.method is selected as bootstrap and N is NULL, no
    bootstrap standard errors will be produced.
  }
  \item{C}{
    Size of subset used to determine an upper bound on the noise level.
    If set to NULL, a non-positive integer, or an integer greater than n,
    C is assumed to be zero. Not used for root-lasso, which does not need
    to estimate the noise.
  }
  \item{sigma.hat}{
    Allows user to supply estimate of the noise variance. Overrides option 
    for C above. Is ignored when method='root-lasso', as the noise level
    is not needed for the square root lasso method.
  }

  \item{\dots}{additional arguments to be passed to the low level
    regression fitting functions (see below).}

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
  be penalized along with other terms; for this reason it is recommended
  that an intercept be excluded from the model.

  Each of five current methods for running high dimensional regression have
  references below. The options 'mc+' and 'scad' call lower level routines
  from the package plus. Similarly, 'lasso' calls a function from the package
  lars. The code for the 'root-lasso' method was adapted from matlab code
  created by Alex Belloni. For 'dantzig' the quantreg package is utilized
  as suggested by Roger Koenker.

}
\value{
  \code{hdlm} returns an object of \code{\link{class}} \code{"hdlm"}.

  The function \code{summary} are used to
  obtain and print a summary and analysis of variance table of the
  results.  The generic accessor functions \code{coefficients},
  \code{effects}, \code{fitted.values} and \code{residuals} extract
  various useful features of the value returned by \code{hdlm}.
}
\references{
  Belloni, A., V. Chernozhukov, and L. Wang. (2011) "Square-Root Lasso: Pivotal
  Recovery of Sparse Signals Via Conic Programming". In: Arxiv preprint
  arXiv:1009.5689.

  Bickel, P.J., Y. Ritov, and A.B. Tsybakov (2009) "Simultaneous analysis of Lasso
  and Dantzig selector". \emph{The Annals of Statistics} 37.4, pp. 1705--1732.

  Buhlmann, P. and S. Van De Geer (2011) \emph{Statistics for High-Dimensional Data:
  Methods, Theory and Applications}. Springer-Verlag New York Inc.

  Candes, E. and T. Tao (2007) "The Dantzig selector: Statistical estimation when p is much
  larger than n". \emph{The Annals of Statistics} 35.6, pp. 2313--2351.

  Chambers, J. M. (1992)
  \emph{Linear models.}
  Chapter 4 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
   (with discussion) \emph{Annals of Statistics}; see also \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.

  Fan, J., Y. Feng, and Y. Wu (2009) "Network exploration via the adaptive LASSO and SCAD penalties".
  \emph{Annals of Applied Statistics} 3.2, pp. 521--541.

  Hastie, Tibshirani and Friedman (2002) Elements of Statistical
  Learning, Springer, NY.

  Wasserman, L., and Roeder, K. (2009), "High Dimensional Variable Selection,"
  The Annals of Statistics, 37, 2178--2201.

  Zhang, C.H. (2010) "Nearly unbiased variable selection under minimax concave
  penalty". \emph{The Annals of Statistics} 38.2, pp. 894--942
}
\author{
  Created by Taylor B. Arnold for point estimation and confidence
  intervals in high-dimensional regression. 
  The design of the function was inspired by the S/R function
  \code{lm} described in Chambers (1992). 
}
\note{
  The reported bootstrap bias and standard errors give the estimated
  bias and errors for the estimator, NOT for the underlying parameters. A distinction
  is necessary  since all of the estimators are biased. The method 'mc+' is chosen as
  a default method because it is 'nearly un-biased', and hence the distinction is negligible
  in most situations. If the user is only interested in point estimates and not standard
  errors, the 'root-lasso' method will generally be preferred as it does not require
  explicit estimation of the noise variance. The Dantzig selector benefits from
  the fastest execution time amongst our set of methods, and theory suggests that it is
  consistent even in the case of very high colinearity in the model matrix.

  This package focuses on methods which both (i) produce sparse estimates and (ii) result
  from a global optimization problem. Users who do not require sparse estimates are directed
  to other methods such as ridge regression, elastic net, and the Bayesian lasso. For
  algorithms which do not fall under point (ii), popular choices are forward stagewise
  regression and OMP. 
}
\section{Future Additions}{
  In the future, we hope to focus on optimizing the code for large datasets (an obviously
  common situation in high-dimensional linear models). This will include using parallel
  computing for the bootstrap method, as well as being more strategic by calculating
  various quantities (the gram matrix X'X, for instance) no more than necessary. 

  Additionally, we hope to include the p-value methodology by Meinshaussen, Meier, and
  Buhlmann, which while more computationally demanding provides a better estimate than
  the single split method. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

  \code{\link{lm}} for fitting traditional low dimensional models.

  \code{\link{lars}}, \code{\link{dselector}}, \code{\link{plus}}, and \code{\link{rootlasso}}
   for lower level functions for the various methods.
}

\examples{
set.seed(1)
x <- matrix(rnorm(40*30), ncol=40, nrow=30)
y <- x[,1] + x[,2] + rnorm(30,sd=0.5)

out <- hdlm(y ~ x - 1, N = 100, method='lasso')
summary(out)

out <- hdlm(y ~ x - 1, method='root-lasso')
summary(out)

out <- hdlm(y ~ x - 1, N = 100, method='mc+')
summary(out)

}
\keyword{regression}
\keyword{high dimensional}
