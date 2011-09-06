% Fits high dimensional models, code adapted from
% File src/library/stats/man/lm.Rd, from R Core Development Team

\name{hdlm}
\alias{hdlm}
\alias{hdlm.fit}
\alias{HDprintCoefmat}
\alias{print.hdlm}
\alias{print.summary.hdlm}
\alias{summary.hdlm}
\title{
Fitting High Dimensional Linear Models
%%  ~~function to do ... ~~
}
\description{
\code{hdlm} is used to fit high dimensional linear models.
  It is used to carry out regression when the number of 
  predictors is greater than (or close to the same order of
  magnitude as) the number of observations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hdlm(formula, data, subset, model = TRUE,
      x = FALSE, y = FALSE, N = 100, C = NULL, ...)
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

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.
  }

  \item{N}{
  Number of bootstrap replicates to run in order to estimate standard errors
  and p-values.
}
  \item{C}{
  Size of subset used to determine an upper bound on the noise level.
  If set to NULL, a non-positive integer, or an integer greater than n,
  C is assumed to be zero.
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

}
\value{
  \code{hdlm} returns an object of \code{\link{class}} \code{"hdlm"}.

  The function \code{summary} are used to
  obtain and print a summary and analysis of variance table of the
  results.  The generic accessor functions \code{coefficients},
  \code{effects}, \code{fitted.values} and \code{residuals} extract
  various useful features of the value returned by \code{lm}.
}
\references{
  Chambers, J. M. (1992)
  \emph{Linear models.}
  Chapter 4 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression"
(with discussion) \emph{Annals of Statistics}; see also \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}.

Hastie, Tibshirani and Friedman (2002) Elements of Statistical
Learning, Springer, NY.
}
\author{
  Created by Taylor B. Arnold for point estimation and confidence
  intervals in high-dimensional regression. Uses lasso algorithm
  described in Efron, Hastie, Johnstone and Tibshirani (2003).
  The design of the function was inspired by the S/R function
  \code{lm} described in Chambers (1992). 
}
\note{
  For more details about internal procedures, see internal comments in source code.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lm}} for fitting traditional low dimensional models.

  \code{\link{lars}} for underlying lasso algorithm.
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
set.seed(1)
x <- matrix(rnorm(40*30), ncol=40, nrow=30)
y <- x[,1] + x[,2] + rnorm(30,sd=0.5)
out <- hdlm(y ~ x - 1)
summary(out)

}
\keyword{regression}
\keyword{high dimensional}
