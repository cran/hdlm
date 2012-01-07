\name{summary.hdlm}
\alias{summary.hdlm}
\alias{summary.hdglm}
\title{summarize coefficients from a "glmnet" object}
\description{
    'summary' method for class '"lm"'.
}
\usage{
\method{summary}{hdlm}(object, ..., level=NULL)
}
\arguments{
  \item{object}{an object of class '"hdlm"', usually, a result of a call to
          'hdlm'.}
  \item{level}{Determines which coefficients to print on regression table.
               Level = 1 gives only those with non-zero coefficents if
               pval.method is equal to 'mean' and only  those with p-value < 1
               otherwise. Level = 2 gives anything with non-zero coefficient
               or non-one p-value, and Level = 3 (or any other choice) gives
               all coefficients.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  When fitting a large model, it can be cumbersome to look at results for all of
  the variables; when there exists high correlation between variables, level = 2
  is often preferable. The different behavior for pval.method equal to 'mean' is
  due to the fact that the mean method (or Bayes method, which sets pval.method to
  'mean') gives many p-values which are close to, but not exactly equal to, one.
}
% R documentation directory.
\keyword{regression}
