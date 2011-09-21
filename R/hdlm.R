# Much of this is taken from the methods for lm found in the recommended
# package stats. 

# We note differences between this file and lm.R found in the directory
# /src/library/stats/R of the standard R source files, expect where
# we have simply left out wlm methods (which we do not yet have analogues for).
hdlm <-
function (formula, data, subset, method=c('mc+', 'root-lasso', 'scad', 'lasso', 'dantzig'),
    p.value.method = c('one-split', 'bootstrap', 'none'),
    model = TRUE, x = FALSE, y = FALSE,  N=NULL, C = NULL, sigma.hat = NULL, ...) 
{

    method <- method[[1]]
    p.value.method <- p.value.method[[1]]
    if(!(method %in% c('mc+', 'scad', 'lasso', 'root-lasso', 'dantzig'))) stop('unsupported method entered')
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))

    x <- model.matrix(mt, mf, contrasts)
    if(!is.matrix(x)) x <- as.matrix(x)
    z <- hdlm.fit(x, y, method, p.value.method, N, C, sigma.hat, ...)

    class(z) <- c("hdlm")
    names(z$coefficients) <- colnames(x)
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    z
}

# fit.hdlm moved to own file, due to space and need to edit for new options easily

print.hdlm <-
function(x, ...){
    digits = max(3, getOption("digits") - 3)
    z <- x
    cat("\nCall:\n", paste(deparse(z$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Parameters:\n", " Observations = ", z$rank[1], ", Variables = ", z$rank[2], 
         ", Bootstrap Replicates = ", z$N, sep="")    
    cat("\n\n")
    invisible(z)
}

summary.hdlm <-
function(object, ...) {
  z <- object
  class(z) <- 'summary.hdlm'
  z
}


print.summary.hdlm <-
function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
      signif.stars= getOption("show.signif.stars"),...)
{
    cat("\nCall:\n", # S has ' ' instead of '\n'
    paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    resid <- x$residuals
    df <- x$rank
    rdf <- x$rank[1]
    if (rdf > 5L) {
        cat("Residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    zz <- zapsmall(quantile(resid), digits + 1)
    rq <- structure(zz, names = nam)
        
    print(rq, digits = digits, ...)
    } else if (rdf > 0L) {
        print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
        cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
    }

    cat("\nCoefficients:\n")
    if(is.null(x$p.value) & is.null(x$standard.error)) {
        MAT <- cbind(x$coefficients)
        colnames(MAT) <- c("Estimate")
        HDprintCoefmat(MAT, digits=digits, signif.stars=FALSE, na.print="NA", pval=FALSE)
    } else if(is.null(x$p.value) & !is.null(x$standard.error)) {
        MAT <- cbind(x$coefficients, x$bias, x$standard.error)
        colnames(MAT) <- c("Estimate", "Bootstrap Bias", "Std. Error")
        HDprintCoefmat(MAT, digits=digits, signif.stars=FALSE, na.print="NA", pval=FALSE)
    } else if(!is.null(x$p.value) & is.null(x$standard.error)) {
        MAT <- cbind(x$coefficients, x$bias, x$standard.error, x$p.value)
        colnames(MAT) <- c("Estimate", "Pr(>|t|)")
        HDprintCoefmat(MAT, digits=digits, signif.stars=signif.stars, na.print="NA")
    } else if(!is.null(x$p.value) & !is.null(x$standard.error)) {
        MAT <- cbind(x$coefficients, x$bias, x$standard.error, x$p.value)
        colnames(MAT) <- c("Estimate", "Bootstrap Bias", "Std. Error", "Pr(>|t|)")
        HDprintCoefmat(MAT, digits=digits, signif.stars=signif.stars, na.print="NA")
    }

 
    cat("\nEstimated sigma:",
    format(signif(x$sigma.hat, digits)), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-squared:", formatC(x$r.squared, digits=digits))
        cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
        "\nF-statistic:", formatC(x$fstatistic[1L], digits=digits),
        "on", x$fstatistic[2L], "and",
        x$fstatistic[3L], "DF,  p-value:",
        format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE), digits=digits),
        "\n")
    }
    cat("\n")#- not in S
    invisible(x)
}


residuals.hdlm <-
    function(object,
             type = c("working","response", "deviance","pearson", "partial"),
             ...)
{
    r <- object$residuals
    res <- r
    res
}

qr.hdlm <- function(x, ...) {
    stop("hdlm objects do not have a qr decomposition; if n <= p, one
          can be generated by function lm.")
}


simulate.hdlm <- function(object, nsim = 1, seed = NULL, ...)
{
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
    if(is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
	set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ftd <- fitted(object)             # == napredict(*, object$fitted)
    nm <- names(ftd)
    n <- length(ftd)
    ntot <- n * nsim
    fam <- "gaussian"
    vars <- deviance(object)/df.residual(object)
    
    val <- ftd + rnorm(ntot, sd = sqrt(vars))

    if(!is.list(val)) {
        dim(val) <- c(n, nsim)
        val <- as.data.frame(val)
    } else class(val) <- "data.frame"
    
    names(val) <- paste("sim", seq_len(nsim), sep="_")
    if (!is.null(nm)) row.names(val) <- nm
    attr(val, "seed") <- RNGstate
    val
}


deviance.hdlm <- function(object, ...) {
    sum(residuals(object)^2, na.rm=TRUE)
}

# Extra df.residual method here for hdlm.
df.residual.hdlm <- function(object, ...) {
    return(object$rank[[1]])
}

formula.hdlm <- function(x, ...)
{
    form <- x$formula
    if( !is.null(form) ) {
        form <- formula(x$terms) # has . expanded
        environment(form) <- environment(x$formula)
        form
    } else formula(x$terms)
}

family.hdlm <- function(object, ...) { gaussian() }


model.frame.hdlm <- function(formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1L]] <- as.name("lm")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
	if (is.null(env)) env <- parent.frame()
        eval(fcall, env, parent.frame())
    }
    else formula$model
}


variable.names.hdlm <- function(object, full = FALSE, ...)
{
    if(object$rank[[1]]) colnames(model.matrix(object))
    else character()
}

case.names.hdlm <- function(object, full = FALSE, ...)
{
    dn <- as.character(1:(object$rank[[1]]))
    dn 
}


anova.hdlm <- function(object, ...)
{
  stop('anova not yet avaliable for hdlm objects')
}

# version of anova.lmlist not yet included

## code originally from John Maindonald 26Jul2000
predict.hdlm <-
    function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
	     interval = c("none", "confidence", "prediction"),
	     level = .95,  type = c("response", "terms"),
	     terms = NULL, na.action = na.pass, weights = 1, ...)
{
    tt <- object
    # Simple / less complete than lm version:
    if(!inherits(object, "hdlm"))
	warning("calling predict.hdlm(<fake-lm-object>) ...")
    if(missing(newdata) || is.null(newdata)) {
	X <- model.matrix(object)
    } else {
        X <- model.matrix(hdlm(object$fitted.values ~ newdata, N=1))
       if(variable.names(object)[[1]] != "(Intercept)") X <- X[,-1]
    }

    beta <- coef(object)
    ynew <- X %*% beta

    return(ynew)
}


effects.hdlm <- function(object, set.sign = FALSE, ...)
{
    stop("'object' has no 'effects' component")
}

## plot.hdlm is written here, unlike with plot.lm
plot.hdlm <- function (x, which = c(1L:3L, 5L), caption = list("Residuals vs Fitted", 
    "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
    expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), 
    panel = if (add.smooth) panel.smooth else points, sub.caption = NULL, 
    main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, 
    qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2), cex.caption = 1) 
{
    if (!inherits(x, "hdlm")) 
        stop("use only with \"hdlm\" objects")
    xlab <- paste('fitted values')
    if(sum(x$rank) < 1e4) {
      plot(fitted(x), residuals(x), pch=19, xlab=xlab,
             ylab='residuals', main='Residuals vs Fitted')
    } else {
      plot(fitted(x), residuals(x), pch='.', xlab=xlab,
             ylab='residuals', main='Residuals vs Fitted')
    }
}

model.matrix.hdlm <- function(object, ...)
{
    if(n_match <- match("x", names(object), 0L)) object[[n_match]]
    else {
        data <- model.frame(object, xlev = object$xlevels, ...)
        NextMethod("model.matrix", data = data,
                   contrasts.arg = object$contrasts)
    }
}

## from base/R/labels.R
labels.hdlm <- function(object, ...)
{
    formula(object)
}







