hdlm <-
function (formula, data, subset,
    model = TRUE, x = FALSE, y = FALSE,  N=100, C = NULL, ...) 
{
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
    z <- hdlm.fit(x, y, N, C, ...)

    class(z) <- c("hdlm")
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

