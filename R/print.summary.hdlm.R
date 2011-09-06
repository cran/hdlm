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
    MAT <- cbind(x$coefficients, x$bias, x$standard.error, x$p.value)
    colnames(MAT) <- c("Estimate", "Bootstrap Bias", "Std. Error", "Pr(>|t|)")
    HDprintCoefmat(MAT, digits=digits, signif.stars=signif.stars, na.print="NA")
    
    
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

