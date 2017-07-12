library(mvtnorm)

fit.kde <- function(x, adjust = 1, bw.fn = bw.SJ)
{
    n <- nrow(x)
    p <- ncol(x)
    factor <- (4 / (n * (p + 2))) ^ (1 / (p + 4))
    bw <- adjust * factor * apply(x, 2, bw.SJ)

    result <- list()
    result$type <- "kde"
    result$x <- x
    result$H <- diag(bw)

    return(result)
}

fit.kde.transformed <- function(x, bounds, adjust = 1, bw.fn = bw.SJ)
{
    result <- list()
    result$type <- "kde.transformed"
    result$transform.bounds <- bounds
    transformed <- mdd.transform_to_unbounded(x, bounds)
    result$kde <- fit.kde(transformed, adjust = adjust, bw.fn = bw.fn)
    return(result)
}

evaluate.kde <- function(fit, x, log = FALSE)
{
    tp <- rep(NA, nrow(x))
    for (i in 1:nrow(x)) {
        p <- dmvnorm(fit$x, x[i,], fit$H)
        if (log) {
            tp[i] <- log(sum(p)) - log(nrow(fit$x))
        } else {
            tp[i] <- sum(p) / nrow(fit$x)
        }
    }
    return(tp)
}
