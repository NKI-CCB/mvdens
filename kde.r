library(mvtnorm)

fit.kde <- function(x, adjust = 1, bw.fn = bw.SJ, diagonal = T, verbose = F, ...)
{
    n <- nrow(x)
    d <- ncol(x)

    result <- list()
    result$type <- "kde"
    result$dim <- ncol(x)
    result$x <- x

    if (diagonal) {
        factor <- n ^ (-1 / (d + 4))
        bw <- adjust * factor * apply(x, 2, bw.fn, ...)
        result$H <- diag(bw)
    } else {
        result$H <- bw.fn(x, ...)
    }
    return(structure(result, class = "mdd.density"))
}

fit.kde.transformed <- function(x, bounds, adjust = 1, bw.fn = bw.SJ, diagonal = T, verbose = F)
{
    result <- list()
    result$type <- "kde.transformed"
    result$transform.bounds <- bounds
    transformed <- mdd.transform_to_unbounded(x, bounds)
    result$kde <- fit.kde(transformed, adjust = adjust, bw.fn = bw.fn, diagonal = diagonal)
    return(structure(result, class = "mdd.density"))
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
