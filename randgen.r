library(mvtnorm)

mdd.rgen <- function(fit, n) {
    stopifnot(class(fit) == "mdd")

    if (fit$type == "kde") {
        s <- sample(nrow(fit$x), n, replace = T)

        if (all(m[!diag(nrow(fit$H))] == 0)) {
            for (i in 1:n) {
                x[i,] <- rnorm(ncol(fit$x), fit$x[s[i],], diag(fit$H))
            }
        } else {
            x <- rmvnorm(n, rep(0, ncol(fit$x)), fit$H)
            x <- x + fit$x[s,]
        }
    } else if (fit$type == "kde.transformed") {
        transformed <- mdd.rgen(fit$kde, n)
        x <- mdd.transform_from_unbounded(transformed, fit$transform.bounds)
    } else if (fit$type == "gmm") {
        s <- sample(fit$k, n, replace = T)
        for (i in 1:n) {
            x[i,] <- rmvnorm(1, fit$centers[s[i],], fit$covariances[[s[i]]])
        }
    } else if (fit$type == "gmm.transformed") {
        transformed <- mdd.rgen(fit$gmm, n)
        x <- mdd.transform_from_unbounded(transformed, fit$transform.bounds)
    } else if (fit$type == "gmm.truncated") {
        if (ncol(x) < 10) {
            method <- "rejection"
        } else {
            method <- "gibbs"
        }
        s <- sample(fit$k, n, replace = T)
        for (i in 1:n) {
            x[i,] <- rtmvnorm(1, fit$centers[s[i],], fit$covariances[[s[i]]], fit$bounds[,1], fit$bounds[,2], algorithm=method)
        }
    } else if (fit$type == "gp") {
        stop("Cannot sample from density function described by gaussian process")
    } else if (fit$type == "vine.copula") {
        sim <- RVineSim(n, fit$RVM)
        x <- reverse.transform.marginals(sim, fit$marginal)
    } else {
        stop("Unknown type")
    }

    return(x)
}
