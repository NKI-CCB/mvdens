library(VineCopula)
#source("marginals.r")

fit.vine.copula <- function(x, marginalfn, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), trunclevel = NA, verbose = F) {
    result <- list()
    result$type <- "vine.copula"
    if (!is.null(marginalfn)) {
        result$marginal <- marginalfn(x, bounds)
        transformed <- transform.marginals(x, result$marginal)
    } else {
        result$marginal <- NULL
        transformed <- x
    }
    result$RVM <- RVineStructureSelect(transformed, cores = 1, trunclevel = trunclevel, progress = verbose)
    result$RVM$names <- colnames(x)
    return(structure(result, class = "mdd.density"))
}

evaluate.vine.copula <- function(fit, x, log = F) {
    if (!is.null(fit$marginal)) {
        transformed <- transform.marginals(x, fit$marginal)
    } else {
        transformed <- x
    }
    p_vc <- RVineLogLik(transformed, fit$RVM, separate = T)$loglik
    if (!is.null(fit$marginal)) {
        p <- marginal.correct.p(fit$marginal, x, p_vc, log = T)
    }
    if (log) {
        return(p)
    } else {
        return(exp(p))
    }
}