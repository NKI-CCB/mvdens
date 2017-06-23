library(VineCopula)
source("marginals.r")

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
    result$RVM <- RVineStructureSelect(transformed, indeptest = T, cores = 1, trunclevel = trunclevel, progress = verbose)
    return(structure(result, class = "mdd.vine.copula"))
}

evaluate.vine.copula <- function(fit, x) {
    if (!is.null(fit$marginal)) {
        transformed <- transform.marginals(x, fit$marginal)
    } else {
        transformed <- x
    }
    p_vc <- RVineLogLik(transformed, fit$RVM, separate = T)$loglik
    if (!is.null(fit$marginal)) {
        p <- marginal.correct.p(fit$marginal, x, p_vc)
    }
    return(p)
}
