library(VineCopula)
#source("marginals.r")

#' Fit a vine copula
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
#' @examples
fit.vine.copula <- function(x, marginalfn, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), trunclevel = NA, indeptest = F, verbose = F) {
    result <- list()
    result$type <- "vine.copula"
    if (!is.null(marginalfn)) {
        result$marginal <- marginalfn(x, bounds)
        transformed <- marginal.transform(x, result$marginal)
    } else {
        result$marginal <- NULL
        transformed <- x
    }
    cat("indep test", indeptest, "\n")
    result$RVM <- RVineStructureSelect(transformed, cores = 1, trunclevel = trunclevel, indeptest = indeptest, progress = verbose)
    result$RVM$names <- colnames(x)
    return(structure(result, class = "mdd.density"))
}

#' Evaluate a vine copula
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
#' @examples
evaluate.vine.copula <- function(fit, x, log = F) {
    if (!is.null(fit$marginal)) {
        transformed <- marginal.transform(x, fit$marginal)
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
