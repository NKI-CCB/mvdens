#' Fit a vine copula
#'
#' Fit a vine copula density function, using one of the marginal function for the marginal transformations, and using VineCopula::RVineStructureSelect for the vine copula function.
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param marginalfn One of the marginal fitting functions.
#' @param bounds Dx2 matrix of boundaries
#' @param trunclevel Relayed to the VineCopula::RVineStructureSelect 
#' @param indeptest Relayed to VineCopula::RVineStructureSelect 
#' @param verbose Relayed to VineCopula::RVineStructureSelect 
#' @export
fit.vine.copula <- function(x, marginalfn, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), trunclevel = NA, indeptest = F, verbose = F) {
  result <- list()
  result$type <- "vine.copula"
    if (!is.null(marginalfn)) {
        result$marginal <- marginalfn(x, bounds)
        transformed <- mvd.marginal.transform(x, result$marginal)
    } else {
        result$marginal <- NULL
        transformed <- x
    }
  result$RVM <- VineCopula::RVineStructureSelect(transformed, cores = 1, trunclevel = trunclevel, indeptest = indeptest, progress = verbose)
  result$RVM$names <- colnames(x)
  return(structure(result, class = "mvd.density"))
}

#' Evaluate a vine copula
#'
#' description
#' @param fit An mvd.density object obtained from fit.vine.copula
#' @param x Matrix or vector of samples at which to evaluate the vine copula. For matrices, rows are samples and columns are variables.
#' @param log Return log probability density
#' @export
evaluate.vine.copula <- function(fit, x, log = F) {
  if (!is.null(fit$marginal)) {
    transformed <- mvd.marginal.transform(x, fit$marginal)
  } else {
    transformed <- x
  }
  p_vc <- VineCopula::RVineLogLik(transformed, fit$RVM, separate = T)$loglik
  if (!is.null(fit$marginal)) {
    p <- .marginal.correct.p(fit$marginal, x, p_vc, log = T)
  }
  if (log) {
      return(p)
    } else {
      return(exp(p))
    }
}
