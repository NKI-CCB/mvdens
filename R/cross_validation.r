#source("kde.r")
#source("gmm.r")
#source("pdf.r")
#source("vine_copula.r")
#source("gp.r")
#source("transform.r")
# 
# source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))
# 
# cwd <- getwd()
# setwd("D:/Research/projects/52-approximate-mc-output-paper/")
# model_simple <- load_sbmlpd_model("1-simple-model", "output_1")
# model_extended <- load_sbmlpd_model("2-extended-model", "output_both")
# setwd(cwd)

#x <- model_extended$posterior$samples[, c(1, 10)]
#bounds <- rbind(c(0, 1), c(-1, 1))
#plot(x)
# 
# subset <- 1:200
# 
# x <- model_extended$posterior$samples[subset,]
# lposterior <- model_extended$posterior$lposterior[subset] + -2.35442
# bounds <- prior_bounds_all(model_extended, q = c(0, 1))
# 
# x <- x
# p <- lposterior
# train.function <- fit.gmm
# K <- 10
# parallel <- F

.cv.iteration <- function(train.function, x.train, x.test, p.train, p.test, log, fit.params)
{
    retval <- list()
    retval$structure <- train.function(x.train, p.train, log, fit.params)
    retval$predicted <- mdd.pdf(retval$structure, x.test, log = log)

    error <- retval$predicted - p.test
    retval$rmse <- sqrt(mean(error * error))

    if (log) {
        error <- retval$predicted - p.test - mdd.marginal.likelihood(retval$structure, x.train, p.train, T)
    } else {
        error <- retval$predicted - p.test / mdd.marginal.likelihood(retval$structure, x.train, p.train, F)
    }
    retval$rmse.corrected <- sqrt(mean(error * error))

    retval$cor.spearman <- cor(retval$predicted, p.test, method = "spearman")

    return(retval)
}

.train.kde <- function(x.train, p.train, log, fit.params) {
    do.call(fit.kde, c(list(x = x.train), fit.params))
}
.train.kde.transformed <- function(x.train, p.train, log, fit.params) {
    do.call(fit.kde.transformed, c(list(x = x.train, p = p.train, log = log), fit.params))
}
.train.gmm <- function(x.train, p.train, log, fit.params) {
    bic <- do.call(gmm.BIC, c(list(x = x.train), fit.params))
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gmm.transformed <- function(x.train, p.train, log, fit.params) {
    bic <- do.call(gmm.transformed.BIC, c(list(x = x.train), fit.params))
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gmm.truncated <- function(x.train, p.train, log, fit.params) {
    bic <- do.call(gmm.truncated.BIC, c(list(x = x.train), fit.params))
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gp <- function(x.train, p.train, log, fit.params) {
    do.call(fit.gp, c(list(x = x.train, p = p.train, log = log), fit.params))
}
.train.vc.ecdf <- function(x.train, p.train, log, fit.params) {
    do.call(fit.vine.copula, c(list(x = x.train, marginalfn = fit.marginal.ecdf), fit.params))
}
.train.vc.ecdf.pareto <- function(x.train, p.train, log, fit.params) {
    do.call(fit.vine.copula, c(list(x = x.train, marginalfn = fit.marginal.ecdf.pareto), fit.params))
}
.train.vc.parametric <- function(x.train, p.train, log, fit.params) {
    do.call(fit.vine.copula, c(list(x = x.train, marginalfn = fit.marginal.parametric), fit.params))
}
.train.vc.mixture <- function(x.train, p.train, log, fit.params) {
    do.call(fit.vine.copula, c(list(x = x.train, marginalfn = fit.marginal.mixture), fit.params))
}

#' Cross validation
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param p
#' @param type
#' @export
#' @examples
mdd.cv <- function(x, p, log, type, nfolds = 10, mcfolds = NULL, mcsize=nrow(x)/10, parallel = F, verbose = F, fit.params=list())
{
    cv <- list()

    train.function <- NULL
    if (type == "kde") {
        train.function <- .train.kde
    } else if (type == "kde.transformed") {
        train.function <- .train.kde.transformed
    } else if (type == "gmm") {
        train.function <- .train.gmm
    } else if (type == "gmm.transformed") {
        train.function <- .train.gmm.transformed
    } else if (type == "gmm.truncated") {
        train.function <- .train.gmm.truncated
    } else if (type == "gp") {
        train.function <- .train.gp
    } else if (type == "vc.ecdf") {
        train.function <- .train.vc.ecdf
    } else if (type == "vc.ecdf.pareto") {
        train.function <- .train.vc.ecdf.pareto
    } else if (type == "vc.parametric") {
        train.function <- .train.vc.parametric
    } else if (type == "vc.mixture") {
        train.function <- .train.vc.mixture
    }

    if (is.null(mcfolds)) {
        reordering <- sample.int(nrow(x))
        holdout_size <- nrow(x) / nfolds

        if (parallel) {
            require(foreach)
            cv$estimates <- foreach(i = 1:nfolds) %dopar% {
                test_ix <- reordering[(i - 1) * holdout_size + 1:holdout_size]
                train_ix <- setdiff(reordering, test_ix)
                return(.cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params))
            }
        } else {
            cv$estimates <- list()
            for (i in 1:nfolds) {
                if (verbose) {
                    cat("Fold =", i, "\n")
                }
                test_ix <- reordering[(i - 1) * holdout_size + 1:holdout_size]
                train_ix <- setdiff(reordering, test_ix)
                cv$estimates[[i]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params)
            }
        }

        cv$predicted <- rep(NA, nrow(x))
        for (i in 1:nfolds) {
            test_ix <- reordering[(i - 1) * holdout_size + 1:holdout_size]
            cv$predicted[test_ix] <- cv$estimates[[i]]$predicted
        }
    } else {
        if (parallel) {
            require(foreach)
            cv$estimates <- foreach(i = 1:mcfolds) %dopar% {
                test_ix <- sample(nrow(x), mcsize)
                train_ix <- setdiff(1:nrow(x), test_ix)
                return(.cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params))
            }
        } else {
            cv$estimates <- list()
            for (i in 1:mcfolds) {
                if (verbose) {
                    cat("Fold =", i, "\n")
                }
                test_ix <- sample(nrow(x), mcsize)
                train_ix <- setdiff(1:nrow(x), test_ix)
                cv$estimates[[i]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params)
            }
        }
    }

    ns <- length(cv$estimates)

    cv$rmse <- rep(NA, ns)
    cv$rmse.corrected <- rep(NA, ns)
    cv$cor.spearman <- rep(NA, ns)
    for (i in 1:ns) {
        cv$rmse[i] <- cv$estimates[[i]]$rmse
        cv$rmse.corrected[i] <- cv$estimates[[i]]$rmse.corrected
        cv$cor.spearman[i] <- cv$estimates[[i]]$cor.spearman
    }
  
  return(cv)
}

cv.compare.plots <- function(cv_list, ...)
{
    require(beeswarm)

    rmse <- lapply(cv_list, function(x) { (x$rmse) })
    beeswarm(rmse, ylab = "root mean squared error", las = 2, ...)

    correlation <- lapply(cv_list, function(x) { x$cor_spearman })
    beeswarm(correlation, ylab = "Spearman correlation", las = 2, ...)
}
