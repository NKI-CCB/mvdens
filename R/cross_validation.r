
.cv.iteration <- function(train.function, x.train, x.test, p.train, p.test, log, fit.params)
{
    retval <- list()
    retval$structure <- train.function(x.train, p.train, log, fit.params)
    retval$predicted <- mvd.pdf(retval$structure, x.test, log = log)

    error <- retval$predicted - p.test
    retval$rmse <- sqrt(mean(error * error))

    if (log) {
        error <- retval$predicted - p.test - mvd.marginal.likelihood(retval$structure, x.train, p.train, T)
    } else {
        error <- retval$predicted - p.test / mvd.marginal.likelihood(retval$structure, x.train, p.train, F)
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
.train.gaussian <- function(x.train, p.train, log, fit.params) {
    do.call(fit.gmm, c(list(x = x.train, K = 1), fit.params))
}
.train.gmm <- function(x.train, p.train, log, fit.params) {
    aic <- do.call(gmm.AIC, c(list(x = x.train), fit.params))
    best <- which.min(aic$AIC)
    return(aic$fits[[best]])
}
.train.gmm.transformed <- function(x.train, p.train, log, fit.params) {
    aic <- do.call(gmm.transformed.AIC, c(list(x = x.train), fit.params))
    best <- which.min(aic$AIC)
    return(aic$fits[[best]])
}
.train.gmm.truncated <- function(x.train, p.train, log, fit.params) {
    aic <- do.call(gmm.truncated.AIC, c(list(x = x.train), fit.params))
    best <- which.min(aic$AIC)
    return(aic$fits[[best]])
}
.train.gp <- function(x.train, p.train, log, fit.params) {
  do.call(fit.gp, c(list(x = x.train, p = p.train), fit.params))
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
.train.factor.mixture <- function(x.train, p.train, log, fit.params) {
  do.call(factor.mixture.AIC, c(list(x = x.train, optimal.only = T), fit.params))
}

#' Estimate the approximation accuracy of a multivariate density using cross validation
#'
#' Perform cross validation to estimate the accuracy of the density approximations. This requires the known probability densities at the sample points, though these probability densities do not have to be normalized. If they are not normalized, then only the correlation is a meaningful performance measure. If they are normalized, the RMSE is also informative.
#' Either K-fold cross validation or Monte Carlo cross validation can be performed.
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param p Vector of posterior probability densities of the samples.
#' @param log Boolean which specifies whether p is in log space.
#' @param type String describing the type of the density approximation to use, can be: kde, kde.transformed, gmm, gmm.transformed, gmm.truncated, gp, vc.ecdf, vc.ecdf.pareto, vc.parametric, vc.mixture
#' @param cvtype String describing the type of cross validation to use, can be either kfoldcv, mccv or fixed
#' @param nfolds The number of cross validation folds (i.e. K when using K-fold cross validation, or the number of iterations in Monte Carlo cross validation).
#' @param mcsize When using Monte Carlo cross validation, this specifies the size of the test set.
#' @param parallel When true, and when parallel workers have been registered, run the cross validation in parallel using foreach. The user should start the parallel workers (see examples).
#' @param verbose Display the progress (when parallel=F).
#' @param fit.params A named list with parameters to be relayed to the density fitting functions.
#' @export
#' @examples
#' mu <- c(0,0)
#' sigma <- rbind(c(1,0.5),c(0.5,1))
#' x <- mvtnorm::rmvnorm(40, mu, sigma)
#' p <- mvtnorm::dmvnorm(x, mu, sigma)
#' res <- mvd.cv(x, p, 'gmm', log=FALSE, cvtype='kfoldcv', nfolds=4, verbose=TRUE)
#' 
#' #library(doParallel)
#' #cl <- parallel::makeCluster(6)
#' #doParallel::registerDoParallel(cl)
#' #parallel::clusterEvalQ(cl, library(mvdens))
#' #res <- mvd.cv(x, p, 'gmm', log=FALSE, cvtype='kfoldcv', nfolds=4, parallel=TRUE)

mvd.cv <- function(x, p, log, type, cvtype, nfolds = 10, mcsize = nrow(x) / 10, parallel = F, verbose = F, fit.params = list()) {
  cv <- list()
  
  train.function <- NULL
  if (type == "kde") {
    train.function <- .train.kde
  } else if (type == "kde.transformed") {
    train.function <- .train.kde.transformed
  } else if (type == "gaussian") {
    train.function <- .train.gaussian
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
  } else if (type == "mfa") {
    train.function <- .train.factor.mixture
  }
  `%dopar%` <- foreach::`%dopar%`
  
  if (cvtype == "kfoldcv") {
    reordering <- sample.int(nrow(x))
    holdout_size <- nrow(x) / nfolds
    
    if (parallel) {
      cv$estimates <- foreach::foreach(i = 1:nfolds) %dopar% {
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
    } else if(cvtype == "mccv") {
        if (parallel) {
            cv$estimates <- foreach::foreach(i = 1:nfolds) %dopar% {
                test_ix <- sample(nrow(x), mcsize)
                train_ix <- setdiff(1:nrow(x), test_ix)
                return(.cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params))
            }
        } else {
          cv$estimates <- list()
          for (i in 1:nfolds) {
            if (verbose) {
              cat("Fold =", i, "\n")
            }
            test_ix <- sample(nrow(x), mcsize)
            train_ix <- setdiff(1:nrow(x), test_ix)
            cv$estimates[[i]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params)
          }
        }
    } else if (cvtype == "fixed") {
      cv$estimates <- list()
      train_ix <- 1:(nrow(x) - mcsize)
      test_ix <- (nrow(x) - mcsize + 1):nrow(x)
      cv$estimates[[1]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], log = log, fit.params)
    } else {
      stop("Invalid cvtype")
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

# mvd.cv.compare.plots <- function(cv_list, ...)
# {
#     require(beeswarm)
# 
#     rmse <- lapply(cv_list, function(x) { (x$rmse) })
#     beeswarm(rmse, ylab = "root mean squared error", las = 2, ...)
# 
#     correlation <- lapply(cv_list, function(x) { x$cor_spearman })
#     beeswarm(correlation, ylab = "Spearman correlation", las = 2, ...)
# }
