source("kde.r")
source("gmm.r")
source("pdf.r")
source("vine_copula.r")
source("gp.r")
source("transform.r")
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

.cv.iteration <- function(train.function, x.train, x.test, p.train, p.test, ...)
{
  retval <- list()
  
  retval$structure <- train.function(x.train, p.train, ...)
  retval$predicted <- mdd.pdf(retval$structure, x.test, log=T)
  
  error <- retval$predicted - p.test
  retval$MSEabs <- mean(error * error)

  error <- (retval$predicted - mean(retval$predicted)) - (p.test - mean(p.test))
  retval$MSE <- mean(error * error)

  retval$cor_pearson <- cor(retval$predicted, p.test, method="pearson")
  retval$cor_spearman <- cor(retval$predicted, p.test, method="spearman")
  
  return(retval)
}

.train.kde <- function(x.train, p.train, ...) {
    return(fit.kde(x.train, ...))
}
.train.kde.transformed <- function(x.train, p.train, ...) {
    return(fit.kde.transformed(x.train, ...))
}
.train.gmm <- function(x.train, p.train, ...) {
    bic <- gmm.BIC(x.train, ...)
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gmm.transformed <- function(x.train, p.train, ...) {
    bic <- gmm.transformed.BIC(x.train, ...)
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gmm.truncated <- function(x.train, p.train, ...) {
    bic <- gmm.truncated.BIC(x.train, ...)
    best <- which.min(bic$BIC)
    return(bic$fits[[best]])
}
.train.gp <- function(x.train, p.train, ...) {
    return(fit.gp(x=x.train, p=p.train, ...))
}
.train.vc.ecdf <- function(x.train, p.train, ...) {
    return(fit.vine.copula(x.train, fit.marginal.ecdf, ...))
}
.train.vc.parametric <- function(x.train, p.train, ...) {
    return(fit.vine.copula(x.train, fit.marginal.parametric, ...))
}
.train.vc.mixture <- function(x.train, p.train, ...) {
    return(fit.vine.copula(x.train, fit.marginal.mixture, ...))
}

mdd.cv <- function(x, p, type, nfolds = 10, mcfolds = NULL, mcsize=nrow(x)/10, parallel = F, verbose = F, ...)
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
                return(.cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], verbose = verbose, ...))
            }
        } else {
            cv$estimates <- list()
            for (i in 1:nfolds) {
                if (verbose) {
                    cat("Fold =", i, "\n")
                }
                test_ix <- reordering[(i - 1) * holdout_size + 1:holdout_size]
                train_ix <- setdiff(reordering, test_ix)
                cv$estimates[[i]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], verbose = verbose, ...)
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
                return(.cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], verbose = verbose, ...))
            }
        } else {
            cv$estimates <- list()
            for (i in 1:mcfolds) {
                if (verbose) {
                    cat("Fold =", i, "\n")
                }
                test_ix <- sample(nrow(x), mcsize)
                train_ix <- setdiff(1:nrow(x), test_ix)
                cv$estimates[[i]] <- .cv.iteration(train.function, x[train_ix,], x[test_ix,], p[train_ix], p[test_ix], verbose = verbose, ...)
            }
        }
    }

    ns <- length(cv$estimates)

    cv$MSEabs <- rep(NA, ns)
    cv$MSE <- rep(NA, ns)
    cv$cor_pearson <- rep(NA, ns)
    cv$cor_spearman <- rep(NA, ns)
    for (i in 1:ns) {
        cv$MSEabs[i] <- cv$estimates[[i]]$MSEabs
        cv$MSE[i] <- cv$estimates[[i]]$MSE
        cv$cor_pearson[i] <- cv$estimates[[i]]$cor_pearson
        cv$cor_spearman[i] <- cv$estimates[[i]]$cor_spearman
    }
  
  return(cv)
}

cv.compare.plots <- function(cv_list)
{
    require(beeswarm)

    mse <- lapply(cv_list, function(x) { log(x$MSE) })
    beeswarm(mse, ylab = "log mean squared error", las = 2)

    correlation <- lapply(cv_list, function(x) { x$cor_spearman })
    beeswarm(correlation, ylab = "Spearman correlation", las = 2)
}

# 
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# result <- mdd.cv(x, lposterior, "gmm", K = 4, parallel = T)
# result_truncated <- mdd.cv(x, lposterior, "gmm.truncated", K = 3, bounds = bounds, parallel = T)
# result_gp <- mdd.cv(x, lposterior - min(lposterior) + 3, "gp", kernel = "se", l = 0.68, parallel = T)
# result_gp_matern <- mdd.cv(x, lposterior - min(lposterior) + 3, "gp", kernel = "matern32", l = 0.93, parallel = T)
# result_gp_exp <- mdd.cv(x, exp(lposterior), "gp", kernel = "se", l = 0.26, parallel = T)
# 
# result_vc_ecdf <- mdd.cv(x, lposterior, "vc.ecdf", bounds = bounds, parallel = T)
# result_vc_parametric <- mdd.cv(x, lposterior, "vc.parametric", bounds = bounds, parallel = T)
# result_vc_mixture <- mdd.cv(x, lposterior, "vc.mixture", bounds = bounds, parallel = T)
# 
# stopCluster(cl)
# 
# par(mfrow=c(2,4))
# plot(lposterior, result$predicted, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20,20), c(-20,20), lty=2)
# plot(lposterior, result_truncated$predicted, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, result_gp$predicted + min(lposterior) - 3, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, result_gp_matern$predicted + min(lposterior) - 3, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, log(result_gp_exp$predicted), xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, result_vc_ecdf$predicted, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, result_vc_parametric$predicted, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# plot(lposterior, result_vc_mixture$predicted, xlim = c(-10, 15), ylim = c(-15, 15))
# lines(c(-20, 20), c(-20, 20), lty = 2)
# par(mfrow = c(1, 1))
# 
# par(mar=c(12,4,4,1))
# beeswarm(lapply(list(result$MSE, result_truncated$MSE, result_gp$MSE, result_gp_matern$MSE, result_gp_exp$MSE, result_vc_ecdf$MSE, result_vc_parametric$MSE, result_vc_mixture$MSE), log),
#          labels = c("GMM", "Truncated GMM", "Gaussian process", "Gaussian process matern", "Gaussian process exp", "Vine Copula - ECDF", "Vine Copula - parametric", "Vine Copula - mixture"),
#          ylab = "log mean squared error", las = 2)
# beeswarm(list(result$cor_spearman, result_truncated$cor_spearman, result_gp$cor_spearman, result_gp_matern$cor_spearman, result_gp_exp$cor_spearman, result_vc_ecdf$cor_spearman, result_vc_parametric$cor_spearman, result_vc_mixture$cor_spearman),
#          labels = c("GMM", "Truncated GMM", "Gaussian process", "Gaussian process matern", "Gaussian process", "Vine Copula - ECDF", "Vine Copula - parametric", "Vine Copula - mixture"),
#          ylab = "Spearman correlation", las = 2)
# par(mar = c(5, 4, 4, 1))
