.kernel.squared.exponential <- function(x1, x2, l) {
    factor <- -1.0 / (2 * l * l)
    if (is.matrix(x1)) {
        sigma <- matrix(0, nrow(x1), nrow(x2))
        for (i in 1:nrow(sigma)) {
            v <- t(x1[i,] - t(x2))
            rsq <- apply(v * v, 1, sum)
            sigma[i,] <- exp(factor * rsq)
            #y <- sqrt(3) * sqrt(rsq) / l
            #sigma[i,j] <- (1 + y) * exp(-y)
        }
    } else {
        sigma <- matrix(0, length(x1), length(x2))
        for (i in 1:nrow(sigma)) {
            v <- x1[i] - x2
            sigma[i,] <- exp(factor * v * v)
        }
    }

    return(sigma)
}

.kernel.se.integral <- function(alpha, l, D) {
    d <- sqrt((2*pi*l*l)^D)
    integral <- sum(d * alpha)
    return(integral)
}

.kernel.matern32 <- function(x1, x2, l) {
    if (is.matrix(x1)) {
        sigma <- matrix(0, nrow(x1), nrow(x2))
        for (i in 1:nrow(sigma)) {
            v <- t(x1[i,] - t(x2))
            rsq <- apply(v * v, 1, sum)
            y <- sqrt(3) * sqrt(rsq) / l
            sigma[i,] <- (1 + y) * exp(-y)
        }
    } else {
        sigma <- matrix(0, length(x1), length(x2))
        for (i in 1:nrow(sigma)) {
            v <- x1[i] - x2
            y <- sqrt(3) * sqrt(v * v) / l
            sigma[i, ] <- (1 + y) * exp(-y)
        }
    }

    return(sigma)
}

.kernel.matern32.integral <- function(alpha, l, D) {
    d <- 2 * (1+D) * (l * sqrt(pi/3)) ^ D * gamma(D) / gamma(D / 2)
    integral <- sum(d * alpha)
    return(integral)
}

.kernel.matern52 <- function(x1, x2, l) {
    if (is.matrix(x1)) {
        sigma <- matrix(0, nrow(x1), nrow(x2))
        for (i in 1:nrow(sigma)) {
            v <- t(x1[i,] - t(x2))
            rsq <- apply(v * v, 1, sum)
            y <- sqrt(5) * sqrt(rsq) / l
            sigma[i,] <- (1 + y + 5 * rsq / (3 * l * l)) * exp(-y)
        }
    } else {
        sigma <- matrix(0, length(x1), length(x2))
        for (i in 1:nrow(sigma)) {
            v <- x1[i] - x2
            y <- sqrt(3) * sqrt(v * v) / l
            sigma[i,] <- (1 + y + 5 * rsq / (3 * l * l)) * exp(-y)
        }
    }

    return(sigma)
}

.gp.log.marginal.likelihood <- function(l, b1, b2, result, verbose) {
    kxx <- result$kernel(result$x, result$x, l) + result$sigman * diag(result$n)
    
    return_value <- try(L <- chol(kxx), silent = !verbose)
    if (inherits(return_value, "try-error")) {
      return(NA)
    }

    if (is.null(result$meanfn)) {
        mean_sub <- rep(NA, length(result$p))
        for (i in 1:length(result$p)) {
            mean_sub[i] <- result$p[i] - b1 - b2 * (result$xcentered[i,] %*% result$xcentered[i,])
        }
        p1 <- sum(backsolve(L, mean_sub, transpose = T) ^ 2)
    } else {
        p1 <- sum(backsolve(L, result$pmean_sub, transpose = T) ^ 2)
    }
    logdet <- 2.0 * sum(log(diag(L)))
    logml <- -0.5 * p1 - 0.5 * logdet - 0.5 * result$n * log(2*pi)
    if (verbose) {
        cat(c("l:", l, "b1:", b1, "b2:", b2, "logml:",logml, "\n"))
    }
    return(logml)
}

.gp.log.marginal.likelihood.optiml <- function(x, result, verbose) {
    return(.gp.log.marginal.likelihood(x, result$b1, result$b2, result, verbose))
}
.gp.log.marginal.likelihood.optimb1 <- function(x, result, verbose) {
    return(.gp.log.marginal.likelihood(result$l, x, result$b2, result, verbose))
}

.gp.cv.optimize <- function(l, result, verbose) {
    nfolds <- 5
    n <- result$n
    if (is.matrix(result$x)) {
      D <- ncol(result$x)
    } else {
      D <- 1
    }

    sse <- rep(NA, nfolds)
    for (fi in 1:nfolds) {
        foldsize <- n / nfolds
        test_ix <- (fi - 1) * foldsize + 1:foldsize
        train_ix <- setdiff(1:n, test_ix)

        if (is.matrix(result$x)) {
          xtrain <- result$x[train_ix,]
          xtest <- result$x[test_ix,]
        } else {
          xtrain <- result$x[train_ix]
          xtest <- result$x[test_ix]
        }
        nt <- length(train_ix)
        K <- result$kernel(xtrain, xtrain, l)

        return_value <- try(L <- chol(K + result$sigman * diag(nt)), silent = !verbose)
        if (inherits(return_value, "try-error")) {
            return(Inf)
        }

        alpha <- backsolve(L, backsolve(L, result$p[train_ix], transpose = T))

        ntest <- length(test_ix)
        ktest <- result$kernel(xtest, xtrain, l)

        integral <- result$kernel.integral(alpha, l, D)
        f <- (ktest %*% alpha)
        diff <- result$p[test_ix] - f
        sse[fi] <- sum(diff ^ 2)
        if (integral < 0) {
            sse[fi] <- sse[fi] * 1e10
        }

    }
    rmse <- sqrt(sum(sse) / n)

    if (verbose) {
         cat(c("l:", l, "rmse:", rmse, "\n"))
    }

    return(rmse)
}

#' Fit a Gaussian process density function
#'
#' Fit a Gaussian process through the posterior density samples. The GP can optionally be normalized such that the predictive function integrates to 1.
#' The GP assumes a zero mean, and can use either a squared exponential or Matern-type covariance function. A single length scale parameter is used, which is
#' optimized by minimizing the RMSE in 5-fold cross validation.
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param p Probability density at the sample locations; does not have to be normalized. Note that the probabilities should not be log transformed.
#' @param kernel Which kernel to use, can be "se", "matern32" or "matern52".
#' @param normalize Boolean specifying whether to normalize the GP such that it integrates to 1.
#' @param sigman Value added to the diagonal of the kernel matrix; can help numerical stability.
#' @param verbose Display progress during the optimization of l.
#' @export
#' @examples
fit.gp <- function(x, p, kernel, l = 1.0, optimize = T, normalize = T, sigman = 1e-12, verbose = F) {
    result <- list()
    result$type <- "gp"

    if (is.matrix(x)) {
        result$n <- nrow(x)
    } else {
        result$n <- length(x)
    }

    result$kernel.name <- kernel
    if (kernel == "squared.exponential" || kernel == "se") {
        result$kernel <- .kernel.squared.exponential
        result$kernel.integral <- .kernel.se.integral
    } else if (kernel == "matern32") {
        result$kernel <- .kernel.matern32
        result$kernel.integral <- .kernel.matern32.integral
    } else if (kernel == "matern52") {
        result$kernel <- .kernel.matern52
        result$kernel.integral <- .kernel.matern52.integral
    }

    result$x <- x
    result$p <- p
    result$sigman <- sigman
    result$log <- F

    if (optimize) {
        #if (isotropic) {
            stopifnot(length(l) > 1)
            opt <- optimize(.gp.cv.optimize, l, result = result, verbose = verbose)
            result$l <- opt$minimum
        #} else {
        #    require(nloptr)
        #    stopifnot(ncol(l) == 3 && nrow(l) == ncol(x))
        #    opt <- sbplx(l[, 1], .gp.cv.optimize, l[, 2], l[, 3], result = result, verbose = verbose)
        #    result$l <- opt$par
        #}
    } else {
        result$l <- l
    }

    result$kxx <- result$kernel(result$x, result$x, result$l) + result$sigman * diag(result$n)
    L <- chol(result$kxx)
    result$alpha <- backsolve(L, backsolve(L, p, transpose = TRUE))

    if (normalize) {
        if (is.matrix(result$x)) {
            D <- ncol(result$x)
        } else {
            D <- 1
        }
        integral <- result$kernel.integral(result$alpha, result$l, D)
        result$s <- 1.0 / integral
    } else {
        result$s <- 1.0
    }

    return(structure(result, class = "mvd.density"))
}

#' Evaluate a Gaussian process density function
#'
#' description
#' @param fit mvd.density object
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
#' @examples
evaluate.gp <- function(fit, x) {
    stopifnot(fit$type == "gp")

    if (fit$log) {
    } else {
        kxsx <- fit$kernel(x, fit$x, fit$l)
        f <- kxsx %*% fit$alpha
        return(fit$s * f)
    }
}
