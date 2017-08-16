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
.gp.log.marginal.likelihood.optimb2 <- function(x, result, verbose) {
    return(.gp.log.marginal.likelihood(result$l, result$b1, x, result, verbose))
}
.gp.log.marginal.likelihood.optim2 <- function(x, result, verbose) {
    return(.gp.log.marginal.likelihood(x[1], x[2], result$b2, result, verbose))
}
.gp.log.marginal.likelihood.optim3 <- function(x, result, verbose) {
    return(.gp.log.marginal.likelihood(x[1], x[2], x[3], result, verbose))
}

#' Fit a Gaussian process density function
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param p
#' @param kernel
#' @export
#' @examples
fit.gp <- function(x, p, kernel, l=1, b1=0, b2=0, meanfn=NULL, sigman=1e-10, verbose=F) {
    result <- list()
    result$type <- "gp"

    if (is.matrix(x)) {
        result$n <- nrow(x)
        result$xmean <- apply(x, 2, mean)
    } else {
        result$n <- length(x)
        result$xmean <- mean(x)
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
    result$xcentered <- t(t(x) - result$xmean)
    result$p <- p
    result$sigman <- sigman
    result$meanfn <- meanfn

    if (!is.null(meanfn)) {
        result$pmean_sub <- result$p - meanfn(result$x)
    }

    num_parameters_range <- (length(l) > 1) + (length(b1) > 1) + (length(b2) > 1)

    if (num_parameters_range == 0) {
        result$l <- l
        result$b1 <- b1
        result$b2 <- b2
    } else if (num_parameters_range == 1) {
        result$l <- l
        result$b1 <- b1
        result$b2 <- b2
        if (length(l) > 1) {
            opt <- optimize(.gp.log.marginal.likelihood.optiml, l, maximum = T, result = result, verbose = verbose)
            result$l <- opt$maximum
        } else if (length(b1) > 1) {
            opt <- optimize(.gp.log.marginal.likelihood.optimb1, b1, maximum = T, result = result, verbose = verbose)
            result$b1 <- opt$maximum
        } else if (length(b1) > 1) {
            opt <- optimize(.gp.log.marginal.likelihood.optimb2, b2, maximum = T, result = result, verbose = verbose)
            result$b2 <- opt$maximum
        }
    } else if (num_parameters_range == 2) {
        if (length(b2) != 1) {
            stop("Not implemented optimization for b2 without l and b1")
        }
        result$b2 <- b2
        ctrl <- list()
        ctrl$fnscale <- -1.0
        opt <- optim(c(1, mean(p)), .gp.log.marginal.likelihood.optim2, result = result, verbose = verbose, method = "Nelder-Mead", control = ctrl)
        result$l <- opt$par[1]
        result$b1 <- opt$par[2]
    } else if (num_parameters_range == 3) {
        ctrl <- list()
        ctrl$fnscale <- -1.0
        opt <- optim(c(1, mean(p), -1), .gp.log.marginal.likelihood.optim3, result = result, verbose = verbose, method = "Nelder-Mead", control = ctrl)
        result$l <- opt$par[1]
        result$b1 <- opt$par[2]
        result$b2 <- opt$par[3]
    } else {
        stop()
    }

    result$kxx <- result$kernel(result$x, result$x, result$l) + result$sigman * diag(result$n)
    L <- chol(result$kxx)

    if (is.null(meanfn)) {
        result$pmean_sub <- rep(NA, length(result$p))
        for (i in 1:length(result$p)) {
            result$pmean_sub[i] <- result$p[i] - result$b1 - result$b2 * (result$xcentered[i,] %*% result$xcentered[i,])
        }
    }

    result$alpha <- backsolve(L, backsolve(L, result$pmean_sub, transpose = TRUE))

    return(structure(result, class = "mdd.density"))
}

.gp.cv.optimize <- function(l, result, verbose, scale_to_1) {
    nfolds <- 5
    n <- result$n
    if (is.matrix(result$x)) {
      D <- ncol(result$x)
    } else {
      D <- 1
    }
    
    K <- result$kernel(result$x, result$x, l)

    return_value <- try(L <- chol(K + result$sigman * diag(n)), silent = !verbose)
    if (inherits(return_value, "try-error")) {
        return(Inf)
    }

    if (scale_to_1) {
        alpha <- backsolve(L, backsolve(L, result$p, transpose = T))
        integral <- result$kernel.integral(alpha, l, D)
        if (integral <= 0) {
            cat(c("l:", l, "integral:", integral, "\n"))
            return(Inf)
        }
        s <- 1.0 / integral
    } else {
        s <- 1.0
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

        f <- (ktest %*% alpha) * s
        diff <- result$p[test_ix] - f
        sse[fi] <- sum(diff ^ 2)
    }
    rmse <- sqrt(sum(sse) / n)

    if (verbose) {
        if (scale_to_1) {
            cat(c("l:", l, "integral:", integral, "rmse:", rmse, "\n"))
        } else {
            cat(c("l:", l, "rmse:", rmse, "\n"))
        }
    }

    return(rmse)
}

#' Fit a Gaussian process density function
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param p
#' @param kernel
#' @export
#' @examples
fit.gp2 <- function(x, p, kernel, l = 1.0, scale_to_1=T, sigman = 1e-12, verbose = F) {
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

    if (length(l) > 1) {
        library(nloptr)
        opt <- neldermead(mean(l), .gp.cv.optimize, lower = l[1], upper=l[2], result = result, verbose = verbose, scale_to_1 = scale_to_1, control = nl.opts(list(xtol_rel = 1e-6, ftol_rel = 1e-6, ftol_abs = 1e-6)))
        #opt <- sbplx(l0, .gp.cv.optimize, lower=0, result = result, verbose = verbose)
        result$l <- opt$par
        #opt <- optimize(.gp.cv.optimize, c(0.1,2), result = result, verbose = verbose)
        #result$l <- opt$minimum
    } else {
        result$l <- l
    }

    result$kxx <- result$kernel(result$x, result$x, result$l) + result$sigman * diag(result$n)
    L <- chol(result$kxx)
    result$alpha <- backsolve(L, backsolve(L, p, transpose = TRUE))

    if (scale_to_1) {
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

    return(structure(result, class = "mdd.density"))
}
#' Fit a Gaussian process density function
#'
#' description
#' @param fit
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
#' @examples
evaluate.gp <- function(fit, x) {
    stopifnot(fit$type == "gp")

    xcenter <- t(t(x) - fit$xmean)
    kxsx <- fit$kernel(x, fit$x, fit$l)
    f <- kxsx %*% fit$alpha
    if (is.null(fit$meanfn)) {
        for (i in 1:nrow(xcenter)) {
            f[i] <- f[i] + fit$b1 + fit$b2 * (xcenter[i,] %*% xcenter[i,])
        }
    } else {
        f <- f + fit$meanfn(x)
    }

    return(f)
}

#' Evaluate a Gaussian process density function
#'
#' description
#' @param fit
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
#' @examples
evaluate.gp2 <- function(fit, x) {
    stopifnot(fit$type == "gp")
    kxsx <- fit$kernel(x, fit$x, fit$l)
    f <- kxsx %*% fit$alpha
    return(fit$s * f)
}

#x <- rnorm(1000)
#p <- dnorm(x)
#xstar <- seq(-10, 10, length.out = 200)

#fit <- fit.gp(x, p, "se", l = c(0.01, 100))
#fstar <- evaluate.gp(fit, xstar)

#plot(xstar, fstar)

#l <- 1.5

#kxx <- .kernel.squared.exponential(x, x, l)
#kxsx <- .kernel.squared.exponential(xstar, x, l)

#fstar <- kxsx %*% solve(kxx + 1e-12 * diag(nrow(kxx))) %*% p
#plot(xstar, fstar, typ="l", ylim=c(-2,2))

#kxxs <- .kernel.squared.exponential(x, xstar, l)
#kxsxs <- .kernel.squared.exponential(xstar, xstar, l)
#covfstar <- kxsxs - kxsx %*% solve(kxx + 1e-8 * diag(nrow(kxx))) %*% kxxs

#polygon(c(xstar, rev(xstar)), c(fstar + diag(covfstar), rev(fstar - diag(covfstar))), col = "grey", border = NA)


#kxsx <- fit$kernel(xstar, fit$x, fit$l)
#fstar <- kxsx %*% solve(fit$kxx) %*% p

#head(sort(eigen(kxx + 1e-13 * diag(nrow(kxx)))$values))