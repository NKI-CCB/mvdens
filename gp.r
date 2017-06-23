
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

.gp.log.marginal.likelihood <- function(l, result, verbose) {
    kxx <- result$kernel(result$x, result$x, l) + result$sigman * diag(result$n)
    L <- chol(kxx)
    p1 <- sum(backsolve(L, result$p, transpose = T) ^ 2)
    logdet <- 2.0 * sum(log(diag(L)))
    logml <- -0.5 * p1 - 0.5 * logdet - 0.5 * result$n * log2(pi)
    if (verbose) {
        cat(c("l:", l, "logml:",logml, "\n"))
    }
    return(logml)
}

fit.gp <- function(x, p, kernel, l=1, sigman=1e-10, verbose=T) {
    result <- list()
    result$type <- "gp"

    if (is.matrix(x)) {
        result$n <- nrow(x)
    } else {
        result$n <- length(x)
    }

    if (kernel == "squared.exponential" || kernel == "se") {
        result$kernel <- .kernel.squared.exponential
    } else if (kernel == "matern32") {
        result$kernel <- .kernel.matern32
    }
    result$x <- x
    result$p <- p
    result$sigman <- sigman

    if (length(l) > 1) {
        # Kernel length range, find optimal marginal likelihood
        opt <- optimize(.gp.log.marginal.likelihood, l, maximum = T, result=result, verbose=verbose)
        result$l <- opt$maximum

        #ctrl <- list()
        #ctrl$fnscale <- -1.0
        #ctrl$reltol <- 1e-6
        #opt <- optim(c(1, 1e-8), .gp.log.marginal.likelihood, result = result, verbose = verbose, method = "Nelder-Mead", lower=c(0.01, 1e-12), upper = c(100,0.1), control = ctrl)

    } else {
        result$l <- l
    }

    result$kxx <- result$kernel(x, x, result$l) + 1e-12 * diag(result$n)
    result$L <- chol(result$kxx)
    result$alpha <- backsolve(result$L, backsolve(result$L, p, transpose = TRUE))

    return(result)
}

evaluate.gp <- function(fit, x) {
    stopifnot(fit$type == "gp")

    kxsx <- fit$kernel(x, fit$x, fit$l)
    f <- kxsx %*% fit$alpha

    return(f)
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