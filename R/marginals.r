#' Marginal function for use in fit.vine.copula: empirical distribution function
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
fit.marginal.ecdf <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), reflect_bounds = T) {
    marginal <- list()
    marginal$type <- "ecdf"

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    marginal$ecdfs <- apply(x, 2, ecdf)
    marginal$varnames <- colnames(x)
    marginal$x <- x
    marginal$x_reflected <- list()
    marginal$correction_factor <- rep(1, ncol(x))
    for (i in 1:ncol(x)) {
        marginal$x_reflected[[i]] <- x[, i]
        if (reflect_bounds) {
            if (bounds[i, 1] != -Inf) {
                reflected <- bounds[i, 1] - (x[, i] - bounds[i, 1])
                marginal$x_reflected[[i]] <- c(marginal$x_reflected[[i]], reflected)
                marginal$correction_factor[i] <- marginal$correction_factor[i] + 1
            }
            if (bounds[i, 2] != Inf) {
                reflected <- bounds[i, 2] + (bounds[i, 2] - x[, i])
                marginal$x_reflected[[i]] <- c(marginal$x_reflected[[i]], reflected)
                marginal$correction_factor[i] <- marginal$correction_factor[i] + 1
            }
        }
    }

    marginal$bw <- rep(NA, ncol(x))
    for (i in 1:ncol(x)) {
        hmax <- 10 * sd(marginal$x_reflected[[i]])
        marginal$bw[i] <- bw.SJ(marginal$x_reflected[[i]], lower = 1e-6 * hmax, upper = hmax)
    }

    return(structure(marginal, class = "mvd.marginal"))
}

#' Marginal function for use in fit.vine.copula: single parametric distribution
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
fit.marginal.parametric <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x)))) {
    library(MASS)

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    marginal <- list()
    marginal$type <- "parametric"
    marginal$dists <- list()
    fitctrl <- list(maxit = 1000)
    for (i in 1:ncol(x)) {
        if (bounds[i, 1] >= bounds[i, 2]) {
            stop("Upper bound must be strictly higher than lower bound")
        }

        if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
            # [0,1] -> beta
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "beta"
            marginal$dists[[i]]$min <- 0
            marginal$dists[[i]]$max <- 1
            dfit <- fitdistr(x[, i], dbeta, list(shape1 = 1, shape2 = 1), control = fitctrl)
            marginal$dists[[i]]$shape1 <- as.numeric(dfit$estimate["shape1"])
            marginal$dists[[i]]$shape2 <- as.numeric(dfit$estimate["shape2"])
        } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
            # [0,inf] -> gamma
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "gamma"
            rescale <- 1 / mean(x[, i])
            dfit <- fitdistr(x[, i], dgamma, list(shape = 2, scale = 1), control = fitctrl)
            dfit <- fitdistr(x[, i] * rescale, dgamma, list(shape = 2, scale = 1), control = fitctrl)
            marginal$dists[[i]]$shape <- as.numeric(dfit$estimate["shape"])
            marginal$dists[[i]]$scale <- as.numeric(dfit$estimate["scale"]) / rescale
        } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
            # [-inf,inf] -> normal
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "normal"
            marginal$dists[[i]]$mu <- mean(x[, i])
            marginal$dists[[i]]$sigma <- sd(x[, i])
        } else {
            # [a,b] -> scaled beta
            a <- bounds[i, 1];
            b <- bounds[i, 2];
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "beta"
            marginal$dists[[i]]$min <- a
            marginal$dists[[i]]$max <- b
            dfit <- fitdistr((x[, i] - a) / (b - a), dbeta, list(shape1 = 1, shape2 = 1), control = fitctrl)
            marginal$dists[[i]]$shape1 <- as.numeric(dfit$estimate["shape1"])
            marginal$dists[[i]]$shape2 <- as.numeric(dfit$estimate["shape2"])
        }

        marginal$dists[[i]]$name <- colnames(x)[i]
    }

    return(structure(marginal, class = "mvd.marginal"))
}

#' Marginal function for use in fit.vine.copula: mixture of parametric distributions
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
fit.marginal.mixture <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x)))) {
    library(mixtools)
    library(betareg)

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    marginal <- list()
    marginal$type <- "mixture"
    marginal$dists <- list()
    for (i in 1:ncol(x)) {
        if (bounds[i, 1] >= bounds[i, 2]) {
            stop("Upper bound must be strictly higher than lower bound")
        }

        if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
            # [0,1] -> beta
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "beta"
            marginal$dists[[i]]$min <- 0
            marginal$dists[[i]]$max <- 1

            d <- data.frame(y = x[, i])
            cf <- coef(lm(qlogis(y) ~ 1, data = d))
            m <- betamix(y ~ 1 | 1, data = d, k = 1:3, control = betareg.control(start = list(mean = 0, precision = 1)), FLXcontrol = list(tolerance = 0.01))
            k <- m$flexmix@k
            if (k == 1) {
                mu <- plogis(coef(m)[1])
                phi <- exp(coef(m)[2])
                marginal$dists[[i]]$p <- 1
            } else {
                mu <- plogis(coef(m)[, 1])
                phi <- exp(coef(m)[, 2])
                marginal$dists[[i]]$p <- m$flexmix@prior
            }
            marginal$dists[[i]]$a <- mu * phi
            marginal$dists[[i]]$b <- (1 - mu) * phi
        } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
            # [0,inf] -> gamma
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "gamma"

            rescale <- 1 / mean(x[, i])
            rescaled <- x[, i] * rescale
            ms <- list()
            BIC <- rep(NA, 3)
            for (k in 1:3) {
                ms[[k]] <- gammamixEM(rescaled, k = k, epsilon = 0.01)
                BIC[k] <- log(nrow(x)) * 2 * k - 2 * ms[[k]]$loglik
            }
            k <- which.min(BIC)
            m <- ms[[k]]

            marginal$dists[[i]]$p <- m$lambda
            marginal$dists[[i]]$shape <- m$gamma.pars[1,]
            marginal$dists[[i]]$scale <- m$gamma.pars[2,] / rescale
        } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
            # [-inf,inf] -> normal
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "normal"

            ms <- list()
            BIC <- rep(NA, 3)
            ms[[1]] <- list()
            ms[[1]]$lambda <- 1
            ms[[1]]$mu <- mean(x[, i])
            ms[[1]]$sigma <- sd(x[, i])
            BIC[1] <- log(nrow(x)) * 2 - 2 * sum(dnorm(x[, i], ms[[1]]$mu, ms[[1]]$sigma, log = T))
            for (k in 2:3) {
                capture.output(return_value <- try(ms[[k]] <- normalmixEM(x[, i], k = k, epsilon = 0.01, arbmean = T, arbvar = T)))
                if (inherits(return_value, "try-error")) {
                    BIC[k] <- NA
                } else {
                    BIC[k] <- log(nrow(x)) * 2 * k - 2 * ms[[k]]$loglik
                }
            }
            k <- which.min(BIC)
            m <- ms[[k]]

            marginal$dists[[i]]$p <- m$lambda
            marginal$dists[[i]]$mu <- m$mu
            marginal$dists[[i]]$sigma <- m$sigma
        } else {
            # [a,b] -> scaled beta
            a <- bounds[i, 1];
            b <- bounds[i, 2];
            marginal$dists[[i]] <- list()
            marginal$dists[[i]]$type <- "beta"
            marginal$dists[[i]]$min <- a
            marginal$dists[[i]]$max <- b

            d <- data.frame(y = (x[, i] - a) / (b - a))
            m <- betamix(y ~ 1, data = d, k = 1:3, FLXcontrol = list(tolerance = 0.01))

            k <- m$flexmix@k
            if (k == 1) {
                mu <- plogis(coef(m)[1])
                phi <- exp(coef(m)[2])
                marginal$dists[[i]]$p <- 1
            } else {
                mu <- plogis(coef(m)[, 1])
                phi <- exp(coef(m)[, 2])
                marginal$dists[[i]]$p <- m$flexmix@prior
            }
            marginal$dists[[i]]$a <- mu * phi
            marginal$dists[[i]]$b <- (1 - mu) * phi
        }

        marginal$dists[[i]]$name <- colnames(x)[i]
    }

    return(structure(marginal, class = "mvd.marginal"))
}

.fit.pareto.tail <- function(x, u, fixed_beta = NA) {
    exceedances <- x[x > u]
    excess <- exceedances - u
    Nu <- length(excess)
    xbar <- mean(excess)

    result <- list()
    if (is.na(fixed_beta)) {
        s2 <- var(excess)
        xi0 <- -0.5 * (((xbar * xbar) / s2) - 1)
        beta0 <- 0.5 * xbar * (((xbar * xbar) / s2) + 1)
        theta <- c(xi0, beta0)

        negloglik <- function(theta, tmp) {
            xi <- theta[1]
            beta <- theta[2]
            cond1 <- beta <= 0
            cond2 <- (xi <= 0) && (max(tmp) > (-beta / xi))
            if (cond1 || cond2) {
                f <- 1e+06
            } else {
                y <- logb(1 + (xi * tmp) / beta)
                y <- y / xi
                f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
            }
            f
        }

        fit <- optim(theta, negloglik, lower = c(0, 0), upper = c(Inf, Inf), tmp = excess, method = "L-BFGS-B")
        result$xi <- fit$par[1]
        result$beta <- fit$par[2]
    } else {
        s2 <- var(excess)
        xi0 <- -0.5 * (((xbar * xbar) / s2) - 1)

        negloglik <- function(xi, tmp, beta) {
            cond1 <- beta <= 0
            cond2 <- (xi <= 0) && (max(tmp) > (-beta / xi))
            if (cond1 || cond2) {
                f <- 1e+06
            } else {
                y <- logb(1 + (xi * tmp) / beta)
                y <- y / xi
                f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
            }
            f
        }

        fit <- optimize(negloglik, lower = 0, upper = 1000, tmp = excess, beta = fixed_beta)
        result$xi <- fit$minimum
        result$beta <- fixed_beta
    }

    return(result)
}

.pareto.cdf <- function(x, u, xi, beta) {
    z <- (x - u) / beta
    return(1 - (1 + xi * z) ^ (-1/xi))
}

.pareto.pdf <- function(x, u, xi, beta, log = T) {
    z <- (x - u) / beta
    if (log) {
        p <- log(1 / beta) - (1 / xi + 1) * log(1 + xi * z)
    } else {
        p <- (1 / beta) * (1 + xi * z) ^ (-1 / xi - 1)
    }
    return(p)
}

#' Marginal function for use in fit.vine.copula: empirical distribution function with pareto tails
.pareto.quantile <- function(x, u, xi, beta) {
    return(u + beta * ((1.0 - x) ^ (-xi) - 1.0) / xi);
}

#' Marginal ecdf+pareto tail
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @export
fit.marginal.ecdf.pareto <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), reflect_bounds = T, pareto_threshold = 0.1, ecdf_bounds_threshold = 1e-3) {
    marginal <- list()
    marginal$type <- "ecdf.pareto"

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    marginal$ecdf <- fit.marginal.ecdf(x, bounds, reflect_bounds)

    use.pareto.tail <- matrix(T, ncol(x), 2)
    for (i in 1:ncol(x)) {
        if (bounds[i, 1] != -Inf) {
            dp <- sum(dnorm(bounds[i, 1], marginal$ecdf$x_reflected[[i]], marginal$ecdf$bw[i])) * marginal$ecdf$correction_factor[i]
            dp <- dp / length(marginal$ecdf$x_reflected[[i]])
            use.pareto.tail[i, 1] <- dp < ecdf_bounds_threshold
        }
        if (bounds[i, 2] != Inf) {
            dp <- sum(dnorm(bounds[i, 2], marginal$ecdf$x_reflected[[i]], marginal$ecdf$bw[i])) * marginal$ecdf$correction_factor[i]
            dp <- dp / length(marginal$ecdf$x_reflected[[i]])
            use.pareto.tail[i, 2] <- dp < ecdf_bounds_threshold
        }
    }

    for (i in 1:ncol(x)) {
        if (use.pareto.tail[i, 1]) {
            u <- as.numeric(quantile(x[, i], pareto_threshold))

            dp <- sum(dnorm(u, marginal$ecdf$x_reflected[[i]], marginal$ecdf$bw[i])) * marginal$ecdf$correction_factor[i]
            dp <- dp / length(marginal$ecdf$x_reflected[[i]])

            marginal$lower.tails[[i]] <- .fit.pareto.tail(-x[, i], - u, pareto_threshold / dp)
            marginal$lower.tails[[i]]$u <- u
            marginal$lower.tails[[i]]$q <- pareto_threshold
            marginal$lower.tails[[i]]$d <- dp
        } else {
            marginal$lower.tails[[i]] <- list()
        }
        if (use.pareto.tail[i, 2]) {
            u <- as.numeric(quantile(x[, i], 1 - pareto_threshold))

            dp <- sum(dnorm(u, marginal$ecdf$x_reflected[[i]], marginal$ecdf$bw[i])) * marginal$ecdf$correction_factor[i]
            dp <- dp / length(marginal$ecdf$x_reflected[[i]])

            marginal$upper.tails[[i]] <- .fit.pareto.tail(x[, i], u, pareto_threshold / dp)
            marginal$upper.tails[[i]]$u <- u
            marginal$upper.tails[[i]]$q <- pareto_threshold
            marginal$upper.tails[[i]]$d <- dp
        } else {
            marginal$upper.tails[[i]] <- list()
        }
    }

    return(structure(marginal, class = "mvd.marginal"))
}

#' Transform variables to U[0,1] using a fitted marginal function
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param marginal mvd.marginal object obtained from one of the marginal fitting functions.
#' @export
mvd.marginal.transform <- function(x, marginal) {
    stopifnot(class(marginal) == "mvd.marginal")

    if (!is.matrix(x)) {
        x <- t(as.matrix(x))
    }
    transformed <- matrix(nrow = nrow(x), ncol = ncol(x))

    if (marginal$type == "ecdf") {
        for (i in 1:ncol(x)) {
            transformed[, i] <- marginal$ecdfs[[i]](x[, i])
        }
    } else if (marginal$type == "ecdf.pareto") {
        for (i in 1:ncol(x)) {
            is_body <- rep(T, nrow(x))
            if (length(marginal$lower.tails[[i]]) > 0) {
                is_tail <- x[, i] < marginal$lower.tails[[i]]$u
                transformed[is_tail, i] <- marginal$lower.tails[[i]]$q - marginal$lower.tails[[i]]$q * .pareto.cdf(-x[is_tail, i], - marginal$lower.tails[[i]]$u, marginal$lower.tails[[i]]$xi, marginal$lower.tails[[i]]$beta)
                is_body[is_tail] <- F
            }
            if (length(marginal$upper.tails[[i]]) > 0) {
                is_tail <- x[, i] > marginal$upper.tails[[i]]$u
                transformed[is_tail, i] <- (1 - marginal$upper.tails[[i]]$q) + marginal$upper.tails[[i]]$q * .pareto.cdf(x[is_tail, i], marginal$upper.tails[[i]]$u, marginal$upper.tails[[i]]$xi, marginal$upper.tails[[i]]$beta)
                is_body[is_tail] <- F
            }
            transformed[is_body, i] <- marginal$ecdf$ecdfs[[i]](x[is_body, i])
        }
    } else if (marginal$type == "parametric") {
        for (i in 1:ncol(x)) {
            if (marginal$dists[[i]]$type == "beta") {
                transformed[, i] <- pbeta((x[, i] - marginal$dists[[i]]$min) / (marginal$dists[[i]]$max - marginal$dists[[i]]$min), marginal$dists[[i]]$shape1, marginal$dists[[i]]$shape2)
            } else if (marginal$dists[[i]]$type == "normal") {
                transformed[, i] <- pnorm(x[, i], marginal$dists[[i]]$mu, marginal$dists[[i]]$sigma)
            } else if (marginal$dists[[i]]$type == "gamma") {
                transformed[, i] <- pgamma(x[, i], shape = marginal$dists[[i]]$shape, scale = marginal$dists[[i]]$scale)
            }
        }
    } else if (marginal$type == "mixture") {
        transformed <- matrix(0, nrow(x), ncol(x))
        for (i in 1:ncol(x)) {
            margin <- marginal$dists[[i]]
            if (margin$type == "beta") {
                for (j in 1:length(margin$p)) {
                    transformed[, i] <- transformed[, i] + margin$p[j] * pbeta((x[, i] - margin$min) / (margin$max - margin$min), margin$a[j], margin$b[j])
                }
            } else if (marginal$dists[[i]]$type == "normal") {
                for (j in 1:length(margin$p)) {
                    transformed[, i] <- transformed[, i] + margin$p[j] * pnorm(x[, i], margin$mu[j], margin$sigma[j])
                }
            } else if (marginal$dists[[i]]$type == "gamma") {
                for (j in 1:length(margin$p)) {
                    transformed[, i] <- transformed[, i] + margin$p[j] * pgamma(x[, i], shape = margin$shape[j], scale = margin$scale[j])
                }
            }
        }
    } else {
        stop("Unknown marginal type")
    }

    return(transformed)
}

#' Reverse transform variables from U[0,1]
#'
#' description
#' @param transformed Matrix or vector of transformed samples. For matrices, rows are samples and columns are variables.
#' @export
reverse.transform.marginals <- function(transformed, marginal) {
    stopifnot(class(marginal) == "mvd.marginal")

    if (!is.matrix(transformed)) {
        transformed <- as.matrix(transformed)
    }
    x <- matrix(nrow = nrow(transformed), ncol = ncol(transformed))

    if (marginal$type == "ecdf") {
    } else if (marginal$type == "ecdf.pareto") {
        for (i in 1:ncol(transformed)) {
            is_body <- rep(T, nrow(x))
            if (length(marginal$lower.tails[[i]]) > 0) {
                tail <- marginal$lower.tails[[i]]
                is_tail <- transformed[, i] < tail$q
                x[is_tail, i] <-  -.pareto.quantile(1.0 - (transformed[is_tail, i] / tail$q), - tail$u, tail$xi, tail$beta);
                is_body[is_tail] <- F
            }
            if (length(marginal$upper.tails[[i]]) > 0) {
                tail <- marginal$upper.tails[[i]]
                is_tail <- transformed[, i] > 1.0 - tail$q
                x[is_tail, i] <- .pareto.quantile((transformed[is_tail, i] - (1.0 - tail$q)) / tail$q, tail$u, tail$xi, tail$beta);
                is_body[is_tail] <- F
            }
            x[is_body, i] <- quantile(marginal$ecdf$ecdfs[[i]], transformed[is_body, i], type = 4)
        }
    } else if (marginal$type == "parametric") {
    } else if (marginal$type == "mixture") {
        for (i in 1:ncol(transformed)) {
            margin <- marginal$dists[[i]]
            ix <- sample(length(margin$p), nrow(x), replace = T, prob = margin$p)
            if (margin$type == "beta") {
                x[, i] <- margin$min + qbeta(transformed[, i], margin$a[ix], margin$b[ix]) * (margin$max - margin$min)
            } else if (margin$type == "beta") {
                x[, i] <- qnorm(transformed[, i], margin$mu[ix], margin$sigma[ix])
            } else if (margin$type == "beta") {
                x[, i] <- qbeta(transformed[, i], shape = margin$shape[ix], scale = margin$scale[ix])
            }
        }
    } else {
        stop("Unknown marginal type")
    }

    return(x)
}

#' Probability density function of a marginal distribution fitted with one of the marginal distribution functions.
#'
#' description
#' @param marginal mvd.marginal object obtained from one of the marginal fitting functions.
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param log Whether to return the density in log
#' @export
marginal.pdf <- function(marginal, x, log = T) {
    stopifnot(class(marginal) == "mvd.marginal")
    stopifnot(log)

    p <- matrix(NA, nrow(x), ncol(x))

    if (marginal$type == "ecdf") {
        stopifnot(ncol(x) == length(marginal$bw))
        for (i in 1:ncol(x)) {
            for (j in 1:nrow(x)) {
                dp <- sum(dnorm(x[j, i], marginal$x[[i]], marginal$bw[i])) * marginal$correction_factor[i]
                p[j, i] <- log(dp) - log(length(marginal$x[[i]]))
            }
        }
    } else if (marginal$type == "ecdf.pareto") {
        for (i in 1:ncol(x)) {
            is_body <- rep(T, nrow(x))
            if (length(marginal$lower.tails[[i]]) > 0) {
                is_tail <- x[, i] < marginal$lower.tails[[i]]$u
                p[is_tail, i] <- log(marginal$lower.tails[[i]]$q) + .pareto.pdf(-x[is_tail, i], - marginal$lower.tails[[i]]$u, marginal$lower.tails[[i]]$xi, marginal$lower.tails[[i]]$beta, log = T)
                is_body[is_tail] <- F
            }
            if (length(marginal$upper.tails[[i]]) > 0) {
                is_tail <- x[, i] > marginal$upper.tails[[i]]$u
                p[is_tail, i] <- log(marginal$upper.tails[[i]]$q) + .pareto.pdf(x[is_tail, i], marginal$upper.tails[[i]]$u, marginal$upper.tails[[i]]$xi, marginal$upper.tails[[i]]$beta, log = T)
                is_body[is_tail] <- F
            }
            for (j in which(is_body)) {
                dp <- sum(dnorm(x[j, i], marginal$ecdf$x[[i]], marginal$ecdf$bw[i])) * marginal$ecdf$correction_factor[i]
                p[j, i] <- log(dp) - log(length(marginal$ecdf$x[[i]]))
            }
        }
    } else if (marginal$type == "parametric") {
        for (i in 1:ncol(x)) {
            if (marginal$dists[[i]]$type == "beta") {
                a <- marginal$dists[[i]]$min
                b <- marginal$dists[[i]]$max
                p[, i] <- dbeta((x[, i] - a) / (b - a), marginal$dists[[i]]$shape1, marginal$dists[[i]]$shape2, log = T) - log(b - a)
            } else if (marginal$dists[[i]]$type == "normal") {
                p[, i] <- dnorm(x[, i], marginal$dists[[i]]$mu, marginal$dists[[i]]$sd, log = T)
            } else if (marginal$dists[[i]]$type == "gamma") {
                p[, i] <- dgamma(x[, i], shape = marginal$dists[[i]]$shape, scale = marginal$dists[[i]]$scale, log = T)
            }
        }
    } else if (marginal$type == "mixture") {
        for (i in 1:ncol(x)) {
            margin <- marginal$dists[[i]]
            dp <- rep(0, nrow(x))
            if (margin$type == "beta") {
                for (j in 1:length(margin$p)) {
                    dp <- dp + margin$p[j] * dbeta((x[, i] - margin$min) / (margin$max - margin$min), margin$a[j], margin$b[j])
                }
                p[, i] <- log(dp) - log(margin$max - margin$min)
            } else if (margin$type == "normal") {
                for (j in 1:length(margin$p)) {
                    dp <- dp + margin$p[j] * dnorm(x[, i], margin$mu[j], margin$sigma[j])
                }
                p[, i] <- log(dp)
            } else if (margin$type == "gamma") {
                for (j in 1:length(margin$p)) {
                    dp <- dp + margin$p[j] * dgamma(x[, i], shape = margin$shape[j], scale = margin$scale[j])
                }
                p[, i] <- log(dp)
            }
        }
    } else {
        stop("Unknown marginal type")
    }
    return(p)
}

marginal.correct.p <- function(marginal, x, p, log = T) {
    stopifnot(class(marginal) == "mvd.marginal")
    mp <- marginal.pdf(marginal, x, log)
    for (i in 1:ncol(x)) {
        p <- p + mp[, i]
    }
    return(p)
}
