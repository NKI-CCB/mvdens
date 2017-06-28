

fit.marginal.ecdf <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), reflect_bounds = T) {
    marginal <- list()
    marginal$type <- "ecdf"
    marginal$ecdfs <- apply(x, 2, ecdf)

    marginal$varnames <- colnames(x)
    marginal$x <- list()
    marginal$correction_factor <- rep(1, ncol(x))
    for (i in 1:ncol(x)) {
        marginal$x[[i]] <- x[,i]
        if (reflect_bounds) {
            if (bounds[i,1] != -Inf) {
                reflected <- bounds[i, 1] - (x[, i] - bounds[i, 1])
                marginal$x[[i]] <- c(marginal$x[[i]], reflected)
                marginal$correction_factor[i] <- marginal$correction_factor[i] + 1
            }
            if (bounds[i, 2] != Inf) {
                reflected <- bounds[i, 2] + (bounds[i, 2] - x[, i])
                marginal$x[[i]] <- c(marginal$x[[i]], reflected)
                marginal$correction_factor[i] <- marginal$correction_factor[i] + 1
            }
        }
    }

    marginal$bw <- rep(NA, ncol(x))
    for (i in 1:ncol(x)) {
        hmax <- 10 * sd(marginal$x[[i]])
        marginal$bw[i] <- bw.SJ(marginal$x[[i]], lower = 1e-6 * hmax, upper = hmax)
    }

    return(structure(marginal, class="mdd.marginal"))
}

fit.marginal.parametric <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x)))) {
    library(MASS)

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
            marginal$dists[[i]]$mean <- mean(x[, i])
            marginal$dists[[i]]$sd <- sd(x[, i])
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

    return(structure(marginal, class = "mdd.marginal"))
}

fit.marginal.mixture <- function(x, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x)))) {
    library(mixtools)
    library(betareg)

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
            BIC[1] <- log(nrow(x)) * 2 - 2 * sum(dnorm(x[,i], ms[[1]]$mu, ms[[1]]$sigma, log = T))
            for (k in 2:3) {
                ms[[k]] <- normalmixEM(x[,i], k = k, epsilon = 0.01, arbmean = T, arbvar = T)
                BIC[k] <- log(nrow(x)) * 2 * k - 2 * ms[[k]]$loglik
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

    return(structure(marginal, class = "mdd.marginal"))
}

transform.marginals <- function(x, marginal) {
    transformed <- matrix(nrow = nrow(x), ncol = ncol(x))

    if (marginal$type == "ecdf") {
        for (i in 1:ncol(x)) {
            transformed[, i] <- 1e-10 + (1-2e-10) * marginal$ecdfs[[i]](x[, i])
        }
    } else if (marginal$type == "parametric") {
        for (i in 1:ncol(x)) {
            if (marginal$dists[[i]]$type == "beta") {
                transformed[, i] <- pbeta((x[, i] - marginal$dists[[i]]$min) / (marginal$dists[[i]]$max - marginal$dists[[i]]$min), marginal$dists[[i]]$shape1, marginal$dists[[i]]$shape2)
            } else if (marginal$dists[[i]]$type == "normal") {
                transformed[, i] <- pnorm(x[, i], marginal$dists[[i]]$mean, marginal$dists[[i]]$sd)
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
    }

    return(transformed)
}

marginal.correct.p <- function(marginal, x, p, log = T) {
    stopifnot(log)

    if (marginal$type == "ecdf") {
        stopifnot(ncol(x) == length(marginal$bw))
        for (i in 1:ncol(x)) {
            for (j in 1:nrow(x)) {
                dp <- sum(dnorm(x[j, i], marginal$x[[i]], marginal$bw[i])) * marginal$correction_factor[i]
                p[j] <- p[j] + log(dp) - log(length(marginal$x[[i]]))
            }
        }
    } else if (marginal$type == "parametric") {
        for (i in 1:ncol(x)) {
            if (marginal$dists[[i]]$type == "beta") {
                a <- marginal$dists[[i]]$min
                b <- marginal$dists[[i]]$max
                p <- p + dbeta((x[, i] - a) / (b - a), marginal$dists[[i]]$shape1, marginal$dists[[i]]$shape2, log = T) - log(b - a)
            } else if (marginal$dists[[i]]$type == "normal") {
                p <- p + dnorm(x[, i], marginal$dists[[i]]$mean, marginal$dists[[i]]$sd, log = T)
            } else if (marginal$dists[[i]]$type == "gamma") {
                p <- p + dgamma(x[, i], shape = marginal$dists[[i]]$shape, scale = marginal$dists[[i]]$scale, log = T)
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
                p <- p + log(dp) - log(margin$max - margin$min)
            } else if (margin$type == "normal") {
                for (j in 1:length(margin$p)) {
                    dp <- dp + margin$p[j] * dnorm(x[, i], margin$mu[j], margin$sigma[j])
                }
                p <- p + log(dp)
            } else if (margin$type == "gamma") {
                for (j in 1:length(margin$p)) {
                    dp <- dp + margin$p[j] * dgamma(x[, i], shape = margin$shape[j], scale = margin$scale[j])
                }
                p <- p + log(dp)
            }
        }
    } else {
        stop("Unknown marginal type")
    }

    return(p)
}

# setClass('mdd.marginal')
# setMethod(summary, 'mdd.marginal', function(object, ...) {
#     cat("mddens marginal distribution of type:", object$type, "\n")
#     if (object$type == "ecdf") {
#         cat("  Number of samples:", nrow(object$x), "\n")
#         cat("  Bandwidths:\n")
#         for (i in 1:length(object$bw)) {
#             cat("   ", object$varnames[i], "-", object$bw[i], "\n")
#         }
#     } else if (object$type == "parametric") {
#         for (i in 1:length(object$dists)) {
#             cat("   ", object$dists[[i]]$name, "-", object$dists[[i]]$type, "\n")
#         }
#     } else if (object$type == "mixture") {
#         for (i in 1:length(object$dists)) {
#             cat("   ", object$dists[[i]]$name, "-", object$dists[[i]]$type, "- p:", paste(object$dists[[i]]$p, collapse = ","), "\n")
#         }
#     }
# })
