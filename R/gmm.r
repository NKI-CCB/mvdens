library(mvtnorm)
library(tmvtnorm)

#source("utils.r")

.assign_kmeanspp <- function(x, K) {
    stopifnot(K >= 1)
    stopifnot(is.matrix(x))

    result <- list()
    result$k <- K

    center_ix <- rep(NA, K)
    center_ix[1] <- sample.int(nrow(x), 1)
    result$centers <- matrix(NA, K, ncol(x))
    result$centers[1,] <- x[center_ix[1],]

    if (K > 1) {
        for (i in 2:K) {
            min_dist_sq <- rep(NA, nrow(x))
            for (xi in 1:nrow(x)) {
                dist_sq <- rep(NA, i - 1)
                for (j in 1:(i - 1)) {
                    dist_sq[j] <- sum((x[xi,] - result$centers[j,]) ^ 2)
                }
                min_dist_sq[xi] <- min(dist_sq)
            }

            total <- sum(min_dist_sq)
            center_ix[i] <- sample.int(nrow(x), 1, prob = min_dist_sq / total)
            result$centers[i,] <- x[center_ix[i],]
        }
    }

    result$weights <- matrix(0, nrow(x), K)

    for (xi in 1:nrow(x)) {
        dist_sq <- apply(result$centers, 1, function(v) { return(sum((x[xi,] - v) ^ 2)) })
        assign_ix <- which.min(dist_sq)
        result$weights[xi, assign_ix] <- 1
    }

    return(result)
}

.mixture_expectation_step <- function(x, fit, truncated=T) {
    log_l <- matrix(NA, nrow(x), fit$k)
    for (ki in 1:fit$k) {
        if (truncated) {
            log_l[, ki] <- dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]],
                                lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = T) + log(fit$component_weights[ki])
        } else {
            log_l[, ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = T) + log(fit$component_weights[ki])
        }
    }
    fit$logl <- sum(apply(log_l, 1, .logsum))

    fit$weights <- exp(log_l)
    total_weights <- apply(fit$weights, 1, sum)
    fit$weights <- fit$weights / total_weights

    return(fit)
}

.mixture_maximization_step <- function(x, fit) {
    fit$component_weights <- rep(NA, fit$k)
    for (ki in 1:fit$k) {
        fit$component_weights[ki] <- sum(fit$weights[, ki]) / nrow(x)
        weighted <- cov.wt(x, fit$weights[, ki])
        fit$centers[ki,] <- weighted$center
        fit$covariances[[ki]] <- weighted$cov
    }
    return(fit)
}

.fit.gmm.internal <- function(x, K, truncated, bounds, epsilon, maxsteps, verbose) {
    singular <- F
    for (j in 1:10) {
        fit <- .assign_kmeanspp(x, K)
        fit$bounds <- bounds
        fit$covariances <- list()
        fit <- .mixture_maximization_step(x, fit)
        
        previous_logl <- -Inf
        singular <- F
        for (i in 1:maxsteps) {
            return_value <- try(fit <- .mixture_expectation_step(x, fit, truncated), silent = !verbose)
            if (inherits(return_value, "try-error")) {
                singular <- T
                break
            }

            if (any(is.na(fit$weights))) {
                singular <- T
                break
            }
            if (any(apply(fit$weights, 2, sum, na.rm=T) == 0)) {
                singular <- T
                break
            }

            if (verbose) {
                cat("Log likelihood: ", fit$logl, "\n")
            }
            if (fit$logl - previous_logl < epsilon) {
                return(fit)
            } else {
                previous_logl <- fit$logl
            }

            fit <- .mixture_maximization_step(x, fit)
        }

        if (!singular) {
            break
        }
    }

    if (singular) {
        fit$logl <- NA
        warning("singular solutions for 10 restarts")
    } else {
        warning("not converged after maxsteps")
    }
    return(fit)
}

fit.gmm <- function(x, K, epsilon = 0.01, maxsteps = 100, verbose = F) {
    fit <- .fit.gmm.internal(x, K, truncated = F, bounds = NA, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)

    result <- list()
    result$type <- "gmm"
    result$K <- K
    result$centers <- fit$centers
    result$covariances <- fit$covariances
    result$proportions <- fit$component_weights
    result$log_likelihood <- fit$logl
    nparam <- K * ncol(x) * (ncol(x) + 1) / 2 + K
    result$BIC <- log(nrow(x)) * nparam - 2 * fit$logl
    result$assignment <- apply(fit$weights, 1, which.max)
    return(structure(result, class = "mdd.density"))
}

fit.gmm.transformed <- function(x, K, bounds, epsilon = 0.01, maxsteps = 100, verbose = F) {
    result <- list()
    result$type <- "gmm.transformed"
    result$transform.bounds <- bounds
    transformed <- mdd.transform_to_unbounded(x, bounds)
    result$gmm <- fit.gmm(transformed, K, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
    return(structure(result, class = "mdd.density"))
}

fit.gmm.truncated <- function(x, K, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), epsilon = 0.01, maxsteps = 100, verbose = F) {
    fit <- .fit.gmm.internal(x, K, truncated = T, bounds = bounds, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)

    result <- list()
    result$type <- "gmm.truncated"
    result$K <- K
    result$centers <- fit$centers
    result$covariances <- fit$covariances
    result$proportions <- fit$component_weights
    result$bounds <- bounds
    result$log_likelihood <- fit$logl
    nparam <- K * ncol(x) * (ncol(x) + 1) / 2 + K
    result$BIC <- log(nrow(x)) * nparam - 2 * fit$logl
    result$assignment <- apply(fit$weights, 1, which.max)
    return(structure(result, class = "mdd.density"))
}

gmm.BIC <- function(x, K = 1:9, epsilon = 0.01, maxsteps = 100, verbose = F) {
    result <- list()
    result$K <- K
    result$BIC <- rep(NA, length(K))
    result$fits <- list()
    for (k in 1:length(K)) {
        if (verbose) {
            cat("Fitting K=", K[k], "\n")
        }
        result$fits[[k]] <- fit.gmm(x, K[k], epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
        result$BIC[k] <- result$fits[[k]]$BIC
    }
    return(result)
}

gmm.transformed.BIC <- function(x, K = 1:9, bounds, epsilon = 0.01, maxsteps = 100, verbose = F) {
  result <- list()
  result$K <- K
  result$BIC <- rep(NA, length(K))
  result$fits <- list()
  for (k in 1:length(K)) {
    if (verbose) {
      cat("Fitting K=", K[k], "\n")
    }
    result$fits[[k]] <- fit.gmm.transformed(x, K[k], bounds, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
    result$BIC[k] <- result$fits[[k]]$gmm$BIC
  }
  return(result)
}

gmm.truncated.BIC <- function(x, K = 1:9, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), epsilon = 0.01, maxsteps = 100, verbose = F) {
    result <- list()
    result$K <- K
    result$BIC <- rep(NA, length(K))
    result$fits <- list()
    for (k in 1:length(K)) {
        if (verbose) {
            cat("Fitting K=", K[k], "\n")
        }
        result$fits[[k]] <- fit.gmm.truncated(x, K[k], bounds = bounds, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
        result$BIC[k] <- result$fits[[k]]$BIC
    }
    return(result)
}