#' Evaluate the cumulative distribution function of a density approximation
#'
#' Currently only implemented for Gaussian mixtures.
#' @param fit An mvd.density object obtained from one of the density fitting functions.
#' @param x Matrix or vector of positions at which to evaluate the density function.
#' @param log Boolean which specifies whether to return the cumulative probability in log space.
#' @export
mvd.cdf <- function(fit, x, log) {
    stopifnot(class(fit) == "mvd.density")

    if (fit$type == "gmm") {
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(fit$centers))

            p <- matrix(NA, nrow(x), fit$K)
            for (ki in 1:fit$K) {
                if (log) {
                    p[, ki] <- pmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[, ki] <- pmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(fit$centers))

            p <- rep(NA, fit$K)
            for (ki in 1:fit$K) {
                if (log) {
                    p[ki] <- pmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[ki] <- pmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (fit$type == "gmm.truncated") {
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(fit$centers))

            p <- matrix(NA, nrow(x), fit$K)
            for (ki in 1:fit$K) {
                if (log) {
                    p[, ki] <- ptmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[, ki] <- ptmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(fit$centers))

            p <- rep(NA, fit$K)
            for (ki in 1:fit$K) {
                if (log) {
                    p[ki] <- ptmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[ki] <- ptmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else {
        stop("Unknown type or not implemented")
    }
}
