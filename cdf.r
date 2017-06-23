source("utils.r")

mdd.cdf <- function(mddens_fit, x, log) {
    if (mddens_fit$type == "gmm") {
        require(mvtnorm)

        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(mddens_fit$centers))

            p <- matrix(NA, nrow(x), mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[, ki] <- pmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[, ki] <- pmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(mddens_fit$centers))

            p <- rep(NA, mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[ki] <- pmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[ki] <- pmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (mddens_fit$type == "gmm.truncated") {
        require(tmvtnorm)
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(mddens_fit$centers))

            p <- matrix(NA, nrow(x), mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[, ki] <- ptmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[, ki] <- ptmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(mddens_fit$centers))

            p <- rep(NA, mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[ki] <- ptmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[ki] <- ptmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else {
        stop("Unknown type")
    }
}
