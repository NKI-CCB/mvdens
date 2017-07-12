source("utils.r")
source("gp.r")
source("vine_copula.r")

mdd.pdf <- function(mddens_fit, x, log = FALSE) {
    if (mddens_fit$type == "kde") {
        return(evaluate.kde(mddens_fit, x, log))
    } else if (mddens_fit$type == "kde.transformed") {
        transformed <- mdd.transform_to_unbounded(x, mddens_fit$transform.bounds)
        p <- evaluate.kde(mddens_fit$kde, transformed, log = log)
        return(mdd.correct_p_for_transformation(x, mddens_fit$transform.bounds, p, log = log))
    } else if (mddens_fit$type == "gmm") {
        require(mvtnorm)

        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(mddens_fit$centers))

            p <- matrix(NA, nrow(x), mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[, ki] <- mvtnorm::dmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[, ki] <- mvtnorm::dmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = FALSE) * mddens_fit$proportions[ki]
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
                    p[ki] <- mvtnorm::dmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[ki] <- mvtnorm::dmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (mddens_fit$type == "gmm.transformed") {
        transformed <- mdd.transform_to_unbounded(x, mddens_fit$transform.bounds)
        p <- mdd.pdf(mddens_fit$gmm, transformed, log = log)
        return(mdd.correct_p_for_transformation(x, mddens_fit$transform.bounds, p, log = log))
    } else if (mddens_fit$type == "gmm.truncated") {
        require(tmvtnorm)
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(mddens_fit$centers))

            p <- matrix(NA, nrow(x), mddens_fit$K)
            for (ki in 1:mddens_fit$K) {
                if (log) {
                    p[, ki] <- dtmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[, ki] <- dtmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = FALSE) * mddens_fit$proportions[ki]
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
                    p[ki] <- dtmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = TRUE) + log(mddens_fit$proportions[ki])
                } else {
                    p[ki] <- dtmvnorm(x, mddens_fit$centers[ki,], mddens_fit$covariances[[ki]], lower = mddens_fit$bounds[, 1], upper = mddens_fit$bounds[, 2], log = FALSE) * mddens_fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (mddens_fit$type == "gp") {
        return(evaluate.gp(mddens_fit, x))
    } else if (mddens_fit$type == "vine.copula") {
        return(evaluate.vine.copula(mddens_fit, x))
    } else {
        stop("Unknown type")
    }
}
