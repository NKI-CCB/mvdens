#' Evaluate the probability density function of a density approximation
#'
#' description
#' @param fit An mvd.density object obtained from one of the density fitting functions.
#' @param x Matrix or vector of positions at which to evaluate the density function.
#' @param log Boolean which specifies whether to return the probability density in log space.
#' @export
mvd.pdf <- function(fit, x, log = FALSE) {
    stopifnot(class(fit) == "mvd.density")

    if (fit$type == "kde") {
        return(.evaluate.kde(fit, x, log))
    } else if (fit$type == "kde.transformed") {
        transformed <- mvd.transform_to_unbounded(x, fit$transform.bounds)
        p <- .evaluate.kde(fit$kde, transformed, log = log)
        return(mvd.correct_p_for_transformation(x, fit$transform.bounds, p, log = log))
    } else if (fit$type == "gmm") {
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(fit$centers))

            p <- matrix(NA, nrow(x), fit$num_components)
            for (ki in 1:fit$num_components) {
                if (log) {
                    p[, ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[, ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(fit$centers))

            p <- rep(NA, fit$num_components)
            for (ki in 1:fit$num_components) {
                if (log) {
                    p[ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (fit$type == "gmm.transformed") {
        transformed <- mvd.transform_to_unbounded(x, fit$transform.bounds)
        p <- mvd.pdf(fit$gmm, transformed, log = log)
        return(mvd.correct_p_for_transformation(x, fit$transform.bounds, p, log = log))
    } else if (fit$type == "gmm.truncated") {
        if (is.matrix(x)) {
            stopifnot(ncol(x) == ncol(fit$centers))

            p <- matrix(NA, nrow(x), fit$num_components)
            for (ki in 1:fit$num_components) {
                if (log) {
                    p[, ki] <- tmvtnorm::dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[, ki] <- tmvtnorm::dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(apply(p, 1, .logsum))
            } else {
                return(apply(p, 1, sum))
            }
        } else {
            stopifnot(length(x) == ncol(fit$centers))

            p <- rep(NA, fit$num_components)
            for (ki in 1:fit$num_components) {
                if (log) {
                    p[ki] <- tmvtnorm::dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = TRUE) + log(fit$proportions[ki])
                } else {
                    p[ki] <- tmvtnorm::dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = FALSE) * fit$proportions[ki]
                }
            }

            if (log) {
                return(.logsum(p))
            } else {
                return(sum(p))
            }
        }
    } else if (fit$type == "gp") {
        return(evaluate.gp(fit, x, log))
    } else if (fit$type == "vine.copula") {
        return(evaluate.vine.copula(fit, x, log))
    } else if (fit$type == "mfa") {
      return(evaluate.factor.mixture(fit, x, log))
    } else if (fit$type == "mfa.transformed") {
      transformed <- mvd.transform_to_unbounded(x, fit$transform.bounds)
      p <- mvd.pdf(fit$mfa, transformed, log = log)
      return(mvd.correct_p_for_transformation(x, fit$transform.bounds, p, log = log))
    } else {
        stop("Unknown type")
    }
}
