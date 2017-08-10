.logit <- function(x) { log(x / (1 - x)) }
.dlogit <- function(x) { 1 / (x - x * x) }
.logistic <- function(x) { 1 / (1 + exp(-x)) }
.logit_scale <- function(x, a, b) { log((a - x) / (x - b)) }
.dlogit_scale <- function(x, a, b) {(b - a) / ((a - x) * (x - b)) }
.logistic_scale <- function(x, a, b) { ex <- exp(x); (a + b * ex) / (ex + 1) }

mdd.transform_to_unbounded <- function(x, bounds) {
    transformed <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (i in 1:ncol(x)) {
        if (bounds[i, 1] >= bounds[i, 2]) {
            stop("Upper bound must be strictly higher than lower bound")
        }

        if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
            # [0,1] -> logit
            transformed[, i] <- .logit(x[, i])
        } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
            # [0,inf] -> log
            transformed[, i] <- log(x[, i])
        } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
            # [-inf,inf] -> no transform
            transformed[, i] <- x[, i]
        } else {
            # [a,b] -> scaled logit
            transformed[, i] <- .logit_scale(x[, i], bounds[i, 1], bounds[i, 2])
        }
    }
    return(transformed)
}

mdd.transform_from_unbounded <- function(transformed_x, bounds) {
    x <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (i in 1:ncol(x)) {
        if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
            # [0,1] -> logit -> inverse = logistic
            x[, i] <- .logistic(transformed_x[, i])
        } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
            # [0,inf] -> log -> inverse = exp
            x[, i] <- exp(transformed_x[, i])
        } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
            # [-inf,inf] -> no transform
            x[, i] <- transformed_x[, i]
        } else {
            # [a,b] -> scaled logit -> inverse = scaled logistic
            x[, i] <- .logistic_scale(transformed_x[, i], bounds[i, 1], bounds[i, 2])
        }
    }
    return(x)
}

mdd.correct_p_for_transformation <- function(x, bounds, p, log = T) {
    if (log) {
        for (i in 1:ncol(x)) {
            if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
                # [0,1] -> logit -> derivative = dlogit
                p <- p + log(.dlogit(x[, i]))
            } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
                # [0,inf] -> log -> derivative = 1/log
                p <- p - log(x[, i])
            } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
                # [-inf,inf] -> no transform
            } else {
                # [a,b] -> scaled logit -> derivative = dlogit_scale
                p <- p + log(.dlogit_scale(x[, i], bounds[i, 1], bounds[i, 2]))
            }
        }
    } else {
        for (i in 1:ncol(x)) {
            if (bounds[i, 1] == 0 && bounds[i, 2] == 1) {
                # [0,1] -> logit -> derivative = dlogit
                p <- p * .dlogit(x[, i])
            } else if (bounds[i, 1] == 0 && bounds[i, 2] == Inf) {
                # [0,inf] ->log(x) -> 1/x
                p <- p / x[, i]
            } else if (bounds[i, 1] == -Inf && bounds[i, 2] == Inf) {
                # [-inf,inf] -> no transform
            } else {
                # [a,b] -> scaled logit -> derivative = dlogit_scale
                p <- p * .dlogit_scale(x[, i], bounds[i, 1], bounds[i, 2])
            }
        }
    }
    return(p)
}
