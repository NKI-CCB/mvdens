library(mvtnorm)
library(tmvtnorm)
library(VineCopula)
library(mixtools)
library(ks)

.logsum <- function(x) {
    xmax <- which.max(x)
    log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
}

setClass('mvd.marginal')
setMethod('summary', 'mvd.marginal', function(object, ...) {
    cat("mvdens marginal distribution of type:", object$type, "\n")
    if (object$type == "ecdf") {
        cat("  Bandwidths:\n")
        for (i in 1:length(object$bw)) {
            cat("   ", object$varnames[i], "-", object$bw[i], "\n")
        }
    } else if (object$type == "parametric") {
        for (i in 1:length(object$dists)) {
            cat("   ", object$dists[[i]]$name, "-", object$dists[[i]]$type, "\n")
        }
    } else if (object$type == "mixture") {
        for (i in 1:length(object$dists)) {
            cat("   ", object$dists[[i]]$name, "-", object$dists[[i]]$type, "- p:", paste(object$dists[[i]]$p, collapse = ","), "\n")
        }
    }
})

setClass('mvd.density')
setMethod('summary', 'mvd.density', function(object, ...) {
    cat("mvdens density approximation of type:", object$type, "\n")
    if (object$type == "kde") {
        cat("  Bandwidth matrix:")
        str(object$H)
        cat("\n")
    } else if (object$type == "kde.transformed") {
        cat("  Transformations:\n")
        for (i in 1:length(object$transform.bounds)) {
            cat("   ", i, "-", mvd.transform_name(object$transform.bounds, i), "\n")
        }
        cat("  Bandwidth matrix:")
        str(object$H)
        cat("\n")
    } else if (object$type == "gmm") {
    } else if (object$type == "gmm.transformed") {
    } else if (object$type == "gmm.truncated") {
    } else if (object$type == "gp") {
        cat("  kernel=", object$kernel, "l = ", object$l, ", s := ", object$s, " \n")
    } else if (object$type == "vine.copula") {
    }
})
