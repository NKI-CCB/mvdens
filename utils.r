.logsum <- function(x) {
    xmax <- which.max(x)
    log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
}