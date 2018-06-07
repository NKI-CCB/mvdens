#' Calculate a marginal likelihood estimate from a density fit
#'
#' Estimates the marginal likelihood from a set of Monte Carlo samples using a fitted posterior density approximation.
#' For Gaussian process fits, the marginal likelihood estimation only makes sense if the GP was normalized during the fitting.
#' @param fit An mvd.density object obtained from one of the density fitting functions.
#' @param x Matrix or vector of positions at which to evaluate the density function.
#' @param p Unnormalized posterior probability values.
#' @param log Boolean which specifies whether p is in log scale. If true, the marginal likelihood is also return in log scale.
#' @export
mvd.marginal.likelihood <- function(fit, x, p, log)
{
    if (fit$type == "gp") {
        if (fit$log == T) {
            stop("Cannot calculate marginal likelihood estimate for Gaussian process on log scale")
        } else {
            if (log) {
                return(-log(fit$s))
            } else {
                return(1/fit$s)
            }
        }
    } else {
        fitp <- mvd.pdf(fit, x, log = TRUE)
        if (log) {
            res <- lm(p - fitp ~ 1)
            logml <- as.numeric(res$coefficients["(Intercept)"])
            return(logml)
        } else {
            res <- lm(log(p) - fitp ~ 1)
            logml <- as.numeric(res$coefficients["(Intercept)"])
            return(exp(logml))
        }
    }
}
