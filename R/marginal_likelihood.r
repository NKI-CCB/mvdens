#source('pdf.r')

#' Calculate a marginal likelihood estimate from a density fit
#'
#' description
#' @param fit 
#' @param x
#' @param p
#' @param log
#' @export
#' @examples
mdd.marginal.likelihood <- function(fit, x, p, log)
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
        fitp <- mdd.pdf(fit, x, log = TRUE)
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
