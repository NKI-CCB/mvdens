#source('pdf.r')

mdd.marginal.likelihood <- function(mddens_fit, x, logp)
{
    fitp <- mdd.pdf(mddens_fit, x, log = TRUE)
    res <- lm(logp - fitp ~ 1)
    return(as.numeric(res$coefficients["(Intercept)"]))
}
