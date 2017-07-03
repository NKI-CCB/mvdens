source('pdf.r')

mdd.marginal.likelihood <- function(mddens_fit, x, logp)
{
    fitp <- mdd.pdf(mddens_fit, x, log = TRUE)
    res <- lm(logp - fitp ~ 1)
    return(as.numeric(res$coefficients["(Intercept)"]))
}

source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))

cwd <- getwd()
setwd("D:/Research/projects/52-approximate-mc-output-paper/")
model <- load_model("lotka-volterra/lotka_volterra_trapperonly", "output_pt_16_1000", "lotka_volterra_trapperonly.xml")
setwd(cwd)

bounds <- prior_bounds_all(model, q = c(0, 1))

source("gmm.r")
source("vine_copula.r")

gmm <- fit.gmm(model$posterior$samples, 4)
gmmt <- fit.gmm.truncated(model$posterior$samples, 6, bounds = bounds)

vc_ecdf <- fit.vine.copula(model$posterior$samples, fit.marginal.ecdf, bounds = bounds)
vc_mixture <- fit.vine.copula(model$posterior$samples, fit.marginal.mixture, bounds = bounds)

mdd.marginal.likelihood(gmm, model$posterior$samples, model$posterior$lposterior)
mdd.marginal.likelihood(gmmt, model$posterior$samples, model$posterior$lposterior)
mdd.marginal.likelihood(vc_ecdf, model$posterior$samples, model$posterior$lposterior)
mdd.marginal.likelihood(vc_mixture, model$posterior$samples, model$posterior$lposterior)


lp <- mdd.pdf(vc_mixture, model$posterior$samples, log = TRUE)
plot(model$posterior$lposterior, lp)


test <- model$posterior$samples[, 1]
plot(test)
plot(density(test, bw = "SJ"))

marginal <- fit.marginal.ecdf(test)
transformed <- transform.marginals(test, marginal)
reverse <- quantile(marginal$ecdfs[[1]], transformed, type = 4)
plot(density(reverse, bw = "SJ"))

library(MASS)
rlm(model$posterior$lposterior - lp ~ 1)