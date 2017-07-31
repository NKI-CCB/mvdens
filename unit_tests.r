
# Marginals

source("marginals.r")

# 1-dimensional

x <- rnorm(100, 1, 2)
marginal <- fit.marginal.ecdf(x)
summary(marginal)
marginal <- fit.marginal.parametric(x)
summary(marginal)
marginal <- fit.marginal.mixture(x)
summary(marginal)

# n-dimensional

x <- cbind(rnorm(100, 1, 2), rnorm(100, 5, 2))
marginal <- fit.marginal.ecdf(x)
summary(marginal)
marginal <- fit.marginal.parametric(x)
summary(marginal)
marginal <- fit.marginal.mixture(x)
summary(marginal)