library(ellipse)

source("gmm.r")
source("pdf.r")
source("transform.r")
source("gp.r")
source("marginals.r")
source("vine_copula.r")

source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep=""))

cwd <- getwd()
setwd("D:/Research/projects/52-approximate-mc-output-paper/")
model_simple <- load_sbmlpd_model("1-simple-model", "output_1")
model_extended <- load_sbmlpd_model("2-extended-model", "output_both")
setwd(cwd)

#x <- model_extended$posterior$samples[, c(1, 10)]
#bounds <- rbind(c(0, 1), c(-1, 1))
#plot(x)

x <- model_extended$posterior$samples
lposterior <- model_extended$posterior$lposterior + -3.3429
bounds <- prior_bounds_all(model_extended, q = c(0, 1))

#x <- model_simple$posterior$samples
#lposterior <- model_simple$posterior$lposterior + -4.08898
#bounds <- prior_bounds_all(model_simple, q=c(0,1))

# GMM
gmm <- fit.gmm(x, 3, verbose = T)
gmmbic <- gmm.BIC(x, 1:9, verbose = T)

plot(gmmbic$K, gmmbic$BIC)

plot(x, col=gmm$assignment+1)
for (i in 1:gmm$K) {
    lines(ellipse(gmm$covariances[[i]], centre=gmm$centers[i,], level = 0.95), col=i+1, lwd=2, lty=1)
}

p <- mdd.pdf(gmm, x, T)
plot(lposterior, p, xlim = c(-5, 15), ylim = c(-5, 15))
cor.test(p, lposterior)

# Truncated GMM

gmmt <- fit.gmm.truncated(x, 3, bounds = bounds, verbose = T)
gmmtbic <- gmm.truncated.BIC(x, 1:5, verbose = T)

plot(gmmtbic$K, gmmtbic$BIC)

plot(x, col = gmmt$assignment + 1)
for (i in 1:gmmt$K) {
    pts <- ellipse(gmmt$covariances[[i]], centre = gmmt$centers[i,], level = 0.95)
    for (j in 1:nrow(gmmt$bounds)) {
        pts[pts[, j] < gmmt$bounds[j, 1], j] <- gmmt$bounds[j, 1]
        pts[pts[, j] > gmmt$bounds[j, 2], j] <- gmmt$bounds[j, 2]
    }
    lines(pts, col = i + 1, lwd = 2, lty = 1)
}

p <- mdd.pdf(gmmt, x, T)
plot(lposterior, p, xlim = c(-5, 15), ylim = c(-5, 15))
cor.test(p, lposterior)

for (i in 1:10) {
    xt <- x[sample.int(nrow(x), 180),]

    gmmt <- fit.gmm.truncated(xt, 3, bounds = bounds, verbose = T)
}




# Transformed GMM

xform <- mdd.transform_to_unbounded(x, bounds)
gmm_xform <- fit.gmm(xform, 6, verbose = T)
p <- mdd.pdf(gmm_xform, xform, T)
p <- mdd.correct_p_for_transformation(x, bounds, p)
plot(lposterior, p, xlim=c(-5,15), ylim=c(-5,15))
cor.test(p, lposterior)


# Gaussian process

gp <- fit.gp(x, lposterior, "se", l = c(0.01, 100))
gp <- fit.gp(x, lposterior - min(lposterior) + 3, "se", l = c(0.01, 10))
gp <- fit.gp(x, lposterior - min(lposterior) + 3, "matern32", l = c(0.01, 10))
gp <- fit.gp(x, exp(lposterior), "se", l = c(0.01, 100))


y <- evaluate.gp(gp, x)
plot(y, lposterior)


# Vine copula

source("marginals.r")
marginal <- fit.marginal.ecdf(x, bounds)
summary(marginal)
marginal <- fit.marginal.parametric(x, bounds)
summary(marginal)
marginal <- fit.marginal.mixture(x, bounds)
summary(marginal)

transformed <- transform.marginals(x, marginal)
plot(transformed[,1])

source("vine_copula.r")

vc <- fit.vine.copula(x[1:100,], fit.marginal.ecdf, bounds = bounds)
p <- evaluate.vine.copula(vc, x[1:100,])
plot(p, lposterior[1:100])

vc <- fit.vine.copula(x[1:100,], fit.marginal.parametric, bounds = bounds)
p <- evaluate.vine.copula(vc, x[1:100,])
plot(p, lposterior[1:100])

vc <- fit.vine.copula(x[1:100,], fit.marginal.mixture)
p <- evaluate.vine.copula(vc, x[1:100,])
plot(p, lposterior[1:100])

vc <- fit.vine.copula(x, fit.marginal.parametric, trunclevel = 5, bounds = bounds)


transformed <- transform.marginals(x, marginal)
vc <- fit.vine.copula(1e-6+0.9999*transformed[1:101,], NULL, trunclevel = 5, verbose = T)
pvc <- mdd.pdf(vc, transformed)
p <- marginal.correct.p(marginal, x, pvc)

p <- evaluate.vine.copula(vc, x)

plot(p, lposterior)








library(mclust)
mc <- Mclust(x, modelNames = c("VVI", "VVV"))
plot(mc)

pc <- prcomp(x)
plot(pc)
plot(pc$x[, 1:2])

gmm_pc <- fit.gmm(pc$x, 6, verbose = T)
p <- mdd.pdf(gmm_pc, pc$x, T)
plot(lposterior, p, xlim = c(-5, 15), ylim = c(-5, 15))

gmm_pc12 <- fit.gmm(pc$x[, 1:8], 6, verbose = T)
plot(pc$x[, 1:2], col = gmm_pc12$assignment + 1)
for (i in 1:gmm_pc12$K) {
    lines(ellipse(gmm_pc12$covariances[[i]], centre = gmm_pc12$centers[i,], level = 0.95), col = i + 1, lwd = 2, lty = 1)
}

p <- mdd.pdf(gmm_pc12, pc$x[,1:8], T)
plot(lposterior, p, xlim = c(-5, 15), ylim = c(-5, 15))






















library(GPfit)

mins <- apply(x, 2, min)
maxs <- apply(x, 2, max)
scaled <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
for (i in 1:ncol(x)) {
    scaled[, i] <- (x[, i] - mins[i]) / (maxs[i] - mins[i])
}

gpfit <- GP_fit(scaled[1:10,], lposterior[1:10], trace=T)
pred <- predict(gpfit, x[501:1000,])
plot(lposterior[501:1000], pred$Y_hat)




library(gptk)
opts <- gpOptions()
opts$learnScales=T
gpfit <- gpCreate(ncol(x), 1, x[1:20,], as.matrix(lposterior[1:20]), options=opts)

tmp <- gpPosteriorMeanVar(gpfit, x[501:1000,])

plot(lposterior[501:1000], tmp, xlim = c(-5, 15), ylim = c(-5, 15))


gpPlot(gpfit, x[21:40,], tmp$mu, tmp$varsigma)


library(mlegp)
gpfit <- mlegp(x[1:200,], lposterior[1:200], constantMean=0)
summary(gpfit)
plot(gpfit)
par(mfrow=c(1,1))
pred <- predict(gpfit, x[501:1000,])
plot(lposterior[501:1000], pred, xlim = c(-5, 15), ylim = c(-5, 15))


kernel <- function(x1, x2, l=1) {
  factor <- -1.0 / (2*l*l)
  if (is.matrix(x1)) {
    sigma <- matrix(0, nrow(x1), nrow(x2))
    for (i in 1:nrow(sigma)) {
      v <- t(x1[i,] - t(x2))
      rsq <- apply(v*v, 1, sum)
      sigma[i,] <- exp(factor*rsq)
      #y <- sqrt(3) * sqrt(rsq) / l
      #sigma[i,j] <- (1 + y) * exp(-y)
    }
  } else {
    sigma <- matrix(0, length(x1), length(x2))
    for (i in 1:nrow(sigma)) {
      v <- x1[i] - x2
      sigma[i,] <- exp(factor*v*v)
    }
  }
  
  return(sigma)
}

train <- 1:800
test <- 801:1000

kxx <- kernel(x[train,], x[train,], 0.2)
kxsx <- kernel(x[test,],x[train,], 0.2)
y <- exp(lposterior[train])
fstar <- kxsx %*% solve(kxx, y)
fstar[fstar <= 0] <- 1e-20
plot(lposterior[test], log(fstar), xlim = c(-5, 15), ylim = c(-5, 15))

y <- exp(lposterior[train])

logml <- function(l) {
  kxx <- kernel(x[train,], x[train,], l)
  L <- chol(kxx)
  p1 <- sum(backsolve(L, y, transpose=T)^2)
  logdet <- 2.0 * sum(log(diag(L)))
  logml <- -0.5 * p1 - 0.5 * logdet - 0.5 * length(y) * log2(pi)
  print(c(l, logml))
  return(logml)
}

optimize(logml, c(0.01, 10), maximum=T)

kxx <- kernel(x[train,], x[train,], 0.254)
kxsx <- kernel(x[test,],x[train,], 0.253)
fstar <- kxsx %*% solve(kxx) %*% exp(lposterior[train])
plot(lposterior[test], log(fstar), xlim = c(-5, 15), ylim = c(-5, 15))


y <- lposterior[train] + 20

logml <- function(l) {
  kxx <- kernel(x[train,], x[train,], l)
  L <- chol(kxx)
  p1 <- sum(backsolve(L, y, transpose=T)^2)
  logdet <- 2.0 * sum(log(diag(L)))
  logml <- -0.5 * p1 - 0.5 * logdet - 0.5 * length(y) * log2(pi)
  print(c(l, logml))
  return(logml)
}

optimize(logml, c(0.01, 10), maximum=T)

kxx <- kernel(x[train,], x[train,], 0.73)
kxsx <- kernel(x[test,],x[train,], 0.73)
fstar <- kxsx %*% solve(kxx) %*% y
plot(lposterior[test], fstar - 20, xlim = c(-5, 15), ylim = c(-5, 15))
