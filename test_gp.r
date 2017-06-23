library(rgl)
library(mvtnorm)

x <- rnorm(100)
p <- dnorm(x)

plot(x, p)
plot(x, log(p))

library(mlegp)
gpfit <- mlegp(x, p)
plot(gpfit)
par(mfrow=c(1,1))

xstar <- seq(-10, 10, length.out=200)
pred <- predict(gpfit, as.matrix(xstar))
plot(xstar, pred)




library(gptk)
opts <- gpOptions()
opts$learnScales=T
gpfit <- gpCreate(1, 1, as.matrix(x), as.matrix(p), options=opts)



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

xstar <- seq(-10, 10, length.out=200)

kxx <- kernel(x,x)
kxxs <- kernel(x,xstar)
kxsx <- kernel(xstar,x)
kxsxs <- kernel(xstar,xstar)

fstar <- kxsx %*% solve(kxx + 1e-8 * diag(nrow(kxx))) %*% (log(p)+20)
plot(xstar, fstar-20, typ="l")
points(x, log(p))
points(xstar, dnorm(xstar, log=T), col="red")

plot(xstar, exp(fstar-20), typ="l")
points(x, p)
points(xstar, dnorm(xstar, log=F), col="red")

fstar <- kxsx %*% solve(kxx + 1e-8 * diag(nrow(kxx))) %*% p
plot(xstar, fstar, typ="l")
points(x, p)
points(xstar, dnorm(xstar, log=), col="red")


# 2d
mu <- c(1,2)
sigma <- rbind(c(1,0.5), c(0.5,2.0))
x <- rmvnorm(100, mean=mu, sigma=sigma)
plot(x)
p <- dmvnorm(x, mean=mu, sigma=sigma)
plot3d(cbind(x, p))

xstar <- matrix(runif(1000,-5,5), ncol=2)

kxx <- kernel(x,x)
kxxs <- kernel(x,xstar)
kxsx <- kernel(xstar,x)
kxsxs <- kernel(xstar,xstar)

fstar <- kxsx %*% solve(kxx + 1e-8 * diag(nrow(kxx))) %*% p
plot3d(cbind(x, p), col="red")
plot3d(cbind(xstar, fstar), add=T)
points(xstar, dnorm(xstar, log=), col="red")

# N-d
D <- 16
logscale <- F
l <- 1

x <- rmvnorm(1000, mean=rep(0,D))
p <- dmvnorm(x, mean=rep(0,D))

#xstar <- matrix(runif(1000*D,-5,5), ncol=D)
xstar <- rmvnorm(1000, mean=rep(0,D))

kxx <- kernel(x,x, l)
kxsx <- kernel(xstar,x, l)

if (logscale) {
  fstar <- kxsx %*% solve(kxx) %*% log(p)
} else {
  fstar <- kxsx %*% solve(kxx) %*% p
}
ftrue <- dmvnorm(xstar, mean=rep(0,D), log=logscale)

plot(fstar, ftrue)

plot(log(fstar), log(ftrue))
