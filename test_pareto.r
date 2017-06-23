library(spd)

source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))

cwd <- getwd()
setwd("D:/Research/projects/52-approximate-mc-output-paper/")
model_simple <- load_sbmlpd_model("1-simple-model", "output_1")
model_extended <- load_sbmlpd_model("2-extended-model", "output_both")
setwd(cwd)


fig1_samples <- model_extended$posterior$samples[, c(1, 10)]
bounds <- rbind(c(0, 1), c(-1, 1))

plot(fig1_samples)
plot(fig1_samples[, 1])
plot(fig1_samples[, 2])

spds <- apply(fig1_samples, 2, spdfit, upper = 0.9, lower = 0.1)

#plot(spds[[2]])



x <- fig1_samples[, 2]
plot(density(x, bw = "SJ"))

reflected <- 1 + (1 - x)
xr <- c(x, reflected)

plot(density(xr, bw = "SJ", from = -1, to = 1))

spd <- spdfit(x, upper = 0.95, lower = 0.05, type = "pwm")
plot(spd, which = 2, add = T)

rdy <- density(xr, bw = "SJ", from = -1, to = 1)
lines(rdy$x, rdy$y * 2, col = "green", lwd = 2)


u <- quantile(x, 0.9)
my_gpdtailfit <- function(x, u) {
    exceedances <- x[x > u]
    excess <- exceedances - u
    Nu <- length(excess)
    xbar <- mean(excess)

    s2 <- var(excess)
    xi0 <- -0.5 * (((xbar * xbar) / s2) - 1)
    beta0 <- 0.5 * xbar * (((xbar * xbar) / s2) + 1)
    theta <- c(xi0, beta0)

    negloglik <- function(theta, tmp) {
        cat(theta, "\n")
        xi <- theta[1]
        beta <- theta[2]
        cond1 <- beta <= 0
        cond2 <- (xi <= 0) && (max(tmp) > (-beta / xi))
        if (cond1 || cond2) {
            f <- 1e+06
        } else {
            y <- logb(1 + (xi * tmp) / beta)
            y <- y / xi
            f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
        }
        f
    }

    fit <- optim(theta, negloglik, tmp = excess)
    #fit <- optim(theta, negloglik, lower = c(-Inf, 1e-6), upper = c(Inf, Inf), tmp = excess, method = "L-BFGS-B")
    #fit <- optimize(negloglikul, lower=0, upper=100, tmp = excess)
    #return(c(-upperlimit / fit$minimum, fit$minimum / upperlimit))
    return(fit$par)
}

.dgpd <- function(x, mu = 0, beta = 1, xi = 0, log = FALSE) {
    stopifnot(min(beta) > 0)
    stopifnot(length(xi) == 1)
    # Density:
    d <- (x - mu) / beta
    nn <- length(d)
    beta <- rep(beta, length.out = nn)
    #index <- (d > 0 & ((1 + xi * d) > 0)) | is.na(d)
    if (xi == 0) {
        d[index] <- log(1 / beta[index]) - d[index]
        d[!index] <- -Inf
    } else {
        d[index] <- log(1 / beta[index]) - (1 / xi + 1) * log(1 + xi * d[index])
        d[!index] <- -Inf
    }

    # Log:
    if (!log) d <- exp(d)
    # Return Value:
    return(d)
}

dgpd <- function(x, xi, mu = 0, beta = 1) {
    (beta ^ (-1)) * (1 + (xi * (x - mu)) / beta) ^ ((-1 / xi) - 1)
}

res <- my_gpdtailfit(x, quantile(x, 0.9))

u <- quantile(x, 0.1)
res <- my_gpdtailfit(-x, -u)

nx <- seq(u, 1, length.out = 100)
lines(nx, dgpd(nx, xi = res[1], mu = u, beta = res[2]) * 0.05, typ = "l", lwd = 3, col = "cyan")

nx <- seq(-1, u, length.out = 100)
lines(nx, dgpd(-nx, xi = res[1], mu = -u, beta = res[2]) * 0.1, typ = "l", lwd = 3, col = "cyan")

dgpd(-nx, xi = res[1], mu = -u, beta = res[2])







x <- seq(0, 1, length.out = 200)
plot(x, dgpd(x, xi = -0.933, mu = 0, beta = 0.058), typ = "l")

x <- seq(-1, 6, length.out = 200)
plot(x, dgpd(x, xi = -1.5, mu = 0, beta = 5), typ = "l")


test <- function(x, xi = 1, mu = 0, sigma = 1) {
    z <- (x - mu) / sigma
    if (xi == 0) {
        return((1 / sigma) * exp(-z))
    } else {
        return((1 / sigma) * (1 + xi * z) ^ (-1 / xi - 1))
    }
}
par(mfrow = c(1, 1))
x <- seq(-1, 5, length.out = 200)
plot(x, test(x, xi = -0.2, mu = 0, sigma = 1), typ = "l", ylim = c(0, 5))
lines(x, test(x, xi = -0.5, mu = 0, sigma = 2.5), col = "red")
lines(x, test(x, xi = -1.0, mu = 0, sigma = 5), col = "red")
lines(x, test(x, xi = -2.0, mu = 0, sigma = 10), col = "red")
lines(x, test(x, xi = -0.1, mu = 0, sigma = 0.5), col = "blue")
lines(x, test(x, xi = -0.05, mu = 0, sigma = 0.25), col = "blue")

plot(x, test(x, xi = -5, mu = 0, sigma = 25), typ = "l")

gpdFit