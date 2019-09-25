

#' Fit a mixture of factor analyzers with a specific number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param num_factors Number of factors in the factor analyzers.
#' @param num_components Integer specifying the number of components.
#' @param epsilon For the EM algorithm, stop when the relative difference in log likelihood is less than this epsilon.
#' @param maxsteps Maximum number of steps to take in the EM algorithm. When the maximum number is reached, the current fit will be returned.
#' @export
fit.factor.mixture <- function(x, num_factors, num_components, epsilon = 1e-5, maxsteps = 100)
{
  nparam <- num_components * num_factors * (num_factors + 1) + ncol(x) * (ncol(x) + 1) + num_components - 1
  if (nparam >= nrow(x)) {
    warning("More parameters than samples, lower the number of factors or components")
    return(NULL);
  }
  
  t <- as.numeric(Sys.time())
  seed <- 1e8 * (t - floor(t))
  
  result <- list()
  result$type <- "fma"
  result$fma_res <- FactMixtAnalysis::fma(x, num_components, num_factors, eps=epsilon, it=maxsteps, seed=seed)
  result$num_components <- num_components
  result$num_factors <- num_factors
  result$weights <- result$fma_res$w[,1] # Not sure why this is a matrix of weights? all columns appear to be identical
  result$factor_loading <- result$fma_res$H
  result$factor_means <- result$fma_res$Beta[,,1]
  result$factor_covariances <- result$fma_res$sigma
  result$variable_mean <- apply(x, 2, mean)
  result$variable_covariance <- result$fma_res$psi
  result$log_likelihood <- tail(result$fma_res$lik, n=1)
  result$AIC <- result$fma_res$aic
  
  # Make covariance matrices symmetric (numerical issue?)
  if (result$num_factors > 1) {
    for (i in 1:result$num_components) {
      for (j in 1:(result$num_factors-1)) {
        for (k in (j+1):result$num_factors) {
          cov <- mean(result$factor_covariances[i,j,k], result$factor_covariances[i,j,k])
          result$factor_covariances[i,j,k] <- cov
          result$factor_covariances[i,k,j] <- cov
        }
      }
    }
  }
  
  return(structure(result, class = "mvd.density"))
}

#' Calculate AIC of a mixture of factor analyzers across a range of number of components and factors.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param factors Vector specifying the number of factors to test
#' @param components Vector specifying the number of components to test
#' @param optimal.only If TRUE, directly return only the GMM with optimal number of components; otherwise return a structure with fits for all tested number of components/factors.
#' @param epsilon See fit.factor.mixture
#' @param maxsteps See fit.factor.mixture
#' @param verbose Display the fitting progress.
#' @export
factor.mixture.AIC <- function(x, factors = 1:5, components = 1:5, optimal.only = F, epsilon = 1e-5, maxsteps = 100, verbose = F) {
  result <- list()
  result$factors <- factors
  result$components <- components
  result$AIC <- matrix(NA, length(factors), length(components))
  result$fits <- list()
  for (i in 1:length(factors)) {
    result$fits[[i]] <- list()
    
    for (j in 1:length(components)) {
      if (verbose) {
        cat("Fitting factors =", factors[i], " components =", components[j], "\n")
      }
      
      if (factors[i] < ncol(x)/2) {
        fit <- fit.factor.mixture(x, factors[i], components[j], epsilon = epsilon, maxsteps = maxsteps)
        if(!is.null(fit)) {
          result$fits[[i]][[j]] <- fit
          result$AIC[i,j] <- fit$AIC
        }
      } else {
        result$fits[[i]][[j]] <- NULL
        result$AIC[i,j] <- NA
      }
    }
  }
  if (optimal.only) {
    ixs <- which(result$AIC == min(result$AIC, na.rm=T), arr.ind = T)
    return(result$fits[[ixs[1]]][[ixs[2]]])
  } else {
    return(result)
  }
}


#' Evaluate a mixture of factor analyzers
#'
#' description
#' @param fit An mvd.density object obtained from fit.factor.mixture
#' @param x Matrix or vector of samples at which to evaluate the mixture of factor analyzers. For matrices, rows are samples and columns are variables.
#' @param log Return log probability density
#' @export
evaluate.factor.mixture <- function(fit, x, log = F) {
  ty <- scale(x, fit$variable_mean, FALSE)
  tz <- t(solve(t(fit$factor_loading) %*% solve(fit$variable_covariance) %*% fit$factor_loading) %*% t(fit$factor_loading) %*% solve(fit$variable_covariance) %*% t(ty))
  
  p <- rep(NA, nrow(x))
  for (i in 1:nrow(x)) {
    compp <- rep(NA, fit$num_components)
    if (fit$num_factors == 1) {
      for (j in 1:fit$num_components) {
        compp[j] <- fit$weights[j] * dnorm(tz[i,], fit$factor_means[j], fit$factor_covariances[j,,], log=log)
      }
    } else {
      for (j in 1:fit$num_components) {
        compp[j] <- fit$weights[j] * mvtnorm::dmvnorm(t(tz[i,]), fit$factor_means[j,], fit$factor_covariances[j,,], log=log)
      }
    }
    p[i] <- mvtnorm::dmvnorm(t(x[i,]), fit$factor_loading %*% tz[i,] + fit$variable_mean, fit$variable_covariance, log=log)
    if (log) {
      p[i] <- p[i] + .logsum(compp)
    } else {
      p[i] <- p[i] * sum(compp)
    }
  }
  return(p)
}
