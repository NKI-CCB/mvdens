
.assign_kmeanspp <- function(x, K) {
  stopifnot(K >= 1)
  stopifnot(is.matrix(x))

  result <- list()
  result$k <- K

  center_ix <- rep(NA, K)
  center_ix[1] <- sample.int(nrow(x), 1)
  result$centers <- matrix(NA, K, ncol(x))
  result$centers[1,] <- x[center_ix[1],]

  if (K > 1) {
      for (i in 2:K) {
          min_dist_sq <- rep(NA, nrow(x))
          for (xi in 1:nrow(x)) {
              dist_sq <- rep(NA, i - 1)
              for (j in 1:(i - 1)) {
                  dist_sq[j] <- sum((x[xi,] - result$centers[j,]) ^ 2)
              }
              min_dist_sq[xi] <- min(dist_sq)
          }

          total <- sum(min_dist_sq)
          center_ix[i] <- sample.int(nrow(x), 1, prob = min_dist_sq / total)
          result$centers[i,] <- x[center_ix[i],]
      }
  }
  
  result$weights <- matrix(0, nrow(x), K)
  
  for (xi in 1:nrow(x)) {
    dist_sq <- apply(result$centers, 1, function(v) { return(sum((x[xi,] - v) ^ 2)) })
    assign_ix <- which.min(dist_sq)
    result$weights[xi, assign_ix] <- 1
  }
  
  return(result)
}

.mixture_expectation_step <- function(x, fit, truncated, tdist=F, tdistdf=3) {
  log_l <- matrix(NA, nrow(x), fit$k)
  for (ki in 1:fit$k) {
    if (tdist) {
      if (truncated) {
        log_l[, ki] <- tmvtnorm::dtmvt(x, fit$centers[ki,], fit$covariances[[ki]], tdistdf,
                                          lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = T) + log(fit$component_weights[ki])
      } else {
        log_l[, ki] <- mvtnorm::dmvt(x, fit$centers[ki,], fit$covariances[[ki]], tdistdf, log = T) + log(fit$component_weights[ki])
      }
    } else {
      if (truncated) {
        log_l[, ki] <- tmvtnorm::dtmvnorm(x, fit$centers[ki,], fit$covariances[[ki]],
                                lower = fit$bounds[, 1], upper = fit$bounds[, 2], log = T) + log(fit$component_weights[ki])
      } else {
        log_l[, ki] <- mvtnorm::dmvnorm(x, fit$centers[ki,], fit$covariances[[ki]], log = T) + log(fit$component_weights[ki])
      }
    }
  }
  fit$logl <- sum(apply(log_l, 1, .logsum))
  
  fit$weights <- exp(log_l)
  total_weights <- apply(fit$weights, 1, sum)
  fit$weights <- fit$weights / total_weights
  
  return(fit)
}

.mixture_maximization_step <- function(x, fit) {
  fit$component_weights <- rep(NA, fit$k)
  for (ki in 1:fit$k) {
    fit$component_weights[ki] <- sum(fit$weights[, ki]) / nrow(x)
    weighted <- cov.wt(x, fit$weights[, ki])
    if (!is.null(fit$min_cov)) {
      weighted$cov <- pmax(weighted$cov, fit$min_cov)
    }
    fit$centers[ki,] <- weighted$center
    fit$covariances[[ki]] <- weighted$cov
  }
  return(fit)
}

.fit.gmm.internal <- function(x, K, truncated, bounds, min_cov, epsilon, maxsteps, numkmeans=20, verbose, tdist, tdistdf = NA) {
  fit <- NULL
  
  for (j in 1:10) {
    if (verbose) {
      cat("Running k-means++..\n")
    }
    
    kmeansfits <- list()
    kmeansfits_logl <- rep(NA, numkmeans)
    for (i in 1:numkmeans) {
      fit <- .assign_kmeanspp(x, K)
      fit$bounds <- bounds
      fit$min_cov <- min_cov
      fit$covariances <- list()
      fit <- .mixture_maximization_step(x, fit)
      
      return_value <- try(fit <- .mixture_expectation_step(x, fit, truncated, tdist, tdistdf), silent = !verbose)
      if (inherits(return_value, "try-error")) {
        next
      }
      if (any(is.na(fit$weights))) {
        next
      }
      if (any(apply(fit$weights, 2, sum, na.rm=T) == 0)) {
        next
      }
      
      kmeansfits[[i]] <- fit
      kmeansfits_logl[i] <- fit$logl
    }
    
    if (sum(!is.na(kmeansfits_logl)) == 0) {
      warning("Could not find initial clustering with k-means++")
      return(NULL)
    }
    
    best_initial_fit <- which.max(kmeansfits_logl)
    fit <- kmeansfits[[best_initial_fit]]
    
    if (verbose) {
      cat("Best k-means likelihood:", max(kmeansfits_logl), "\n")
      cat("Starting EM..\n")
    }
    
    fit <- .mixture_maximization_step(x, fit)
    
    previous_logl <- -Inf
    singular <- F
    converged <- F
    for (i in 1:maxsteps) {
      return_value <- try(fit <- .mixture_expectation_step(x, fit, truncated, tdist, tdistdf), silent = !verbose)
      if (inherits(return_value, "try-error")) {
        singular <- T
        break
      }
      
      if (any(is.na(fit$weights))) {
        singular <- T
        break
      }
      if (any(apply(fit$weights, 2, sum, na.rm=T) == 0)) {
        singular <- T
        break
      }
      
      if (verbose) {
        cat("Log likelihood: ", fit$logl, "\n")
      }
      if (fit$logl - previous_logl < epsilon) {
        converged <- T
        break;
      } else {
        previous_logl <- fit$logl
      }
      
      fit <- .mixture_maximization_step(x, fit)
    }
    
    if (!singular) {
      break
    }
  }
  
  if (singular) {
    warning("singular solutions for 10 restarts")
    return(NULL)
  } else {
    if (!converged) {
      warning("Solution not converged in", maxsteps, ", returning optimal found solution")
    }
    return(fit)
  }
}

#' Fit a Gaussian mixture with a specific number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Integer specifying the number of components of the Gaussian mixture.
#' @param epsilon For the EM algorithm, stop when the difference in log likelihood is less than this epsilon.
#' @param maxsteps Maximum number of steps to take in the EM algorithm. When the maximum number is reached, the current fit will be returned and a warning will be issued.
#' @param verbose Display the fitting progress by showing the likelihood at every iteration.
#' @export
fit.gmm <- function(x, K, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  nparam <- K * (ncol(x) + ncol(x) * (ncol(x) + 1) / 2) + K - 1
  if (nparam >= nrow(x)) {
    warning("More parameters than samples, consider lowering K")
  }
  
  fit <- .fit.gmm.internal(x, K, truncated = F, bounds = NA, min_cov = NULL, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose, tdist = F)
  if (is.null(fit)) {
    return(NULL)
  }
  
  result <- list()
  result$type <- "gmm"
  result$K <- K
  result$centers <- fit$centers
  result$covariances <- fit$covariances
  result$proportions <- fit$component_weights
  result$log_likelihood <- fit$logl
  result$AIC <- 2 * nparam - 2 * fit$logl
  result$BIC <- log(nrow(x)) * nparam - 2 * fit$logl
  #result$assignment <- apply(fit$weights, 1, which.max)
  return(structure(result, class = "mvd.density"))
}

#' Fit a transformed Gaussian mixture with a specific number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Integer specifying the number of components of the Gaussian mixture.
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param epsilon For the EM algorithm, stop when the difference in log likelihood is less than this epsilon.
#' @param maxsteps Maximum number of steps to take in the EM algorithm. When the maximum number is reached, the current fit will be returned and a warning will be issued.
#' @param verbose Display the fitting progress by showing the likelihood at every iteration.
#' @export
fit.gmm.transformed <- function(x, K, bounds, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  result <- list()
  result$type <- "gmm.transformed"
  result$transform.bounds <- bounds
  transformed <- mvd.transform_to_unbounded(x, bounds)
  result$gmm <- fit.gmm(transformed, K, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
  if (is.null(result$gmm)) {
    return(NULL)
  }
  return(structure(result, class = "mvd.density"))
}

#' Fit a truncated Gaussian mixture with a specific number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Integer specifying the number of components of the Gaussian mixture.
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param epsilon For the EM algorithm, stop when the difference in log likelihood is less than this epsilon.
#' @param maxsteps Maximum number of steps to take in the EM algorithm. When the maximum number is reached, the current fit will be returned and a warning will be issued.
#' @param verbose Display the fitting progress by showing the likelihood at every iteration.
#' @export
fit.gmm.truncated <- function(x, K, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), min_cov = NULL, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  fit <- .fit.gmm.internal(x, K, truncated = T, bounds = bounds, min_cov = min_cov, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose, tdist = F)
  if (is.null(fit)) {
    return(NULL)
  }
  
  result <- list()
  result$type <- "gmm.truncated"
  result$K <- K
  result$centers <- fit$centers
  result$covariances <- fit$covariances
  result$proportions <- fit$component_weights
  result$bounds <- bounds
  result$log_likelihood <- fit$logl
  nparam <- K * ncol(x) * (ncol(x) + 1) / 2 + K
  result$AIC <- 2 * nparam - 2 * fit$logl
  result$BIC <- log(nrow(x)) * nparam - 2 * fit$logl
  result$assignment <- apply(fit$weights, 1, which.max)
  return(structure(result, class = "mvd.density"))
}

#' Fit a truncated t-distribution mixture with a specific number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Integer specifying the number of components of the mixture.
#' @param df Degrees of freedom of the t-distrubtions.
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param epsilon For the EM algorithm, stop when the difference in log likelihood is less than this epsilon.
#' @param maxsteps Maximum number of steps to take in the EM algorithm. When the maximum number is reached, the current fit will be returned and a warning will be issued.
#' @param verbose Display the fitting progress by showing the likelihood at every iteration.
#' @export
fit.tmixture.truncated <- function(x, K, df = 3, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), min_cov = NULL, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  fit <- .fit.gmm.internal(x, K, truncated = T, bounds = bounds, min_cov = min_cov, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose, tdist = T, tdistdf = df)
  if (is.null(fit)) {
    return(NULL)
  }
  
  result <- list()
  result$type <- "gmm.truncated"
  result$K <- K
  result$centers <- fit$centers
  result$covariances <- fit$covariances
  result$proportions <- fit$component_weights
  result$bounds <- bounds
  result$log_likelihood <- fit$logl
  nparam <- K * ncol(x) * (ncol(x) + 1) / 2 + K
  result$AIC <- 2 * nparam - 2 * fit$logl
  result$BIC <- log(nrow(x)) * nparam - 2 * fit$logl
  result$assignment <- apply(fit$weights, 1, which.max)
  return(structure(result, class = "mvd.density"))
}

#' Calculate AIC of a Gaussian mixture across a range of number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Vector specifying the number of components to test
#' @param optimal.only If TRUE, directly return only the GMM with optimal number of components; otherwise return a structure with fits for all tested number of components
#' @param epsilon See fit.gmm
#' @param maxsteps See fit.gmm
#' @param verbose Display the fitting progress.
#' @export
gmm.AIC <- function(x, K = 1:6, optimal.only = F, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  result <- list()
  result$K <- K
  result$AIC <- rep(NA, length(K))
  result$fits <- list()
  for (k in 1:length(K)) {
    if (verbose) {
      cat("Fitting K=", K[k], "\n")
    }
    
    nparam <- K[k] * (ncol(x) + ncol(x) * (ncol(x) + 1) / 2) + K[k] - 1
    if (k > 1 && nparam >= nrow(x)) {
      result$fits[[k]] <- NULL
      result$AIC[k] <- NA
    } else {
      fit <- fit.gmm(x, K[k], epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
      if (!is.null(fit)) {
        result$fits[[k]] <- fit
        result$AIC[k] <- result$fits[[k]]$AIC
      }
    }
  }
  if (optimal.only) {
    if (sum(!is.na(result$AIC)) == 0) {
      return(NULL)
    } else {
      ix <- which.min(result$AIC)
      return(result$fits[[ix]])
    }
  } else {
    return(result)
  }
}

#' Calculate AIC of a transformed Gaussian mixture across a range of number of components.
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Vector specifying the number of components to test
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param optimal.only If TRUE, directly return only the GMM with optimal number of components; otherwise return a structure with fits for all tested number of components
#' @param epsilon See fit.gmm.transformed
#' @param maxsteps See fit.gmm.transformed
#' @param verbose Display the fitting progress.
#' @export
gmm.transformed.AIC <- function(x, K = 1:6, bounds, optimal.only = F, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
  result <- list()
  result$K <- K
  result$AIC <- rep(NA, length(K))
  result$fits <- list()
  for (k in 1:length(K)) {
    if (verbose) {
      cat("Fitting K=", K[k], "\n")
    }
    nparam <- K[k] * (ncol(x) + ncol(x) * (ncol(x) + 1) / 2) + K[k] - 1
        if (k > 1 && nparam >= nrow(x)) {
            result$fits[[k]] <- NULL
            result$AIC[k] <- NA
        } else {
            result$fits[[k]] <- fit.gmm.transformed(x, K[k], bounds, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
            if (!is.null(result$fits[[k]])) {
              result$AIC[k] <- result$fits[[k]]$gmm$AIC
            }
        }
    }
    if (optimal.only) {
      if (sum(!is.na(result$AIC)) == 0) {
        return(NULL)
      } else {
        ix <- which.min(result$AIC)
        return(result$fits[[ix]])
      }
    } else {
        return(result)
    }
}

#' Calculate AIC of a truncated Gaussian mixture across a range of number of components
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param K Vector specifying the number of components to test
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param optimal.only If TRUE, directly return only the GMM with optimal number of components; otherwise return a structure with fits for all tested number of components
#' @param epsilon See fit.gmm.truncated
#' @param maxsteps See fit.gmm.truncated
#' @param verbose Display the fitting progress.
#' @export
gmm.truncated.AIC <- function(x, K = 1:6, bounds = cbind(rep(-Inf, ncol(x)), rep(Inf, ncol(x))), optimal.only = F, min_cov = NULL, epsilon = 1e-5, maxsteps = 1000, verbose = F) {
    result <- list()
    result$K <- K
    result$AIC <- rep(NA, length(K))
    result$fits <- list()
    for (k in 1:length(K)) {
        if (verbose) {
            cat("Fitting K=", K[k], "\n")
        }
        nparam <- K[k] * (ncol(x) + ncol(x) * (ncol(x) + 1) / 2) + K[k] - 1
        if (k > 1 && nparam >= nrow(x)) {
            result$fits[[k]] <- NULL
            result$AIC[k] <- NA
        } else {
            result$fits[[k]] <- fit.gmm.truncated(x, K[k], bounds = bounds, min_cov = min_cov, epsilon = epsilon, maxsteps = maxsteps, verbose = verbose)
            if (!is.null(result$fits[[k]])) {
              result$AIC[k] <- result$fits[[k]]$AIC
            }
        }
    }
    if (optimal.only) {
      if (sum(!is.na(result$AIC)) == 0) {
        return(NULL)
      } else {
        ix <- which.min(result$AIC)
        return(result$fits[[ix]])
      }
    } else {
        return(result)
    }
}
