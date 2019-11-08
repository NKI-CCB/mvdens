.bw.scott.product <- function(x) {
  sds <- apply(x, 2, sd)
  n <- nrow(x)
  d <- ncol(x)
  H <- diag(sds) * n ^ (-1 / (d + 4))
  return(H)
}

#' Fit a kernel density estimate
#'
#' description
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param adjust Scalar multiplication of the bandwidth
#' @param bw.fn Function used to calculate the bandwidth. If diagonal = T, bw.fn should accept a vector of 1-dimensional values. If diagonal = F, bw.fn should accept a matrix of samples.
#' @param diagonal If true, estimate a diagonal bandwidth matrix; otherwise estimate a full bandwidth matrix.
#' @param verbose No effect.
#' @param ... Further arguments are passed to bw.fn
#' @export
#' @examples
#' x <- mvtnorm::rmvnorm(50, c(0, 0), rbind(c(1, 0.5), c(0.5, 1)))
#' fit <- fit.kde(x)
#' fit <- fit.kde(x, bw.fn = ks::Hpi, diagonal = FALSE)
fit.kde <- function(x, adjust = 1, diagonal = ncol(x) > 6,
                    bw.fn = if(ncol(x) < 6) { if(diagonal) ks::Hpi.diag else ks::Hpi } else { if(diagonal) .bw.scott.product else NULL},
                    verbose = FALSE, ...)
{
  if (is.null(bw.fn)) {
    if (!diagonal && ncol(x) > 6) {
      warning("A non-diagonal KDE is selected for >6 dimensions, for which a bandwidth selection method is not implemented. Please supply a custom bandwidth selection method or choose a diagonal bandwidth.")
    } else {
      warning("NULL bandwidth function")
    }
    return(NULL)
  }
  n <- nrow(x)
  d <- ncol(x)
  
  result <- list()
  result$type <- "kde"
  result$dim <- ncol(x)
  result$x <- x
  
  if (diagonal) {
    result$H <- .bw.scott.product(x)
  } else {
    result$H <- adjust * bw.fn(x, ...)
  }
  return(structure(result, class = "mvd.density"))
}

#' Fit a kernel density estimate after transforming the variables to an unbounded domain.
#'
#' The following transformations are used
#' [0,inf] -> log(x)
#' [-inf,0] -> log(-x)
#' [0,1] -> logit(x) = log(x / (1 - x))
#' [a,b] -> scaled logit(x) = log((a - x) / (x - b)) 
#' @param x Matrix or vector of samples. For matrices, rows are samples and columns are variables.
#' @param bounds Dx2 matrix specifying the lower and upper bound for each variable.
#' @param adjust Scalar multiplication of the bandwidth
#' @param bw.fn Function used to calculate the bandwidth. If diagonal = T, bw.fn should accept a vector of 1-dimensional values. If diagonal = F, bw.fn should accept a matrix of samples.
#' @param diagonal If true, estimate a diagonal bandwidth matrix; otherwise estimate a full bandwidth matrix.
#' @param verbose No effect.
#' @export
#' @examples
#' x <- exp(mvtnorm::rmvnorm(50, c(0, 0), rbind(c(1, 0.5), c(0.5, 1))))
#' fit <- fit.kde(x, rbind(c(0, -Inf), c(0, Inf)))
#' fit <- fit.kde(x, rbind(c(0, -Inf), c(0, Inf)), bw.fn = ks::Hpi, diagonal = FALSE)
fit.kde.transformed <- function(x, bounds, adjust = 1, bw.fn = bw.SJ, diagonal = TRUE, verbose = FALSE)
{
  result <- list()
  result$type <- "kde.transformed"
  result$transform.bounds <- bounds
  transformed <- mvd.transform_to_unbounded(x, bounds)
  result$kde <- fit.kde(transformed, adjust = adjust, bw.fn = bw.fn, diagonal = diagonal)
  return(structure(result, class = "mvd.density"))
}

.evaluate.kde <- function(fit, x, log = FALSE)
{
  tp <- rep(NA, nrow(x))
  for (i in 1:nrow(x)) {
    p <- mvtnorm::dmvnorm(fit$x, x[i,], fit$H)
    if (log) {
      tp[i] <- log(sum(p)) - log(nrow(fit$x))
    } else {
      tp[i] <- sum(p) / nrow(fit$x)
    }
  }
  return(tp)
}
