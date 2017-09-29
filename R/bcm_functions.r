
.write_variables_to_xml <- function(model, xml_root) {
    for (i in 1:model$nvar) {
        var <- xmlNode("variable")
        logspace <- model$prior$variable_attrs[[i]]["logspace"]
        if (is.na(logspace)) {
            xmlAttrs(var) <- c(name = model$variables[i])
        } else {
            xmlAttrs(var) <- c(name = model$variables[i], logspace = paste(logspace))
        }
        xml_root$children[[length(xml_root$children) + 1]] <- var
    }
    return(xml_root)
}

.write_bounds_to_xml <- function(model, xml_root) {
    for (i in 1:model$nvar) {
        var <- xmlNode("variable_bound")
        bounds <- prior_bounds(model, i, c(0,1))
        xmlAttrs(var) <- c(name = model$variables[i], lower = bounds[1], upper = bounds[2])
        xml_root$children[[length(xml_root$children) + 1]] <- var
    }
    return(xml_root)
}

.write_transformations_to_xml <- function(model, xml_root) {
    for (i in 1:model$nvar) {
        var <- xmlNode("variable_transformation")

        distribution <- model$prior$variable_attrs[[i]]["distribution"]
        if (distribution == "normal") {
            transformation <- "none"
        } else if (distribution == "uniform") {
            a = as.numeric(model$prior$variable_attrs[[i]]["lower"]);
            b = as.numeric(model$prior$variable_attrs[[i]]["upper"]);
            if (a == 0 && b == 1) {
                transformation <- "logit"
            } else {
                transformation <- "logit_scale"
            }
        } else if (distribution == "beta") {
            transformation <- "logit"
        } else if (distribution == "exponential") {
            transformation <- "log"
        } else if (distribution == "gamma") {
            transformation <- "log"
        }

        xmlAttrs(var) <- c(name = paste(model$variables[i]),
                       transform = transformation)

        if (distribution == "uniform" && transformation == "logit_scale") {
            a = as.numeric(model$prior$variable_attrs[[i]]["lower"]);
            b = as.numeric(model$prior$variable_attrs[[i]]["upper"]);
            xmlAttrs(var) <- c(a = a, b = b)
        }

        xml_root$children[[length(xml_root$children) + 1]] <- var
    }
    return(xml_root)
}

#' Export a mvdens density object as xml that can be used as prior in the BCM software.
#'
#' See http://ccb.nki.nl/software/bcm
#' @param bcm.model A BCM model results object obtained from one of BCM's load.model functions
#' @param fit The mvdens fit object to export
#' @param outfn Output filename
#' @export
#' @examples
#' # Assumes that BCM is installed, an environment variable BCM_ROOT is specified, and there is a model named "model" with output directory "output_dir" in the current working folder
#' source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))
#' model <- load_sbmlpd_model("model", "output_dir")
#' x <- model$posterior$samples[sample_ix,]
#' p <- model$posterior$lposterior[sample_ix]
#' bounds <- prior_bounds_all(model, q = c(0, 1))
#' gmm <- gmm.BIC(x, K = 1:6, optimal.only = T)
#' mvd.export.bcm(model, gmm, "posterior_gmm.xml")
mvd.export.bcm <- function(bcm.model, fit, outfn)
{
    xml_root <- xmlNode("prior")
    xml_root <- .write_variables_to_xml(bcm.model, xml_root)
    xml_root <- .write_bounds_to_xml(bcm.model, xml_root)

    if (fit$type == "kde") {
        xmlAttrs(xml_root) <- c(type = "KDE")
        comp <- xmlNode("kde")
        xmlAttrs(comp) <- c(samples = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                            H = paste(apply(fit$H, 2, paste, collapse = ","), collapse = ";"))
        xml_root$children[[length(xml_root$children) + 1]] <- comp
    } else if (fit$type == "gmm") {
        xmlAttrs(xml_root) <- c(type = "GMM")
        for (i in 1:fit$K) {
            comp <- xmlNode("component")
            xmlAttrs(comp) <- c(weight = fit$proportions[i],
                      mean = paste(fit$centers[i,], collapse = ";"),
                      covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
            xml_root$children[[length(xml_root$children) + 1]] <- comp
        }
    } else if (fit$type == "gmm.transformed") {
        xmlAttrs(xml_root) <- c(type = "GMM")
        xml_root <- .write_transformations_to_xml(bcm.model, xml_root)
        for (i in 1:fit$gmm$K) {
            comp <- xmlNode("component")
            xmlAttrs(comp) <- c(weight = fit$gmm$proportions[i],
                      mean = paste(fit$gmm$centers[i,], collapse = ";"),
                      covariance = paste(apply(fit$gmm$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
            xml_root$children[[length(xml_root$children) + 1]] <- comp
        }
    } else if (fit$type == "gmm.truncated") {
        xmlAttrs(xml_root) <- c(type = "GMMtruncated")
        bounds <- prior_bounds_all(bcm.model, c(0, 1))
        for (i in 1:fit$K) {
            mu <- fit$centers[i,]
            sigma <- fit$covariances[[i]]
            pt <- dtmvnorm(mu, mu, sigma, lower = bounds[, 1], upper = bounds[, 2], log = T)
            pn <- dmvnorm(mu, mu, sigma, log = T)
            nc <- pt - pn

            comp <- xmlNode("component")
            xmlAttrs(comp) <- c(weight = fit$proportions[i],
                      mean = paste(fit$centers[i,], collapse = ";"),
                      covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"),
                      normalizing_factor = nc)
            xml_root$children[[length(xml_root$children) + 1]] <- comp
        }
    } else if (fit$type == "vine.copula") {
        xmlAttrs(xml_root) <- c(type = "VineCopula")
        for (i in 1:bcm.model$nvar) {
            margin_node <- xmlNode("marginal")
            if (fit$marginal$type == "ecdf") {
                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                           type = "ecdf",
                                           bw = fit$marginal$bw[i],
                                           x = paste(sort(fit$marginal$x[[i]]), collapse = ";"))
            } else if (fit$marginal$type == "ecdf.pareto") {

                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                           type = "ecdf-pareto",
                                           bw = fit$marginal$ecdf$bw[i],
                                           x = paste(sort(fit$marginal$ecdf$x[[i]]), collapse = ";"))
                if (length(fit$marginal$lower.tails[[i]]) > 0) {
                    xmlAttrs(margin_node) <- c(xmlAttrs(margin_node),
                                                lxi = fit$marginal$lower.tails[[i]]$xi,
                                                lbeta = fit$marginal$lower.tails[[i]]$beta,
                                                lu = fit$marginal$lower.tails[[i]]$u,
                                                lq = fit$marginal$lower.tails[[i]]$q,
                                                ld = fit$marginal$lower.tails[[i]]$d)
                }
                if (length(fit$marginal$upper.tails[[i]]) > 0) {
                    xmlAttrs(margin_node) <- c(xmlAttrs(margin_node),
                                           uxi = fit$marginal$upper.tails[[i]]$xi,
                                           ubeta = fit$marginal$upper.tails[[i]]$beta,
                                           uu = fit$marginal$upper.tails[[i]]$u,
                                           uq = fit$marginal$upper.tails[[i]]$q,
                                           ud = fit$marginal$upper.tails[[i]]$d)
                }
            } else if (fit$marginal$type == "parametric") {
                dist <- fit$marginal$dists[[i]]
                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i], type = dist$type)
                if (dist$type == "normal") {
                    xmlAttrs(margin_node) <- c(mu = dist$mean, sigma = dist$sd)
                } else if (dist$type == "beta") {
                    xmlAttrs(margin_node) <- c(a = dist$shape1, b = dist$shape2, min = dist$min, max = dist$max)
                } else if (dist$type == "gamma") {
                    xmlAttrs(margin_node) <- c(shape = dist$shape, scale = dist$scale)
                }
            } else if (fit$marginal$type == "mixture") {
                dist <- fit$marginal$dists[[i]]
                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i], type = dist$type, p = paste(dist$p, collapse = ";"))
                if (dist$type == "normal") {
                    xmlAttrs(margin_node) <- c(mu = paste(dist$mu, collapse = ";"),
                                               sigma = paste(dist$sigma, collapse = ";"))
                } else if (dist$type == "beta") {
                    xmlAttrs(margin_node) <- c(a = paste(dist$a, collapse = ";"),
                                               b = paste(dist$b, collapse = ";"),
                                               min = dist$min,
                                               max = dist$max)
                } else if (dist$type == "gamma") {
                    xmlAttrs(margin_node) <- c(shape = paste(dist$shape, collapse = ";"),
                                               scale = paste(dist$scale, collapse = ";"))
                }
            }
            xml_root$children[[length(xml_root$children) + 1]] <- margin_node
        }
        xmlvar_rvm <- xmlNode("RVineCopula")
        normalized_RVM <- RVineMatrixNormalize(fit$RVM)
        xmlAttrs(xmlvar_rvm) <- c(names = paste(normalized_RVM$names, collapse = ";"),
                                  matrix = paste(normalized_RVM$Matrix, collapse = ";"),
                                  family = paste(normalized_RVM$family, collapse = ";"),
                                  par = paste(normalized_RVM$par, collapse = ";"),
                                  par2 = paste(normalized_RVM$par2, collapse = ";"),
                                  maxmat = paste(normalized_RVM$MaxMat, collapse = ";"),
                                  codirect = paste(as.numeric(normalized_RVM$CondDistr$direct), collapse = ";"),
                                  coindirect = paste(as.numeric(normalized_RVM$CondDistr$indirect), collapse = ";"))

        xml_root$children[[length(xml_root$children) + 1]] <- xmlvar_rvm
    } else if (fit$type == "gp") {
        xmlAttrs(xml_root) <- c(type = "GaussianProcess")
        
        gp_params <- xmlNode("GaussianProcess")
        xmlAttrs(gp_params) <- c(kernel = fit$kernel.name,
                                 l = fit$l,
                                 s = fit$s,
                                 x = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                                 p = paste(fit$p, collapse=";"))
        xml_root$children[[length(xml_root$children) + 1]] <- gp_params
    } else if (fit$type == "resample") {
        xmlAttrs(xml_root) <- c(type = "resample")
        comp <- xmlNode("samples")
        xmlAttrs(comp) <- c(samples = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                            p = paste(fit$p, collapse = ";"))
        xml_root$children[[length(xml_root$children) + 1]] <- comp
    }

    saveXML(xml_root, outfn)
}
