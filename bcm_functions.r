
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

mdd.export.bcm <- function(bcm.model, fit, outfn)
{
    xml_root <- xmlNode("prior")
    xml_root <- .write_variables_to_xml(bcm.model, xml_root)
    xml_root <- .write_bounds_to_xml(bcm.model, xml_root)

    if (fit$type == "kde") {
        xmlAttrs(xml_root) <- c(type = "KDE")
        for (i in 1:fit$K) {
            comp <- xmlNode("component")
            xmlAttrs(comp) <- c(weight = fit$proportions[i],
                      mean = paste(fit$centers[i,], collapse = ";"),
                      covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
            xml_root$children[[length(xml_root$children) + 1]] <- comp
        }
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
    } else if (fit$type == "vine.copula") {
        xmlAttrs(xml_root) <- c(type = "VineCopula")
        for (i in 1:bcm.model$nvar) {
            margin_node <- xmlNode("marginal")
            if (fit$marginal$type == "ecdf") {
                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                           type = "ecdf",
                                           bw = fit$marginal$bw[i],
                                           x = paste(sort(fit$marginal$x[[i]]), collapse = ";"))
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
                    xmlAttrs(margin_node) <- c(mu = paste(dist$mean, collapse = ";"),
                                               sigma = paste(dist$sd, collapse = ";"))
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
                                 b1 = fit$b1,
                                 b2 = fit$b2,
                                 x = paste(apply(fit$x, 2, paste, collapse = ","), collapse = ";"),
                                 p = paste(fit$p, collapse=";"))
        xml_root$children[[length(xml_root$children) + 1]] <- gp_params
    }

    saveXML(xml_root, outfn)
}

#source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))

#cwd <- getwd()
#setwd("D:/Research/projects/52-approximate-mc-output-paper/")
#model <- load_model("lotka-volterra/lotka_volterra_trapperonly", "output_pt_16_1000", "lotka_volterra_trapperonly.xml")
#model_t0 <- load_model("lotka-volterra/lotka_volterra_trapperonly", "output_smc_t0", "lotka_volterra_trapperonly.xml")
#setwd(cwd)

#bounds <- prior_bounds_all(model, q=c(0,1))

#source("gmm.r")
#gmm <- fit.gmm(model$posterior$samples, 4)
#gmmbic <- gmm.BIC(model$posterior$samples)
#plot(gmmbic$BIC)
#mdd.export.bcm(model, gmmbic$fits[[5]], "trapper_posterior_gmm5.xml")

#source("vine_copula.r")

#vc <- fit.vine.copula(model$posterior$samples, bounds = bounds, marginalfn = fit.marginal.ecdf)
#mdd.export.bcm(model, vc, "trapper_posterior_vc_ecdf.xml")

#vc_parametric <- fit.vine.copula(model$posterior$samples, bounds = bounds, marginalfn = fit.marginal.parametric)
#mdd.export.bcm(model, vc_parametric, "trapper_posterior_vc_parametric.xml")

#vc_mixture <- fit.vine.copula(model$posterior$samples, bounds = bounds, marginalfn = fit.marginal.mixture)
#mdd.export.bcm(model, vc_mixture, "trapper_posterior_vc_mixture.xml")

#source("gp.r")

#subset <- 1:1000
#gp <- fit.gp(model$posterior$samples[subset,], model$posterior$lposterior[subset], "se", l = c(0.01, 10), b1 = min(model$posterior$lposterior[subset]) - 5, b2 = 0, verbose=T)
#mdd.export.bcm(model, gp, "trapper_posterior_gp_se.xml")

#p2 <- evaluate.gp(gp, model$posterior$samples[subset,])
#plot(model$posterior$samples[subset, 6], p2, xlim = c(0, 2))

#testx <- matrix(NA, 500, 9)
#for (i in 1:500) {
    #testx[i,] <- gp$x[i,]
#}
#testx[, 6] <- seq(0, 2, length.out = 500)

#p2 <- evaluate.gp(gp, testx)
#plot(testx[, 6], p2)


#validate <- evaluate.gp(gp, model$posterior$samples[501:600,])
#plot(model$posterior$lposterior[501:600], validate)

#t0ix <- sample(which(!is.infinite(model_t0$posterior$llikelihood)), size = 100)
#x <- rbind(model$posterior$samples[subset,], model_t0$posterior$samples[t0ix,])
#p <- c(model$posterior$lposterior[subset], model_t0$posterior$lprior[t0ix] + model_t0$posterior$llikelihood[t0ix])

#gp <- fit.gp(x, p, "se", l=c(0.01, 10), b1=min(p), b2=0, sigman=1e-10, verbose=T)
#gp <- fit.gp(x, p, "se", l = c(0.01, 10), b1 = c(-4000, 0), b2 = c(-10,0), sigman = 1e-10, verbose = T)
#mdd.export.bcm(model, gp, "trapper_posterior_gp_se_added.xml")

#gpfull <- fit.gp(model$posterior$samples,
                 #model$posterior$lposterior,
                 #"se",
                 #l = gp$l,
                 #b1 = gp$b1,
                 #b2 = gp$b2)
#gpfull <- fit.gp(rbind(model$posterior$samples, model_t0$posterior$samples[t0ix,]),
                 #c(model$posterior$lposterior, model_t0$posterior$lprior[t0ix] + model_t0$posterior$llikelihood[t0ix]),
                 #"se",
                 #l = gp$l,
                 #b1 = gp$b1,
                 #b2 = gp$b2)
#mdd.export.bcm(model, gpfull, "trapper_posterior_gp_se.xml")