
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

    if (fit$type == "gmm") {
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
        for (i in 1:fit$K) {
            comp <- xmlNode("component")
            xmlAttrs(comp) <- c(weight = fit$proportions[i],
                      mean = paste(fit$centers[i,], collapse = ";"),
                      covariance = paste(apply(fit$covariances[[i]], 2, paste, collapse = ","), collapse = ";"))
            xml_root$children[[length(xml_root$children) + 1]] <- comp
        }
    } else if (fit$type == "vine.copula") {
        xmlAttrs(xml_root) <- c(type = "VineCopula")
        for (i in 1:model$nvar) {
            margin_node <- xmlNode("marginal")
            if (fit$marginal$type == "ecdf") {
                xmlAttrs(margin_node) <- c(name = bcm.model$variables[i],
                                           type = "ecdf",
                                           bw = fit$marginal$bw[i],
                                           x = paste(sort(fit$marginal$x[[i]]), collapse = ";"))
            }
            xml_root$children[[length(xml_root$children) + 1]] <- margin_node
        }
        xmlvar_rvm <- xmlNode("RVineCopula")
        xmlAttrs(xmlvar_rvm) <- c(names = paste(fit$RVM$names, collapse = ";"),
                                  matrix = paste(fit$RVM$Matrix, collapse = ";"),
                                  family = paste(fit$RVM$family, collapse = ";"),
                                  par = paste(fit$RVM$par, collapse = ";"),
                                  par2 = paste(fit$RVM$par2, collapse = ";"),
                                  maxmat = paste(fit$RVM$MaxMat, collapse = ";"),
                                  codirect = paste(as.numeric(fit$RVM$CondDistr$direct), collapse = ";"),
                                  coindirect = paste(as.numeric(fit$RVM$CondDistr$indirect), collapse = ";"))

        xml_root$children[[length(xml_root$children) + 1]] <- xmlvar_rvm

    }

    saveXML(xml_root, outfn)
}

#source("gmm.r")
#source(paste(Sys.getenv("BCM_ROOT"), "/scripts/plots_functions.r", sep = ""))

#cwd <- getwd()
#setwd("D:/Research/projects/52-approximate-mc-output-paper/")
#model <- load_model("lotka-volterra/lotka_volterra_trapperonly", "output_pt_16_1000", "lotka_volterra_trapperonly.xml")
#setwd(cwd)

#bounds <- prior_bounds_all(model, q=c(0,1))

#gmm <- fit.gmm(model$posterior$samples, 4)
#gmmbic <- gmm.BIC(model$posterior$samples)
#plot(gmmbic$BIC)
#mdd.export.bcm(model, gmmbic$fits[[5]], "trapper_posterior_gmm5.xml")

source("vine_copula.r")
vc <- fit.vine.copula(model$posterior$samples, marginalfn = fit.marginal.ecdf)

mdd.export.bcm(model, vc, "trapper_posterior_vc_ecdf.xml")

