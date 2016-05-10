#' Builds classes
#'
#' @description Defines classes
#' @param class.name name of model class
#' @param parameter list of parameter values
#' @param prior optional list of prior parameters
#' @param start optional list of starting values
#' @param b.fun drift function b
#' @param s.fun variance function s
#' @param h.fun jump high function h
#' @param sT.fun variance function \eqn{\widetilde{s}}
#' @param y0.fun function for the starting point, if dependent on parameter
#' @param fun regression function
#' @param Lambda intensity rate of Poisson process
#' @param priorRatio list of functions for prior ratios, only for jumpDiffusion, is missing: non-informative estimation
#' @return class
#' @examples
#' set.to.class("jumpDiffusion")
#' cl_jd <- set.to.class("jumpDiffusion", 
#'              parameter = list(theta = 0.1, phi = 0.01, gamma2 = 0.1, xi = 3))
#' summary(class.to.list(cl_jd))
#'

#' @export
set.to.class <- function(class.name = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression", "mixedRegression"),
                         parameter, prior, start, b.fun, s.fun, h.fun, sT.fun, y0.fun, fun, Lambda, priorRatio){
  class.name <- match.arg(class.name)

  if(missing(parameter)){
    print("parameter has to be a list of ")
    return(defaults(class.name))
  }


  if(missing(prior)) prior <- getPrior(parameter, class.name)
  if(missing(start)) start <- parameter

  if(missing(b.fun)) b.fun <- function(phi, t, y) phi
  if(missing(s.fun)) s.fun <- function(gamma2, t, y) sqrt(gamma2)
  if(missing(h.fun)) h.fun <- function(theta, t, y) theta

  if(missing(sT.fun)) sT.fun <- function(t, y) 1
  if(missing(y0.fun)) y0.fun <- function(phi, t) 1
  if(missing(fun)){
    if(class.name %in% c("Regression", "mixedRegression")) fun1 <- function(phi, t) phi
    if(class.name == "jumpRegression") fun2 <- function(t, N, theta) theta
  }else{
    if(class.name %in% c("Regression", "mixedRegression")) fun1 <- fun
    if(class.name == "jumpRegression") fun2 <- fun
  }

  if(missing(Lambda)) Lambda <- function(t, xi) xi*t


  if(class.name == "jumpDiffusion"){
    if(missing(priorRatio)){
      priorRatio <- list(
        phi = function(phi_drawn, phi_old) 1,
        theta = function(theta_drawn, theta_old) 1,
        gamma2 = function(gamma2_drawn, gamma2_old) 1,
        xi = function(xi_drawn, xi_old) 1)
    }
    return(new(Class = class.name, theta = parameter$theta, phi = parameter$phi, gamma2 = parameter$gamma2, xi = parameter$xi,
               b.fun = b.fun, s.fun = s.fun, h.fun = h.fun, Lambda = Lambda, priorRatio = priorRatio,
               prior = prior, start = start))
  }
  if(class.name == "Merton"){
    if(missing(priorRatio)) priorRatio <- function(xi_drawn, xi_old) 1
    return(new(Class = class.name, thetaT = parameter$thetaT, phi = parameter$phi, gamma2 = parameter$gamma2, xi = parameter$xi,
               Lambda = Lambda, prior = prior, start = start, priorRatio = priorRatio))
  }
  if(class.name == "Diffusion"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2,
                 b.fun = b.fun, sT.fun = sT.fun,
                 prior = prior, start = start))
  }
  if(class.name == "mixedDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
                 y0.fun = y0.fun, b.fun = b.fun, sT.fun = sT.fun,
                 prior = prior, start = start))
  }
  if(class.name == "hiddenDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2, sigma2 = parameter$sigma2,
                 b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun,
                 prior = prior, start = start))
  }
  if(class.name == "hiddenmixedDiffusion"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
                 sigma2 = parameter$sigma2, b.fun = b.fun, sT.fun = sT.fun, y0.fun = y0.fun,
                 prior = prior, start = start))
  }
  if(class.name == "jumpRegression"){
    return(new(Class = class.name, theta = parameter$theta, gamma2 = parameter$gamma2, xi = parameter$xi, fun = fun2, Lambda = Lambda,
               sT.fun = sT.fun, prior = prior, start = start))
  }
  if(class.name == "NHPP"){
    if(missing(priorRatio)) priorRatio <- function(xi_drawn, xi_old) 1
    return(new(Class = class.name, xi = parameter$xi, Lambda = Lambda, start = parameter$xi, priorRatio = priorRatio))
  }
  if(class.name == "Regression"){
    return(new(Class = class.name, phi = parameter$phi, gamma2 = parameter$gamma2,
               fun = fun1, sT.fun = sT.fun,
               prior = prior, start = start))
  }
  if(class.name == "mixedRegression"){
    return(new(Class = class.name, phi = parameter$phi, mu = parameter$mu, Omega = parameter$Omega, gamma2 = parameter$gamma2,
               fun = fun1, sT.fun = sT.fun,
               prior = prior, start = start))
  }

}

defaults <- function(class.name = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression", "mixedRegression")){

  class.name <- match.arg(class.name)
  if(class.name == "jumpDiffusion"){
    name.vec <- c("theta", "phi", "gamma2", "xi")
  }
  if(class.name == "Merton"){
    name.vec <- c("thetaT", "phi", "gamma2", "xi")
  }
  if(class.name == "Diffusion"){
    name.vec <- c("phi", "gamma2")
  }
  if(class.name == "mixedDiffusion"){
    name.vec <- c("phi", "mu", "Omega", "gamma2")
  }
  if(class.name == "hiddenDiffusion"){
    name.vec <- c("phi", "gamma2", "sigma2")
  }
  if(class.name == "hiddenmixedDiffusion"){
    name.vec <- c("phi", "mu", "Omega", "gamma2", "sigma2")
  }
  if(class.name == "jumpRegression"){
    name.vec <- c("theta", "gamma2", "xi")
  }
  if(class.name == "NHPP"){
    name.vec <- "xi"
  }
  if(class.name == "Regression"){
    name.vec <- c("phi", "gamma2")
  }
  if(class.name == "mixedRegression"){
    name.vec <- c("phi", "mu", "Omega", "gamma2")
  }
  name.vec
}

#' Builds a list from class
#'
#' @description Class to list
#' @param cl class
#' @return list
#' @export

class.to.list <- function(cl){
  
  sN <- slotNames(cl)
  res <- lapply(sN, function(name) slot(cl, name))
  names(res) <- sN
  res
}
