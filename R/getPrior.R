#' Builds list of prior parameters
#'
#' @description Creation of list of parameters conditional on true values \eqn{\mu} and \eqn{\gamma^2}
#' @param parameter list of parameters
#' @param model name of model
#' @return list of prior values
#' @export

getPrior <- function(parameter, model = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion",
                                          "jumpRegression", "NHPP", "Regression", "mixedRegression")){
  model <- match.arg(model)

  if(model=="jumpDiffusion"){
    prior <- list()
  }
  if(model == "Merton"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi, m.thetaT = parameter$thetaT, v.thetaT = parameter$thetaT,
                  alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="Diffusion"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
  }
  if(model=="mixedDiffusion"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="hiddenDiffusion"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma=3, beta.sigma=parameter$sigma2*2)
  }
  if(model=="hiddenmixedDiffusion"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma = 3, beta.sigma = parameter$sigma2*2)

  }
  if(model =="jumpRegression"){
    prior <- list(m.theta = parameter$theta, v.theta = parameter$theta^2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="NHPP"){
    prior <- list()
  }
  if(model == "Regression"){
    prior <- list(m.phi = parameter$phi, v.phi = parameter$phi^2, alpha.gamma = 3, beta.gamma = 2*parameter$gamma2)
  }
  if(model == "mixedRegression"){
    prior <- list(m.mu = parameter$mu, v.mu = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  prior
}
