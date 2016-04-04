#' Builds list of prior parameters
#'
#' @description Creation of list of parameters conditional on true values \eqn{\mu} and \eqn{\gamma^2}
#' @param parameter list of parameters
#' @param model name of model
#' @return list of prior values
#' @export

getPrior <- function(parameter, model = c("jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion",
                                          "reg_hiddenNHPP", "NHPP", "Regression", "mixedRegression")){
  model <- match.arg(model)

  if(model=="jumpDiffusion"){
    prior <- list()
  }
  if(model == "Merton"){
    prior <- list(mu_phi = parameter$phi, s_phi = parameter$phi, mu_th = parameter$thetaT, s_th = parameter$thetaT,
                  alpha = 3, beta = parameter$gamma2*2)
  }
  if(model=="Diffusion"){
    prior <- list(mu = parameter$phi, Omega = diag(parameter$phi^2, length(parameter$phi), length(parameter$phi)), alpha = 3, beta = 2*parameter$gamma2)
  }
  if(model=="mixedDiffusion"){
    prior <- list( m = parameter$mu, v = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  if(model=="hiddenDiffusion"){
    prior <- list(mu = parameter$phi, Omega = parameter$phi^2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma=3, beta.sigma=parameter$sigma2*2)
  }
  if(model=="hiddenmixedDiffusion"){
    prior <- list( m = parameter$mu, v = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2, alpha.sigma = 3, beta.sigma = parameter$sigma2*2)

  }
  if(model =="reg_hiddenNHPP"){
    prior <- list(mu = parameter$theta, Omega = parameter$theta^2, alpha = 3, beta = parameter$gamma2*2)
  }
  if(model=="NHPP"){
    prior <- list()
  }
  if(model == "Regression"){
    prior <- list(mu = parameter$phi, Omega = diag(parameter$phi^2, length(parameter$phi), length(parameter$phi)), alpha = 3, beta = 2*parameter$gamma2)
  }
  if(model == "mixedRegression"){
    prior <- list( m = parameter$mu, v = parameter$mu^2, alpha.omega = rep(3, length(parameter$mu)),
                   beta.omega = parameter$Omega*2, alpha.gamma = 3, beta.gamma = parameter$gamma2*2)
  }
  prior
}
