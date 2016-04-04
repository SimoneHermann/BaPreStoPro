#' Estimation
#'
#' @description Method for the S4 classes
#' @param model.class class
#' @param ... parameters dependent on the model class
#'
setGeneric("estimate", function(model.class, ...) {
  standardGeneric("estimate")
})


########
#' Estimation for diffusion process
#'
#' @description Bayesian estimation of the parameters of the stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#'
#' @examples
#' cl_diff <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl_diff, t = t, y0 = 0.5, plot.series = TRUE)
#' est_diff <- estimate(cl_diff, t, data, 1000)
#' plot(est_diff)
#' @export
setMethod(f = "estimate", signature = "Diffusion",
          definition = function(model.class, t, data, nMCMC){

    result <- estSDE_single(t = t, X = data, prior = model.class@prior, start = model.class@start, bSDE = model.class@b.fun, sVar = model.class@sT.fun, len = nMCMC)
    he <- matrix(0, ncol(result$phi) + 1, 2)
    he[1, ] <- diagnostic(result$gamma2)
    for(i in 2:(ncol(result$phi)+1)) he[i, ] <- diagnostic(result$phi[,i-1])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.Diffusion", phi = result$phi, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

})



########
#' Estimation for mixed diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' cl <- set.to.class("mixedDiffusion", parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), b.fun = function(phi, t, x) phi*x, sT.fun = function(t, x) x)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(cl, t, data[1:20,], 2000)
#' plot(est)
#' # OU
#' b.fun <- function(phi, t, y) phi[1]-phi[2]*y
#' mu <- c(10, 5); Omega <- c(0.9, 0.01); phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' cl <- set.to.class("mixedDiffusion", parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), b.fun = b.fun, sT.fun = function(t, x) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(cl, t, data[1:20,], 2000)
#' plot(est)
#'
#' @export
setMethod(f = "estimate", signature = "mixedDiffusion",
          definition = function(model.class, t, data, nMCMC) {

    result <- estSDE(t, data, model.class@prior, model.class@start, bSDE = model.class@b.fun, sVar = model.class@sT.fun, len = nMCMC)
    he <- matrix(0, ncol(result$mu) + 1, 2)
    he[1, ] <- diagnostic(result$gamma2)
    for(i in 2:(ncol(result$mu)+1)) he[i, ] <- diagnostic(result$mu[,i-1])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.mixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

})


########
#' Estimation for noisy / hidden diffusion process
#'
#' @description Bayesian estimation of the model,
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param Npart number of particles in the particle Gibbs sampler
#'
#' @examples
#' cl <- set.to.class("hiddenDiffusion", y0.fun = function(phi, t) 0.5, parameter = list(phi = 5, gamma2 = 1, sigma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data$Z, 1000)
#' plot(est)
#' # OU
#' b.fun <- function(phi, t, y) phi[1]-phi[2]*y
#' cl <- set.to.class("hiddenDiffusion", y0.fun = function(phi, t) 0.5, parameter = list(phi = c(10, 5), gamma2 = 1, sigma2 = 0.1), b.fun = b.fun, sT.fun = function(t, x) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data$Z, 1000)
#' plot(est)
#' @export
setMethod(f = "estimate", signature = "hiddenDiffusion",
          definition = function(model.class, t, data, nMCMC, Npart = 100) {

    result <- partFiltering(t, data, prior = model.class@prior, start = model.class@start, len = nMCMC, sigmaTilde = model.class@sT.fun, Npart = Npart,
                            y0.fun = model.class@y0.fun, b.fun = model.class@b.fun, maxIt = 1)
    he <- matrix(0, ncol(result$phi) + 2, 2)
    he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$sigma2)
    for(i in 3:(ncol(result$phi)+2)) he[i, ] <- diagnostic(result$phi[,i-2])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.hiddenDiffusion", phi = result$phi, gamma2 = result$gamma2, sigma2 = result$sigma2, Y.est = result$X,
                  model = class.to.list(model.class), t = t, Z = data, burnIn = burnIn, thinning = thinning)
    return(result)

})



########
#' Estimation for noisy/hidden mixed diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @param Npart number of particles in the particle Gibbs sampler
#'
#' @examples
#' mu <- c(5, 1); Omega <- c(0.9, 0.04); phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' y0.fun <- function(phi, t) phi[2]
#' cl <- set.to.class("hiddenmixedDiffusion", y0.fun = y0.fun, b.fun = function(phi, t, y) phi[1], parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 1, sigma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' \dontrun{
#' est <- estimate(cl, t, data$Z[1:20,], 2000)
#' plot(est)
#' }
#' @export
setMethod(f = "estimate", signature = "hiddenmixedDiffusion",
          definition = function(model.class, t, data, nMCMC, Npart = 100) {

    result <-  partFiltering_mixed(t, data, prior = model.class@prior, start = model.class@start, len = nMCMC, sigmaTilde = model.class@sT.fun,
                                   y0.fun = model.class@y0.fun, b.fun = model.class@b.fun, Npart = Npart, maxIt = 1)
    he <- matrix(0, ncol(result$mu) + 2, 2)
    he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$sigma2)
    for(i in 3:(ncol(result$mu)+2)) he[i, ] <- diagnostic(result$mu[,i-2])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.hiddenmixedDiffusion", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                  sigma2 = result$sigma2, Y.est = result$X,
                  model = class.to.list(model.class), t = t, Z = data, burnIn = burnIn, thinning = thinning)
    return(result)

})





########
#' Estimation for Poisson process
#'
#' @description Bayesian estimation of a nonhomogeneous Poisson process.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain

#' @examples
#' cl <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)), Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est_NHPP <- estimate(cl, t, data$Times, 10000)
#' plot(est_NHPP)
#' @export
setMethod(f = "estimate", signature = "NHPP",
          definition = function(model.class, t, data, nMCMC) {

  if(length(t) != length(data)){
    res <- est_NHPP(jumpTimes = data, Tend = max(t), start = model.class@start, n = nMCMC, Lambda = model.class@Lambda)
    if(is.vector(res)){
      he <- diagnostic(res); burnIn <- he[1]; thinning <- min( he[2], ceiling((nMCMC-burnIn)/100) )
    }else{
      he <- matrix(0, ncol(res), 2)
      for(i in 1:nrow(res)) he[i,] <- diagnostic(res[i,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    }

    result <- new(Class = "est.NHPP", xi = t(res), jumpTimes = as.numeric(data), t = t, N = TimestoN(data, t),
                  model = class.to.list(model.class), burnIn = burnIn, thinning = thinning)

  } else {
    jumpTimes <- dNtoTimes(diff(data), t[-1])
    res <- est_NHPP(jumpTimes, Tend = max(t), start = model.class@start, n = nMCMC, Lambda = model.class@Lambda)
    if(is.numeric(res)){
      he <- diagnostic(res); burnIn <- he[1]; thinning <- min( he[2], ceiling((nMCMC-burnIn)/100) )
    }else{
      he <- matrix(0, ncol(res), 2)
      for(i in 1:ncol(res)) he[i,] <- diagnostic(res[i,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )
    }

    result <- new(Class = "est.NHPP", xi = t(res), t = t, N = data, jumpTimes = dNtoTimes(diff(data), t),
                  model = class.to.list(model.class), burnIn = burnIn, thinning = thinning)
  }
    return(result)

 })



########
#' Estimation for jump diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t + h(\eta,t,Y_t)dN_t}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#'
#' @examples
#' cl <- set.to.class("jumpDiffusion", parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)), Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(cl, t, data, 10000)
#' plot(est)
#' @export
setMethod(f = "estimate", signature = "jumpDiffusion",
          definition = function(model.class, t, data, nMCMC) {

    if(is.list(data)){
      result <- est_JD_Euler(data$Y, N=data$N, t=t, n = nMCMC, start = model.class@start, b = model.class@b.fun,
                             s = model.class@s.fun, h = model.class@h.fun, priorRatio = model.class@priorRatio,
                             Lambda = model.class@Lambda)
      he <- matrix(0, 3 + nrow(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$theta)
      for(i in 4:(3+nrow(result$xi))) he[i,] <- diagnostic(result$xi[i-3,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.jumpDiffusion", theta = result$theta, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,
                    model = class.to.list(model.class), t = t, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)
    } else {
      result <- est_JD_Euler(data, t=t, n = nMCMC, start = model.class@start, b = model.class@b.fun,
                             s = model.class@s.fun, h = model.class@h.fun,
                             priorRatio = model.class@priorRatio, Lambda = model.class@Lambda)
      he <- matrix(0, 3 + nrow(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$theta)
      for(i in 4:(3+nrow(result$xi))) he[i,] <- diagnostic(result$xi[i-3,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.jumpDiffusion", theta = result$theta, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,
                    N.est = result$N, model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }

    return(result)
})


########
#' Estimation for jump diffusion process
#'
#' @description Bayesian estimation of a stochastic process
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#'
#' @examples
#' cl <- set.to.class("Merton", parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' \dontrun{
#' est_hidden <- estimate(cl, t, data$Y, 1000)
#' plot(est_hidden)
#' }
#' @export
setMethod(f = "estimate", signature = "Merton",
          definition = function(model.class, t, data, nMCMC) {
    if(is.list(data)){
      result <- est_Merton(data$Y, data$N, t, n=nMCMC, start = model.class@start, prior = model.class@prior, Lambda=model.class@Lambda)
      he <- matrix(0, 3 + nrow(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$thetaT)
      for(i in 4:(3+nrow(result$xi))) he[i,] <- diagnostic(result$xi[i-3,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.Merton", thetaT = result$thetaT, phi = result$phi, gamma2 = result$gamma2, xi = result$xi,
                    model = class.to.list(model.class), t = t, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)

    } else {
      result <- est_Merton(data, t=t, n=nMCMC, start = model.class@start, prior = model.class@prior, Lambda=model.class@Lambda)
      he <- matrix(0, 3 + nrow(result$xi), 2)
      he[1, ] <- diagnostic(result$gamma2); he[2,] <- diagnostic(result$phi); he[3,] <- diagnostic(result$thetaT)
      for(i in 4:(3+nrow(result$xi))) he[i,] <- diagnostic(result$xi[i-3,])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.Merton", thetaT = result$thetaT, phi = result$phi, gamma2 = result$gamma2,
                    xi = result$xi, N.est = result$N,
                    model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    }
     return(result)

})


########
#' Estimation for regression model dependent on Poisson process
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(t_i, N_i, \theta) + \epsilon_i}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#'
#' @examples
#' t <- seq(0,1, by = 0.01)
#' cl <- set.to.class("reg_hiddenNHPP", fun = function(t, N, theta) theta[1]*t + theta[2]*N, parameter = list(theta = c(1,2), gamma2 = 0.1, xi = 10))
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' \dontrun{
#' est_hid <- estimate(cl, t, data$Y, 1000)
#' plot(est_hid)
#' }
#' @export
setMethod(f = "estimate", signature = "reg_hiddenNHPP",
          definition = function(model.class, t, data, nMCMC) {

    if(is.list(data)){
      result <- est_reg_hiddenNHPP(Y = data$Y, N = data$N, t = t, fun = model.class@fun, n = nMCMC, start = model.class@start,
                                   prior = model.class@prior, Lambda = model.class@Lambda)

      he <- matrix(0, ncol(result$theta) + ncol(result$xi) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$theta)+1)) he[i, ] <- diagnostic(result$theta[,i-1])
      for(i in (ncol(result$theta)+2):(ncol(result$theta) + ncol(result$xi) + 1)) he[i, ] <- diagnostic(result$xi[,i - ncol(result$theta) - 1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.reg_hiddenNHPP", theta = result$theta, gamma2 = result$gamma2, xi = result$xi,
                    model = class.to.list(model.class), t = t, Y = data$Y, N = data$N, burnIn = burnIn, thinning = thinning)

    }else{
      result <- est_reg_hiddenNHPP(data, t=t, fun = model.class@fun, n = nMCMC, start = model.class@start,
                                   prior = model.class@prior, Lambda = model.class@Lambda)
      he <- matrix(0, ncol(result$theta) + ncol(result$xi) + 1, 2)
      he[1, ] <- diagnostic(result$gamma2)
      for(i in 2:(ncol(result$theta)+1)) he[i, ] <- diagnostic(result$theta[,i-1])
      for(i in (ncol(result$theta)+2):(ncol(result$theta) + ncol(result$xi) + 1)) he[i, ] <- diagnostic(result$xi[,i - ncol(result$theta) - 1])
      burnIn <- max(he[, 1])
      thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

      result <- new(Class = "est.reg_hiddenNHPP", theta = result$theta, gamma2 = result$gamma2,
                    xi = result$xi, N.est = result$N,
                    model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)

    }
    return(result)

})


########
#' Estimation for regression model
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#'
#' @examples
#' t <- seq(0,1, by = 0.01)
#' cl <- set.to.class("Regression", fun = function(phi, t) phi[1]*t + phi[2], parameter = list(phi = c(1,2), gamma2 = 0.1))
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' @export
setMethod(f = "estimate", signature = "Regression",
          definition = function(model.class, t, data, nMCMC) {

    result <- estReg_single(t = t, y = data, prior = model.class@prior, start = model.class@start, fODE = model.class@fun, sVar = model.class@sT.fun, len = nMCMC)
    he <- matrix(0, ncol(result$phi) + 1, 2)
    he[1, ] <- diagnostic(result$gamma2)
    for(i in 2:(ncol(result$phi)+1)) he[i, ] <- diagnostic(result$phi[,i-1])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.Regression", phi = result$phi, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

 })


########
#' Estimation for mixed regression model
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(\phi_j, t_i) + \epsilon_i, \phi_j~N(\mu, \Omega)}.
#' @param model.class class of the respective model including all required information, see function set.to.class
#' @param t vector of time points
#' @param data vector or list or matrix of observation variables
#' @param nMCMC length of Markov chain
#' @examples
#' mu <- c(10, 5); Omega <- c(0.9, 0.01); phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' cl <- set.to.class("mixedRegression", parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), fun = function(phi, t) phi[1]*t + phi[2], sT.fun = function(t) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data[1:20,], 2000)
#' plot(est)
#'
#' @export
setMethod(f = "estimate", signature = "mixedRegression",
          definition = function(model.class, t, data, nMCMC) {

    result <- estReg(t, y = data, model.class@prior, model.class@start, fODE = model.class@fun, sVar = model.class@sT.fun, len = nMCMC)
    he <- matrix(0, ncol(result$mu) + 1, 2)
    he[1, ] <- diagnostic(result$gamma2)
    for(i in 2:(ncol(result$mu)+1)) he[i, ] <- diagnostic(result$mu[,i-1])
    burnIn <- max(he[, 1])
    thinning <- min( max(he[, 2]), ceiling((nMCMC-burnIn)/100) )

    result <- new(Class = "est.mixedRegression", phi = result$phi, mu = result$mu, Omega = result$Omega, gamma2 = result$gamma2,
                  model = class.to.list(model.class), t = t, Y = data, burnIn = burnIn, thinning = thinning)
    return(result)

})

