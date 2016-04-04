

########
#' Simulation of diffusion process
#'
#' @description Simulation of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param object class object of parameters: "Diffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param y0 starting point of the process
#' @param mw mesh width for finer Euler approximation
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "Diffusion",
          definition = function(object, nsim = 1, seed = NULL, t, y0, mw = 10, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- drawSDE(object@phi, object@gamma2, t, object@b.fun, fODE = function(phi, t) y0,
                                  sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }
              if(plot.series){
                plot(t, result[1,], type = "l", ylab = "Y", ylim = range(result))
                for(i in 2:nsim) lines(t, result[i,], col = i)
              }
            }else{
              result <- drawSDE(object@phi, object@gamma2, t, object@b.fun, fODE = function(phi, t) y0,
                                sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
              if(plot.series){
                plot(t, result, type = "l", ylab = "Y")
              }
            }

    return(result)
})


########
#' Simulation of diffusion process
#'
#' @description Simulation of a stochastic process
#'   \eqn{dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param object class object of parameters: "mixedDiffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param y0 starting point of the process
#' @param mw mesh width for finer Euler approximation
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedDiffusion", parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), b.fun = function(phi, t, x) phi*x, sT.fun = function(t, x) x)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "mixedDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, y0, mw = 10, plot.series = TRUE) {
      set.seed(seed)
      if(nsim > 1){
        result <- list()
        for(j in 1:nsim){
          Y <- matrix(0, nrow(object@phi), length(t))
          for(i in 1:nrow(object@phi)){
            Y[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, fODE = function(phi, t) y0,
                                  sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
          }
          result[[j]] <- Y
        }
        if(plot.series){
          plot(t, result[[1]][1,], type = "l", ylab = "Y", ylim = range(result[[1]]))
          for(i in 2:nrow(object@phi)) lines(t, result[[1]][i,])
        }
      }else{
        result <- matrix(0, nrow(object@phi), length(t))
        for(i in 1:nrow(object@phi)){
          result[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, fODE = function(phi, t) y0,
                                sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
        }
        if(plot.series){
          plot(t, result[1,], ylim = range(result), type = "l", ylab="Y")
          for(i in 2:nrow(object@phi)) lines(t, result[i,])
        }
      }

      return(result)
})



########
#' Simulation of regression model
#'
#' @description Simulation of the regression model
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i}.
#' @param object class object of parameters: "Diffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Regression", parameter = list(phi = 5, gamma2 = 0.1), fun = function(phi, t) phi*t)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "Regression",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- object@fun(object@phi, t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              }
              if(plot.series){
                plot(t, result[1,], ylab = "Y", ylim = range(result))
                for(i in 2:nsim) points(t, result[i,], col = i)
              }
            }else{
              result <- object@fun(object@phi, t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              if(plot.series){
                plot(t, result, ylab = "Y")
              }

            }
            return(result)
          })


########
#' Simulation of mixed regression model
#'
#' @description Simulation of regression model
#'   \eqn{y_i = f(\phi_j, t_i) + \epsilon_i, \phi_j\sim N(\mu, \Omega)}.
#' @param object class object of parameters: "mixedRegression"
#' @param nsim number of response vectors to simulate. Defaults = 1.
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedRegression", parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), fun = function(phi, t) phi*t, sT.fun = function(t) t)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "mixedRegression",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            if(nsim > 1){
              result <- list()
              for(j in 1:nsim){
                Y <- matrix(0, nrow(object@phi), length(t))
                for(i in 1:nrow(object@phi)){
                  Y[i,] <- object@fun(object@phi[i,], t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
                }
                result[[j]] <- Y
              }
              if(plot.series){
                plot(t, result[[1]][1,], type = "l", ylab = "Y", ylim = range(result[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result[[1]][i,])
              }
            }else{
              result <- matrix(0, nrow(object@phi), length(t))
              for(i in 1:nrow(object@phi)){
                result[i,] <- object@fun(object@phi[i,], t) + rnorm(length(t), 0, sqrt(object@gamma2*object@sT.fun(t)))
              }
              if(plot.series){
                plot(t, result[1,], ylim = range(result), type = "l", ylab="Y")
                for(i in 2:nrow(object@phi)) lines(t, result[i,])
              }
            }

            return(result)
          })



########
#' Simulation of diffusion process
#'
#' @description Simulation of a hidden stochastic process
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param object class object of parameters: "hiddenDiffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param mw mesh width for finer Euler approximation
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("hiddenDiffusion", parameter = list(phi = 0.5, gamma2 = 0.01, sigma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' @export
setMethod(f = "simulate", signature = "hiddenDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, mw = 10, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result <- matrix(0, nsim, length(t))
              for(i in 1:nsim){
                result[i,] <- drawSDE(object@phi, object@gamma2, t, object@b.fun, fODE = object@y0.fun,
                                      sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }
              obs <- apply(result, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2)))
              if(plot.series){
                plot(t, result[1,], type = "l", ylab = "Y", ylim = range(obs))
                for(i in 2:nsim) lines(t, result[i,], col = i)
                for(i in 1:nsim) points(t, obs[i,], col = i)
              }
            }else{
              result <- drawSDE(object@phi, object@gamma2, t, object@b.fun, fODE = object@y0.fun,
                                sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)

              obs <- rnorm(length(t), result, sqrt(object@sigma2))
              if(plot.series){
                plot(t, obs, ylab="simulations"); lines(t, result)
                legend("topleft", "hidden process", col = 1, lty = 1, box.lty = 0, inset = 0.01)
              }
              result <- list(Z = obs, Y = result)

            }

            return(result)
          })


########
#' Simulation of hidden mixed diffusion process
#'
#' @description Simulation of a stochastic process
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param object class object of parameters: "hiddenmixedDiffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param mw mesh width for finer Euler approximation
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' mu <- c(5, 1); Omega <- c(0.9, 0.04); phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' y0.fun <- function(phi, t) phi[2]
#' model <- set.to.class("hiddenmixedDiffusion", y0.fun = y0.fun, b.fun = function(phi, t, y) phi[1], parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 1, sigma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "hiddenmixedDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, mw = 10, plot.series = TRUE) {
            set.seed(seed)
            if(nsim > 1){
              result1 <- list()
              for(j in 1:nsim){
                Y <- matrix(0, nrow(object@phi), length(t))
                for(i in 1:nrow(object@phi)){
                  Y[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, fODE = object@y0.fun,
                                   sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
                }
                result1[[j]] <- Y
              }
              result2 <- lapply(result1, function(Y) apply(Y, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2))))
              if(plot.series){
                par(mfrow = c(2,1))
                plot(t, result1[[1]][1,], type = "l", ylab = "Y", ylim = range(result2[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result1[[1]][i,])
                plot(t, result2[[1]][1,], type = "l", ylab = "Z", ylim = range(result2[[1]]))
                for(i in 2:nrow(object@phi)) lines(t, result2[[1]][i,])
              }
            }else{
              result1 <- matrix(0, nrow(object@phi), length(t))
              for(i in 1:nrow(object@phi)){
                result1[i,] <- drawSDE(object@phi[i,], object@gamma2, t, object@b.fun, fODE = object@y0.fun,
                                      sigmaTilde = object@sT.fun, mw = mw, strictly.positive = FALSE)
              }

              result2 <- apply(result1, 2, function(vec) rnorm(length(vec), vec, sqrt(object@sigma2)))
              if(plot.series){
                par(mfrow = c(2,1))
                plot(t, result1[1,], ylim = range(result2), type = "l", ylab = "Y")
                for(i in 2:nrow(object@phi)) lines(t, result1[i,])
                plot(t, result2[1,], ylim = range(result2), type = "l", ylab = "Z")
                for(i in 2:nrow(object@phi)) lines(t, result2[i,])
              }
              result <- list(Z = result2, Y = result1)

           }

            return(result)
          })




########
#' Simulation of Poisson process
#'
#' @description Simulation of Poisson process.
#' @param object class object of parameters: "NHPP"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)), Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "NHPP",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            result <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)
            if(plot.series){
              if(nsim > 1){
                plot(t, result$N[1,], ylab="N", type="l", ylim = range(result$N))
                for(i in 2:nsim) lines(t, result$N[i,])
              }else{
                plot(t, result$N, ylab="N", type="l")
              }
            }

            return(result)
          })



########
#' Simulation of jump diffusion process
#'
#' @description Simulation of jump diffusion process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t + h(\eta,t,Y_t)dN_t}.
#' @param object class object of parameters: "jumpDiffusion"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param y0 starting point of process
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("jumpDiffusion", parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)), Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5)
#' @export
setMethod(f = "simulate", signature = "jumpDiffusion",
          definition = function(object, nsim = 1, seed = NULL, t, y0, plot.series = TRUE) {
            set.seed(seed)
            N <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)$N
            result <- list(N = N, Y = sim_JD_Euler(t, object@phi, object@theta, object@gamma2, object@b.fun, object@s.fun, object@h.fun, start = y0, N))
            if(plot.series){
              if(nsim > 1){
                par(mfrow = c(2,1))
                plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
                for(i in 2:nsim) lines(t, N[i,])
                plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
                for(i in 1:nsim) lines(t, result$Y[i,])
              }else{
                par(mfrow = c(2,1))
                plot(t, N, type = "l")
                plot(t, result$Y, type = "l", ylab = "Y")
              }
            }
            return(result)
          })





########
#' Simulation of jump diffusion process
#'
#' @description Simulation of jump diffusion process
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param object class object of parameters: "Merton"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param y0 starting point of process
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("Merton", parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5)
#' @export
setMethod(f = "simulate", signature = "Merton",
          definition = function(object, nsim = 1, seed = NULL, t, y0, plot.series = TRUE) {
            set.seed(seed)

            N <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)$N
            result <- list(N = N, Y = simY(t, object@phi, object@thetaT, object@gamma2, start = y0, N))

            if(plot.series){
              if(nsim > 1){
                par(mfrow = c(2,1))
                plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
                for(i in 2:nsim) lines(t, N[i,])
                plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
                for(i in 1:nsim) lines(t, result$Y[i,])
              }else{
                par(mfrow = c(2,1))
                plot(t, N, type = "l")
                plot(t, result$Y, type = "l", ylab = "Y")
              }
            }
            return(result)
          })



########
#' Simulation of regression model dependent on Poisson process
#'
#' @description Simulation of of the regression model
#'   \eqn{y_i = f(t_i, N_i, \theta) + \epsilon_i}.
#' @param object class object of parameters: "reg_hiddenNHPP"
#' @param nsim number of response vectors to simulate. Defaults to 1
#' @param seed optional: seed number for random number generator
#' @param t vector of time points to make predictions for
#' @param plot.series logical(1), if TRUE, simulated series are depicted grafically
#' @examples
#' model <- set.to.class("reg_hiddenNHPP", fun = function(t, N, theta) theta[1]*t + theta[2]*N, parameter = list(theta = c(1,2), gamma2 = 0.1, xi = 10))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t)
#' @export
setMethod(f = "simulate", signature = "reg_hiddenNHPP",
          definition = function(object, nsim = 1, seed = NULL, t, plot.series = TRUE) {
            set.seed(seed)

            N <- simN(t, object@xi, len = nsim, Lambda = object@Lambda)$N
            if(nsim > 1){
              result <- apply(N, 2, function(Nt) object@fun(t, Nt, object@theta) + rnorm(length(t), 0, sqrt(object@gamma2)))
            }else{
              result <- object@fun(t, N, object@theta) + rnorm(length(t), 0, sqrt(object@gamma2))
            }

            result <- list(N = N, Y = result)

            if(plot.series){
              if(nsim > 1){
                par(mfrow = c(2,1))
                plot(t, N[1,], type = "l", ylab = "N", ylim = range(N))
                for(i in 2:nsim) lines(t, N[i,])
                plot(t, result$Y[1,], type = "l", ylab = "Y", ylim = range(result$Y))
                for(i in 1:nsim) lines(t, result$Y[i,])
              }else{
                par(mfrow = c(2,1))
                plot(t, N, type = "l")
                plot(t, result$Y, ylab = "Y")
              }
            }
            return(result)
          })

##################################
##################################


#' Function for simulating the data
#'
#' @description Simulation of variables depend on the model.
#' @param t time points
#' @param cl class with informations of the parameters and the type of process
#' @param start starting point of the process
#' @param mw mesh width (in the case of SDE)
#' @param plot.series if TRUE, data is plotted
#' @return vector or matrix of simulated variables


drawData <- function(t, cl, start, mw = 10, plot.series = FALSE){

  class.name <- class(cl)[1]

  if(class.name == "jumpDiffusion"){
    N <- simN(t, cl@xi, len = 1, Lambda = cl@Lambda)$N
    result <- list(N = N, Y = sim_JD_Euler(t, cl@phi, cl@theta, cl@gamma2, cl@b.fun, cl@s.fun, cl@h.fun, start, N))
    if(plot.series){
      par(mfrow = c(2,1))
      plot(t, N, type="l")
      plot(t, result$Y, type="l", ylab="Y")
    }
  }

  if(class.name == "Merton"){
    N <- simN(t, cl@xi, len = 1, Lambda = cl@Lambda)$N
    result <- list(N = N, Y = simY(t, cl@phi, cl@thetaT, cl@gamma2, start, N))
    if(plot.series){
      par(mfrow = c(2,1))
      plot(t, N, type="l")
      plot(t, result$Y, type="l", ylab="Y")
    }

  }
  if(class.name == "Diffusion"){
    result <- drawSDE(cl@phi, cl@gamma2, t, cl@b.fun, fODE = function(phi, t) start,
                      sigmaTilde = cl@sT.fun, mw = mw, strictly.positive = FALSE)
    if(plot.series){
      plot(t, result, type="l", ylab="Y")
    }
  }
  if(class.name == "mixedDiffusion"){
    result <- matrix(0, nrow(cl@phi), length(t))
    for(i in 1:nrow(cl@phi)){
      result[i,] <- drawSDE(cl@phi[i,], cl@gamma2, t, cl@b.fun, fODE = function(phi, t) start,
                            sigmaTilde = cl@sT.fun, mw = mw, strictly.positive = FALSE)
    }
    if(plot.series){
      plot(t, result[1,], ylim = range(result), type = "l", ylab="Y")
      for(i in 2:nrow(cl@phi)) lines(t, result[i,])
    }
  }
  if(class.name == "hiddenDiffusion"){
    result <- drawSDE(cl@phi, cl@gamma2, t, cl@b.fun, fODE = cl@y0.fun,
                      sigmaTilde = cl@sT.fun, mw = mw, strictly.positive = FALSE)
    obs <- rnorm(length(t), result, sqrt(cl@sigma2))
    if(plot.series){
      plot(t, obs, ylab="simulations"); lines(t, result)
      legend("topleft", "hidden process", col = 1, lty = 1, box.lty = 0, inset = 0.01)
    }
    result <- list(Z = obs, Y = result)
  }
  if(class.name == "hiddenmixedDiffusion"){
    result <- matrix(0, nrow(cl@phi), length(t))
    for(i in 1:nrow(cl@phi)){
      result[i,] <- drawSDE(cl@phi[i,], cl@gamma2, t, cl@b.fun, fODE = cl@y0.fun,
                            sigmaTilde = cl@sT.fun, mw = mw, strictly.positive = FALSE)
    }
    obs <- apply(result, 2, function(vec) rnorm(length(vec), vec, sqrt(cl@sigma2)))
    if(plot.series){
      par(mfrow=c(2,1))
      plot(t, result[1,], ylim=range(obs), type = "l", ylab = "simulations process")
      for(i in 2:nrow(cl@phi)) lines(t, result[i,])
      plot(t, obs[1,], ylim = range(obs), type = "l", ylab = "simulations")
      for(i in 2:nrow(cl@phi)) lines(t, obs[i,])
    }
    result <- list(Z = obs, Y = result)
  }
  if(class.name == "reg_hiddenNHPP"){
    N <- simN(t, cl@xi, len = 1, Lambda = cl@Lambda)$N
    result <- cl@fun(t, N, cl@theta) + rnorm(length(t), 0, sqrt(cl@gamma2))
    result <- list(N = N, Y = result)
    if(plot.series){
      par(mfrow = c(2,1))
      plot(t, N, type="l")
      plot(t, result$Y, ylab="Y")
    }
  }
  if(class.name == "NHPP"){
    result <- simN(t, cl@xi, len = 1, Lambda = cl@Lambda)
    if(plot.series){
      plot(t, result$N, ylab="N", type="l")
    }
  }
  if(class.name == "Regression"){
    result <- drawReg(cl@phi, cl@gamma2, t, cl@fun, sVar = cl@sT.fun)
    if(plot.series){
      plot(t, result, ylab = "Y")
    }
  }
  if(class.name == "mixedRegression"){
    result <- matrix(0, nrow(cl@phi), length(t))
    for(i in 1:nrow(cl@phi)){
      result[i,] <- drawReg(cl@phi[i,], cl@gamma2, t, cl@fun, sVar = cl@sT.fun)
    }
    if(plot.series){
      plot(t, result[1,], ylim = range(result), type = "l", ylab = "Y")
      for(i in 2:nrow(cl@phi)) lines(t, result[i,])
    }

  }
  return(result)
}


#' Function for simulating diffusion process
#'
#' @description Simulation of process defined by \eqn{dYt = b(\phi, t, Y_t)dt + \gamma sigmaTilde(t, Y_t)dW_t}.
#' @param phi parameter \eqn{\phi}
#' @param gamma2 parameter \eqn{\gamma^2}
#' @param t vector of time points
#' @param b drift function
#' @param fODE function for the starting point dependent on phi, for fixed y0: fODE = function(phi, t) y0
#' @param sigmaTilde variance function \eqn{s(\gamma, t, y) = \gamma sigmaTilde(t, y)}
#' @param mw mesh width to simulate the time-continuity
#' @param strictly.positive if TRUE, only positive values for process \eqn{Y_t}
#' @return data series in t


drawSDE <- function(phi, gamma2, t, b, fODE, sigmaTilde, mw = 10, strictly.positive = TRUE){
  lt <- length(t)
  lt2 <- (lt-1)*mw+1
  t2 <- seq(min(t), max(t), length = lt2)
  dt2 <- t2[2] - t2[1]
  X <- numeric(lt2)
  X[1] <- -1

  if(strictly.positive){
    while(any(X < 0)){
      err <- rnorm(lt2-1, 0, sqrt(dt2))
      X[1] <- fODE(phi, t2[1])
      for(j in 2:lt2){
        X[j] <- X[j-1] + b(phi, t2[j-1], X[j-1])*dt2 + sqrt(gamma2)*sigmaTilde(t2[j-1], X[j-1])*err[j-1]
      }
    }
  } else {
    err <- rnorm(lt2-1, 0, sqrt(dt2))
    X[1] <- fODE(phi, t2[1])
    for(j in 2:lt2){
      X[j] <- X[j-1] + b(phi, t2[j-1], X[j-1])*dt2 + sqrt(gamma2)*sigmaTilde(t2[j-1], X[j-1])*err[j-1]
    }
  }

  X[seq(1, lt2, by = mw)]
}



#' Simulation of counting process
#'
#' @description Simulation of counting process and event times.
#' @param t vector of times
#' @param xi parameter vector \eqn{\xi}
#' @param len number of samples to be drawn
#' @param start vector: start[1] starting point time, start[2] starting point for Poisson process
#' @param Lambda intensity rate function
#' @param int if no Lambda: one out of "Weibull" or "Exp" for intensity function
#'
#' @return
#' \item{N}{Poisson process}
#' \item{Times}{event times}
#' @export
simN <- function(t, xi, len, start = c(0,0), Lambda, int = c("Weibull","Exp")){  # start=(T_n,n)
  if(missing(Lambda)){
    int <- match.arg(int)
    if(int == "Weibull"){
      Lambda <- function(t){
        (t/xi[2])^xi[1]
      }
      F_1 <- function(u, Tn_1){
        xi[2]*(Lambda(Tn_1)-log(1-u))^(1/xi[1])
      }
    }else{
      F_1 <- function(u, Tn_1){
        (log(exp(xi[1]*Tn_1+xi[2])-log(1-u))-xi[2])/xi[1]
      }
    }
    drawTn <- function(Tn_1){
      help <- F_1(runif(1), Tn_1)
      if(help < t[1]){
        return(t[1])
      }else{
        return(t[which(min(abs(help-t)) == abs(help-t))[1]])
      }
    }
  } else {

    drawTn <- function(Tn_1){
      prob <- 1-exp(-(Lambda(t, xi)-Lambda(Tn_1, xi)))
      u <- runif(1)
      t[which(abs(u-prob) == min(abs(u-prob)))]
    }
  }
  lt <- length(t)
  drawTn_vec <- function(Tn_1){
    sapply(Tn_1, drawTn)
  }

  Tn <- drawTn_vec(rep(start[1], len))
  times <- Tn

  while(min(Tn, na.rm = T) < t[lt]){
    ind <- which(Tn < t[lt])
    Tn[ind] <- drawTn_vec(Tn[ind])
    Tn[-ind] <- NA
    if(all(is.na(Tn))) break
    times <- cbind(times, Tn)
  }
  if(len > 1){
    N_out <- t(apply(times, 1, TimestoN, t)) + start[2]
    Times <- times
  }else{
    N_out <- TimestoN(times, t) + start[2]
    Times <- as.numeric(times)
  }
  list(Times = Times, N = N_out)
}

#' Simulation of Jump diffusion process
#'
#' @description Simulation of Jump diffusion process.
#' @param t vector of times
#' @param phi parameter \eqn{\phi}
#' @param thetaT parameter \eqn{\widetilde{\theta}}
#' @param gamma2 parameter \eqn{\gamma^2}
#' @param start starting point ofprocess \eqn{y_0}
#' @param N Poisson process variables in t
#'
#' @return matrix or vector
#'

simY <- function(t, phi, thetaT, gamma2, start, N){
  l <- length(t)
  if(is.matrix(N)){
    number <- nrow(N)
    Y <- sapply(1:number, function(i){
      start*exp(phi*t-gamma2*t/2+sqrt(gamma2)*cumsum(c(0, rnorm(l-1, 0, sqrt(diff(t[-l])))))+thetaT*N[i,])
    })
    return(t(Y))
  }
  if(is.numeric(N)){
    Y <- start*exp(phi*t-gamma2*t/2+sqrt(gamma2)*cumsum(c(0, rnorm(l-1, 0, sqrt(diff(t[-l])))))+thetaT*N)
    return(Y)
  }
}

#' Simulation of Jump diffusion process
#'
#' @description Simulation of Jump diffusion process.
#' @param t vector of times
#' @param phi parameter \eqn{\phi}
#' @param theta parameter \eqn{\theta}
#' @param gamma2 parameter \eqn{\gamma^2}
#' @param b.fun drift function
#' @param s.fun variance function
#' @param h.fun jump high function
#' @param start starting point \eqn{y_0}
#' @param N Poisson process variables in t
#'
#'
#' @return matrix or vector
#'

sim_JD_Euler <- function(t, phi, theta, gamma2, b.fun, s.fun, h.fun, start, N){
  l <- length(t)
  dt <- diff(t)
  if(is.matrix(N)){
    number <- nrow(N)
    Y <- matrix(start, number, l)
    for(i in 2:l){
      W <- rnorm(number, 0, sqrt(dt[i-1]))
      Y[,i] <- Y[,i-1] + b.fun(phi, t[i-1], Y[,i-1])*dt[i-1] + s.fun(gamma2, t[i-1], Y[,i-1])*W +
                  h.fun(theta, t[i-1], Y[,i-1])*(N[,i]-N[,i-1])
    }
  }
  if(is.vector(N)){
    Y <- rep(start, l)
    for(i in 2:l){
      W <- rnorm(1, 0, sqrt(dt[i-1]))
      Y[i] <- Y[i-1] + b.fun(phi, t[i-1], Y[i-1])*dt[i-1] + s.fun(gamma2, t[i-1], Y[i-1])*W +
        h.fun(theta, t[i-1], Y[i-1])*(N[i]-N[i-1])
    }
  }
  return(Y)
}


#' Simulation of regression model including the NHPP
#'
#' @description Simulation.
#' @param t vector of times
#' @param N vector of Poisson process
#' @param fun regression function
#' @param theta parameter \eqn{\theta}
#' @param gamma2 parameter\eqn{\gamma^2}
#'
#'
#' @return matrix or vector
#'

sim_reg_hiddenNHPP <- function(t, N, fun, theta, gamma2){
  lt <- length(t)

  if(is.matrix(N)){
    number <- nrow(N)
    result <- sapply(1:number, function(i) fun(t, N[i,], theta) + rnorm(lt, 0, sqrt(gamma2)) )
  }
  if(is.numeric(N)){
    result <- fun(t, N, theta) + rnorm(lt, 0, sqrt(gamma2))
  }
  return(result)
}

