#' Bayesian estimation in mixed stochastic differential equations
#'
#' @description Bayesian estimation of the random effects \eqn{\phi_i} in the mixed SDE
#'  \eqn{dY_i(t) = b(\phi_i, t, Y_i(t))dt + \gamma s(t, Y_i(t)) dW_i(t), \phi_i~N(\mu, \Omega), i=1,...,n} and the parameters
#'  \eqn{\mu, \Omega, \gamma^2}.
#' @param t vector of observation times
#' @param y matrix or list of the n trajectories
#' @param prior list of prior parameters - list(m, v, alpha.omega, beta.omega, alpha.gamma, beta.gamma)
#' @param start list of starting values
#' @param bSDE b(phi, t, x) drift function
#' @param sVar variance function s^2
#' @param ipred which of the n trajectories is the one to be predicted
#' @param cut the index how many of the ipred-th series are used for estimation
#' @param len number of iterations of the MCMC algorithm - chain length
#' @param mod model out of {Gompertz, Richards, logistic, Weibull, Paris, Paris2}, only used instead of bSDE
#' @param propPar proposal standard deviation of phi is |start$mu|*propPar
#'
#' @return
#' \item{phi}{samples from posterior of \eqn{\phi}}
#' \item{mu}{samples from posterior of \eqn{\mu}}
#' \item{Omega}{samples from posterior of \eqn{\Omega}}
#' \item{gamma2}{samples from posterior of \eqn{\gamma^2}}
#'
#' @examples
#'
#' mu <- 0.5; Omega <- 0.1; phi <- matrix(rnorm(100, mu, sqrt(Omega)))
#' cl <- set.to.class("mixedDiffusion", parameter=list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.1)
#' Y <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#'
#' prior <- list(m = 1, v = 1, alpha.omega = 3, beta.omega = 2*Omega, alpha.gamma = 3, beta.gamma = 2*0.01)
#' start <- list(phi = phi, mu = mu, gamma2 = 0.01)
#' est <- estSDE(t, Y, bSDE = cl@@b.fun, prior = prior, start = start, sVar = cl@@sT.fun, len = 1000)
#' par(mfrow = c(1,2))
#' plot(est$mu, type = "l", ylab = expression(mu))
#' plot(est$gamma2, type = "l", ylab = expression(gamma^2))

#' @details
#' Simulation from the posterior distribution of the random effect from n independent trajectories of the SDE (the Brownian motions \eqn{W1,...,Wn} are independent).
#'
#' @keywords estimation
#' @references Hermann et al. (2015)
#' @export

estSDE <- function(t, y, prior, start, bSDE, sVar, ipred = 1, cut, len = 1000, mod = c("Gompertz", "logistic", "Weibull", "Richards", "Paris", "Paris2"), propPar = 0.2){
  mod <- match.arg(mod)
  if(is.matrix(y)){
    if(nrow(y) == length(t)){
      y <- t(y)
    }else{
      if(ncol(y) != length(t)){
        print("length of t has to be equal to the columns of y")
        break
      }
    }
    if(missing(cut)) cut <- length(t)
    t1 <- t
    y1 <- y
    t <- list()
    y <- list()
    for(i in (1:nrow(y1))[-ipred]){
      t[[i]] <- t1
      y[[i]] <- y1[i,]
    }
    t[[ipred]] <- t1[1:cut]
    y[[ipred]] <- y1[ipred, 1:cut]
  }
  if(missing(bSDE)) bSDE <- getFun("SDE", mod)
  if(missing(sVar)) sVar <- function(t, x) 1

#  if(Omega=="diag"){
    postOm <- function(phi, mu){
      postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
    }
#   }else{
#     postOm <- function(phi,mu){
#       postOmega_matrix(prior$alpha.omega,phi,mu)
#     }
#   }
  propSd <- abs(start$mu)*propPar
  postPhii <- function(lastPhi, mu, Omega, gamma2, X, t, propSd){  # X, t vektoren
    lt <- length(t)
    dt <- diff(t)
    phi_old <- lastPhi
    phi_drawn <- phi_old + rnorm(length(mu), 0, propSd)
#    ratio <- dmvnorm(phi_drawn, mu, as.matrix(Omega)) / dmvnorm(phi_old, mu, as.matrix(Omega))
    ratio <- prod(dnorm(phi_drawn, mu, sqrt(Omega)) / dnorm(phi_old, mu, sqrt(Omega)) )
    ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt)))
    if(is.na(ratio)) ratio <- 0
    if(runif(1) <= ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }
  n <- length(y)

  postGamma2 <- function(phi){
    alphaPost <- prior$alpha.gamma + sum(sapply(y, length)-1)/2
    help <- numeric(n)
    for(i in 1:n){
      ni <- length(t[[i]])
      delta <- diff(t[[i]])
      help[i] <- sum( (y[[i]][-1] - y[[i]][-ni] - bSDE(phi[i,], t[[i]][-ni], y[[i]][-ni])*delta)^2/(sVar(t[[i]][-ni], y[[i]][-ni])^2*delta) )
    }
    betaPost <-  prior$beta.gamma + sum(help)/2
    1/rgamma(1, alphaPost, betaPost)
  }

  phi_out <- list()
  mu_out <- matrix(0, len, length(start$mu))
  Omega_out <- matrix(0, len, length(start$mu))
  gamma2_out <- numeric(len)

  phi <- start$phi
  gamma2 <- start$gamma2
  mu <- start$mu
  Omega <- postOm(phi, mu)

  for(count in 1:len){

    for(i in 1:n){
      phi[i,] <- postPhii(phi[i,], mu, Omega, gamma2, y[[i]], t[[i]], propSd)
    }
    mu <- postmu(phi, prior$m, prior$v, Omega)
    Omega <- postOm(phi, mu)
    gamma2 <- postGamma2(phi)

    phi_out[[count]] <- phi
    mu_out[count,] <- mu
    Omega_out[count,] <- Omega
    gamma2_out[count] <- gamma2

    if (count%%50 == 0){
      propSd <- sapply(1:length(phi[1,]), function(i){
        ad.propSd(sapply(phi_out[(count-50+1):count], function(mat) mat[1,i]), propSd[i], count/50) })
#            print(propSd)
    }

  }
  list(phi = phi_out, mu = mu_out, Omega = Omega_out, gamma2 = gamma2_out)
}

#' Bayesian estimation in stochastic differential equations
#'
#' @description Bayesian estimation of the parameters in the SDE
#'  \eqn{dY(t)= b(\phi, t, Y(t))dt + \gamma s(t, Y(t)) dW(t)}.
#' @param t vector of observation times
#' @param X vector of the M trajectories
#' @param prior list of prior parameters - list(mu, Omega, alpha, beta)
#' @param start list of starting values
#' @param bSDE drift function
#' @param sVar variance function
#' @param len number of iterations of the MCMC algorithm
#'
#' @examples
#' cl <- set.to.class("Diffusion", parameter = list(phi = 1, gamma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' Y <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#'
#' prior <- list(mu = 1, Omega = 1, alpha = 3, beta = 2*0.1)
#' start <- list(phi = 1, gamma2 = 0.1)
#' est <- estSDE_single(t, Y, bSDE = cl@@b.fun, prior = prior, start = start, sVar = cl@@sT.fun, len = 1100)
#' par(mfrow = c(1,2))
#' plot(est$phi, type = "l", ylab = expression(phi))
#' plot(est$gamma2, type = "l", ylab = expression(gamma^2))
#' @details
#' Simulation from the posterior distribution of the random effect from n independent trajectories of the SDE (the Brownian motions \eqn{W1,...,Wn} are independent).
#' @return
#' \item{phi}{estimator of \eqn{\phi}}
#' \item{gamma2}{estimator of \eqn{\gamma^2}}
#' @export

estSDE_single <- function(t, X, prior, start, bSDE, sVar, len = 1000){  # p liste von prior-parametern, len=Anzahl von Ziehungen

  propSd <- abs(prior$mu)/5
  lt <- length(t)
  dt <- t[-1] - t[-lt]
  lphi <- length(propSd)

  postPhi <- function(lastPhi, gamma2, propSd){
    phi_old <- lastPhi

    phi_drawn <- phi_old + rnorm(lphi, 0, propSd)
    ratio <- dmvnorm(phi_drawn, prior$mu, as.matrix(prior$Omega)) / dmvnorm(phi_old, prior$mu, as.matrix(prior$Omega))
    ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old, t[-lt], X[-lt])*dt, sqrt(gamma2*sVar(t[-lt], X[-lt])^2*dt)))
    if(is.na(ratio)){ratio <- 0}
    if(runif(1) < ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }
  postGamma2 <- function(phi){
    alphaPost <- prior$alpha + (lt-1)/2
    betaPost <-  prior$beta + sum( (X[-1] - X[-lt] - bSDE(phi, t[-lt], X[-lt])*dt)^2/(sVar(t[-lt], X[-lt])^2*dt) )/2
    1/rgamma(1, alphaPost, betaPost)
  }

  phi_out <- matrix(0, len, lphi)
  gamma2_out <- numeric(len)

  phi <- start$phi
  gamma2 <- start$gamma2

  for(count in 1:len){

    phi <- postPhi(phi, gamma2, propSd)
    gamma2 <- postGamma2(phi)

    phi_out[count, ] <- phi
    gamma2_out[count] <- gamma2

    if (count%%50 == 0){
      propSd <- sapply(1:length(phi), function(i){
        ad.propSd(phi_out[(count-50+1):count, i], propSd[i], count/50) })
#      print(propSd)
    }


  }
  list(phi = phi_out, gamma2 = gamma2_out)
}

