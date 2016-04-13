#' Estimation function
#'
#' @description Bayesian estimation of the parameter of the model \eqn{Y_i = X_{t_i} + \epsilon_i} with
#'   \eqn{dX_t = b(\phi,t,X_t)dt + s(\gamma,t,X_t)dW_t}.
#' @param t vector of time points
#' @param y vector of observation variables
#' @param prior list of prior values
#' @param start list of starting values
#' @param est one out of "Bayes" or "MLE"
#' @param len length of Markov chain
#' @param sigmaTilde variance function \eqn{s(\gamma, t, y) = \gamma sigmaTilde(t, y)}
#' @param Npart number of particles
#' @param y0.fun function for starting point dependent on \eqn{\phi}
#' @param b.fun drift function
#' @param parPropPhi parameter for proposal standard deviation
#' @param maxIt maximal iteration of MH step of \eqn{\phi}
#'
#' @return
#' \item{phi}{Markov chain of \eqn{\phi}}
#' \item{gamma2}{Markov chain of \eqn{\gamma}2}
#' \item{sigma2}{Markov chain of \eqn{\sigma}2}
#' \item{X}{filtered process}

#' @author
#' Simone Hermann and Adeline Leclercq-Samson
partFiltering <- function(t, y, prior, start, est = c("Bayes","MLE"), len = 1000, sigmaTilde, Npart = 10, y0.fun = 1, b.fun = 1, parPropPhi = 5, maxIt = 10){
  est <- match.arg(est)

  if(missing(sigmaTilde)){
    sigmaTilde <- function(t,x) 1
  }
  lt <- length(t)
  lphi <- length(start$phi)

  phi_out <- matrix(0,len,length(start$phi))
  sigma2_out <- numeric(len)
  gamma2_out <- numeric(len)
  X_out <- matrix(0,len,length(t))

  phi <- start$phi
  sigma2 <- start$sigma2
  gamma2 <- start$gamma2

  postSigma2 <- function(alpha, beta, X, y){
    alphaPost <- alpha + length(y)/2
    betaPost <-  beta + sum((y-X)^2)/2

    1/rgamma(1, alphaPost, betaPost)
  }

  postGamma2 <- function(alpha, beta, X, t, phi, b.fun, sigmaTilde){
    ni <- length(t)
    alphaPost <- alpha + (ni-1)/2
    help <- sum( (X[-1] - X[-ni] - b.fun(phi,t[-ni],X[-ni])*(t[-1]-t[-ni]))^2/(sigmaTilde(t[-ni],X[-ni])^2*(t[-1]-t[-ni])) )

    betaPost <-  beta + help/2

    1/rgamma(1, alphaPost, betaPost)
  }

  if(est == "Bayes"){
    prior_phi <- function(phi, k){
      dnorm(phi[k],prior$mu[k],sqrt(prior$Omega[k]))
    }
    est_gamma2 <- function(phi, X){
      postGamma2(prior$alpha.gamma, prior$beta.gamma, X, t, phi, b.fun, sigmaTilde)
    }
    est_sigma2 <- function(X){
      postSigma2(prior$alpha.sigma, prior$beta.sigma, X, y)
    }
  }else{
    if(est == "MLE"){
      prior_phi <- function(phi, k) 1
      est_gamma2 <- function(phi, X){
        mean((X[-1]-X[-lt]-b.fun(phi,t[-lt],X[-lt])*diff(t))^2/(diff(t)*sigmaTilde(t[-lt],X[-lt])^2))
      }
      est_sigma2 <- function(X){
        mean((y-X)^2)
      }
      }
  }
  propSd <- abs(start$phi)/parPropPhi
  estPhi_and_X <- function(X, lastPhi, gamma2, sigma2, B_fixed, propSd){
    likeli <- function(phi,X){
      prod(dnorm(X[-1],X[-lt]+b.fun(phi, t[-lt], X[-lt])*diff(t), sqrt(gamma2*diff(t))*sigmaTilde(t[-lt],X[-lt])))
    }
    # output variables
    phi_out <- lastPhi
    phi <- lastPhi
    for(k in 1:lphi){
      for(count in 1:maxIt){
        phi[k] <- phi_out[k] + rnorm(1, 0, propSd[k])
        ratio <- prior_phi(phi, k)/prior_phi(phi_out, k)
        ratio <- ratio*likeli(phi, X)/likeli(phi_out, X)
        ratio <- ratio*dnorm(y[1], y0.fun(phi, t[1]), sqrt(sigma2))/dnorm(y[1], y0.fun(phi_out,t[1]), sqrt(sigma2))
        if(is.na(ratio)) ratio <- 0

        if(runif(1) <= ratio){
          phi_out[k] <- phi[k]
          break
        }
      }
    }
    res_SMC <- SMC(phi_out, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, X.cond = X, B_fixed = B_fixed)
    lign.B <- A.to.B(res_SMC$parents)
    indice <- sample(1:Npart, 1, prob = res_SMC$W[,lt])
    X <- res_SMC$x[indice,]

    list(phi = phi_out, X = X, B_fixed = lign.B[indice,])
  }
  res_SMC <- SMC(phi, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, conditional = FALSE)
  lign.B <- A.to.B(res_SMC$parents)
  indice <- sample(1:Npart, 1, prob =  res_SMC$W[,lt])
  X <- res_SMC$x[indice,]
  B_fixed <- lign.B[indice,]

  for(count in 1:len){

    # Particle Gibbs
    result <- estPhi_and_X(X, phi, gamma2, sigma2, B_fixed, propSd)
    phi <- result$phi
    X <- result$X
    B_fixed <- result$B_fixed

    gamma2 <- est_gamma2(phi,X)
    sigma2 <- est_sigma2(X)

    phi_out[count,] <- phi
    sigma2_out[count] <- sigma2
    gamma2_out[count] <- gamma2
    X_out[count,] <- X

    if (count%%50 == 0){
      propSd <- sapply(1:lphi, function(i){
        ad.propSd(phi_out[(count-50+1):count, i], propSd[i], count/50) })
#      print(propSd)
    }

    if (count%%1000 == 0) message(paste(count, "iterations done"))

  }
  list(phi = phi_out, sigma2 = sigma2_out, gamma2 = gamma2_out, X = X_out)
}

#' Estimation function
#'
#' @description Bayesian estimation of the parameter of the model \eqn{Y_{ij} = X_{t_{ij}} + \epsilon_{ij}} with
#'   \eqn{dX_t = b(\phi_j, t, X_t)dt + s(\gamma, t, X_t)dW_t}.
#' @param t vector of time points
#' @param y matrix of observation variables
#' @param prior list of prior values
#' @param start list of starting values
#' @param len length of Markov chain
#' @param sigmaTilde variance function \eqn{s(\gamma, t, y) = \gamma sigmaTilde(t, y)}
#' @param y0.fun function for starting point dependent on \eqn{\phi}
#' @param b.fun drift function
#' @param Npart number of particles
#' @param parPropPhi parameter for proposal standard deviation
#' @param maxIt maximal iteration of MH step of \eqn{\phi}
#'
#' @return
#' \item{phi}{Markov chain of \eqn{\phi}}
#' \item{mu}{Markov chain of \eqn{\mu}}
#' \item{Omega}{Markov chain of \eqn{\Omega}}
#' \item{gamma2}{Markov chain of \eqn{\gamma^2}}
#' \item{sigma2}{Markov chain of \eqn{\sigma^2}}
#' \item{X}{filtered process}

#' @author
#' Simone Hermann and Adeline Leclercq-Samson
partFiltering_mixed <- function(t, y, prior, start, len = 1000, sigmaTilde, y0.fun = 1, b.fun = 1, Npart = 10, parPropPhi = 5, maxIt = 10){
  # if y is a matrix and t a vector, they will be transformed to lists
  if(is.matrix(y)){
    if(nrow(y) == length(t)){
      y <- t(y)
    }else{
      if(ncol(y)!=length(t)){
        stop("length of t has to be equal to the columns of y")
      }
    }
    t1 <- t
    y1 <- y
    t <- list()
    y <- list()
    for(i in (1:nrow(y1))){
      t[[i]] <- t1
      y[[i]] <- y1[i,]
    }
  }

  if(missing(sigmaTilde)){
    sigmaTilde <- function(t,x) 1
  }
  lphi <- length(start$mu)

  postSigma2 <- function(alpha, beta, X, y){    # X and y lists with length n, each list entry is one series of length n_i, alpha, beta prior parameters
    alphaPost <- alpha + sum(unlist(lapply(y,length)))/2
    betaPost <-  beta + sum((unlist(y)-unlist(X))^2)/2

    1/rgamma(1, alphaPost, betaPost)
  }

  postGamma2 <- function(alpha, beta, X, t, phi, b.fun, sigmaTilde){ #  phi matrix, X, t lists, alpha, beta prior parameters
    n <- length(X)
    alphaPost <- alpha + sum(sapply(X,length)-1)/2

    help <- numeric(n)
    for(i in 1:n){
      ni <- length(t[[i]])
      dt <- t[[i]][-1]-t[[i]][-ni]
      help[i] <- sum( (X[[i]][-1] - X[[i]][-ni] - b.fun(phi[i,],t[[i]][-ni],X[[i]][-ni])*dt)^2/(sigmaTilde(t[[i]][-ni],X[[i]][-ni])^2*dt) )
    }
    betaPost <-  beta + sum(help)/2

    1/rgamma(1, alphaPost, betaPost)
  }
  propSd <- abs(start$mu)/parPropPhi
  postPhii_Xi <- function(y, t, X, lastPhi, mu, Omega, gamma2, sigma2, B_fixed, propSd){
    lt <- length(t)
    likeli <- function(phi,X){
      prod(dnorm(X[-1], X[-lt]+b.fun(phi, t[-lt], X[-lt])*diff(t), sqrt(gamma2*diff(t))*sigmaTilde(t[-lt],X[-lt])))  # true transition density???
    }
    # output variables
    phi_out <- lastPhi
    phi <- lastPhi
    for(k in 1:lphi){
      for(count in 1:maxIt){
        phi[k] <- phi_out[k] + rnorm(1, 0, propSd[k])
        ratio <- dnorm(phi[k], mu[k], sqrt(Omega[k]))/dnorm(phi_out[k], mu[k], sqrt(Omega[k]))
        ratio <- ratio*likeli(phi, X)/likeli(phi_out, X)
        ratio <- ratio*dnorm(y[1], y0.fun(phi, t[1]), sqrt(sigma2))/dnorm(y[1], y0.fun(phi_out, t[1]), sqrt(sigma2))
        if(is.na(ratio)) ratio <- 0

        if(runif(1) <= ratio){
          phi_out[k] <- phi[k]
          break
        }
      }
    }
    res_SMC <- SMC(phi_out, gamma2, sigma2, Npart, t, y, b.fun, y0.fun, sigmaTilde, X.cond = X, B_fixed = B_fixed)
    lign.B <- A.to.B(res_SMC$parents)
    indice <- sample(1:Npart, 1, prob = res_SMC$W[,lt])
    X <- res_SMC$x[indice,]

    list(phi = phi_out, X = X, B_fixed = lign.B[indice,])

  }

  n <- length(y)

  phi_out <- list()
  mu_out <- matrix(0, len, length(start$mu))
  Omega_out <- matrix(0, len, length(start$mu))
  sigma2_out <- numeric(len)
  gamma2_out <- numeric(len)
  X_out <- list()

  phi <- start$phi
  sigma2 <- start$sigma2
  gamma2 <- start$gamma2
  mu <- start$mu
  Omega <- postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)
  X <- list()
  B_fixed <- list()

  for(i in 1:n){

    if(dnorm(y[[i]][1], y0.fun(phi[i,], t[[i]][1]), sqrt(sigma2)) == 0){
      stop("bad starting values")
    }

    result <- SMC(phi[i,], gamma2, sigma2, Npart, t[[i]], y[[i]], b.fun, y0.fun, sigmaTilde, conditional = FALSE)
    lign.B <- A.to.B(result$parents)
    indice <- sample(1:Npart, 1, prob = result$W[,length(t[[i]])])
    X[[i]] <- result$x[indice,]
    B_fixed[[i]] <- lign.B[indice,]
  }

  for(count in 1:len){

    for(i in 1:n){
      help <- postPhii_Xi(y[[i]], t[[i]], X[[i]], phi[i,], mu, Omega, gamma2, sigma2, B_fixed[[i]], propSd)
      X[[i]] <- help$X
      phi[i,] <- help$phi
      B_fixed[[i]] <- help$B_fixed
    }

    mu <- postmu(phi, prior$m, prior$v, Omega)
    Omega <- postOmega(prior$alpha.omega, prior$beta.omega, phi, mu)

    sigma2 <- postSigma2(prior$alpha.sigma, prior$beta.sigma, X, y)
    gamma2 <- postGamma2(prior$alpha.gamma, prior$beta.gamma, X, t, phi, b.fun, sigmaTilde)

    phi_out[[count]] <- phi
    mu_out[count,] <- mu
    Omega_out[count,] <- Omega
    sigma2_out[count] <- sigma2
    gamma2_out[count] <- gamma2
    X_out[[count]] <- X

    if (count%%50 == 0){
      propSd <- sapply(1:lphi, function(i){
        ad.propSd(phi_ij(phi_out[(count-50+1):count], 1, i), propSd[i], count/50) })
#      print(propSd)
    }

    if (count%%100 == 0) message(paste(count, "iterations done"))
  }
  list(phi = phi_out, mu = mu_out, Omega = Omega_out, sigma2 = sigma2_out, gamma2 = gamma2_out, X = X_out)
}


####################
# helping functions...

A.to.B <- function(A){
  ntimes <- dim(A)[2] + 1;
  Npart <- dim(A)[1];
  B <- matrix(0, Npart, ntimes)
  B[,ntimes] <- 1:Npart;
  for (t in seq(ntimes, 2, by = -1)){
    B[,t-1] <- A[B[,t],t-1]}
  return(B)
}

SMC <- function(phi, gamma2, sigma2, Npart, times, y, b.fun, y0.fun, sigmaTilde, conditional = TRUE, X.cond, B_fixed){

  ntimes <- length(times)
  dt <-  diff(times)
  x <- matrix(y0.fun(phi, times[1]), Npart, ntimes)
  w <- Weights <- matrix(1,Npart,ntimes)
  parents <-  matrix(1, Npart, ntimes-1)

  w[,1] <- dnorm(y[1], mean = x[,1], sd = sqrt(sigma2))
  Weights[,1] <- w[,1]/sum(w[,1]) # or only 1/Npart...
  if(conditional){
    x[B_fixed[1],1] <- X.cond[1]
  }else{
    if(missing(X.cond)) X.cond <- numeric(ntimes)
  }

  for(n in 2:ntimes){
    if(conditional){
      ##############################
      # sampling of A_{n-1}^{-B_n^K}
      ##############################
      set.parents  <- (1:Npart)[-B_fixed[n]]

      On_1 <- rmultinom(1, Npart-1, Weights[,n-1])
      O <- On_1[B_fixed[n-1]] + 1
      he <- sample(set.parents, O-1)

      parents[B_fixed[n], n-1] <- B_fixed[n-1]
      parents[he, n-1] <- B_fixed[n-1]
      parents[-c(B_fixed[n],he), n-1] <- sample(set.parents, Npart-O, replace = TRUE, prob =  Weights[-B_fixed[n],n-1])
    }else{
      set.parents <- 1:Npart
      parents[, n-1] <- sample(1:Npart, Npart, replace = T, prob = Weights[,n-1])
    }

    x[,1:(n-1)] <- x[parents[,n-1], 1:(n-1)]
    x.past <- x[,n-1]

    ##########################################
    # simulation of particles
    ##########################################

    Gpost <- 1/(1/(gamma2*sigmaTilde(times[n-1], x.past)^2*dt[n-1]) + 1/sigma2)
    mpost <- Gpost*((x.past+b.fun(phi, times[n-1], x.past)*dt[n-1])/(gamma2*sigmaTilde(times[n-1], x.past)^2*dt[n-1])
                    + y[n]/sigma2)
    x.new <- rnorm(Npart, mpost, sqrt(Gpost))
    x[set.parents,n] <- x.new[set.parents]
    x[-set.parents,n] <- X.cond[n]
    wn <- dnorm(x[,n], x.past+b.fun(phi, times[n-1], x.past)*dt[n-1], sqrt(gamma2*sigmaTilde(times[n-1], x.past)^2*dt[n-1]))*
      dnorm(y[n], x[,n], sqrt(sigma2))/dnorm(x[,n], mpost, sqrt(Gpost) )

    w[,n] <- wn
    Weights[,n] <- w[,n]/sum(w[,n])
  }
  return(list(W = Weights, x = x , parents = parents) )
}

