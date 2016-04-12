

#' Gibbs sampler
#'
#' @description Bayesian estimation of the parameter of the jump diffusion process
#'   \eqn{X_t = x_0 \exp( \phi t - \gamma2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param X vector of observed variables
#' @param N vector of Poisson process variables, optional: if missing, filtering
#' @param t vector of time points
#' @param n length of Markov chain
#' @param start list of starting values, list(thetaT, gamma2)
#' @param prior list of prior values, list(mu_phi, s_phi, mu_th, s_th, alpha, beta)
#' @param Lambda intensity rate function
#' @param rangeN range of candidates for filtering N
#'
#' @return
#' \item{phi}{estimator of \eqn{\phi}}
#' \item{gamma2}{estimator of \eqn{\gamma}2}
#' \item{thetaT}{estimator of \eqn{log(1+\theta})}
#' \item{xi}{estimator of \eqn{\xi}}
#' \item{N}{estimator of latent variable \eqn{N}, if not observed}

est_Merton <- function(X, N, t, n = 1000, start, prior, Lambda, rangeN = 2){

  Delta <- diff(t)
  t.l <- t

  dlX <- diff(log(X))
  x0 <- X[1]
  X <- X[-1]
  t <- t[-1]
  if(!missing(N)){
    dN <- diff(N)
    N <- N[-1]
  }
  logX <- log(X)
  logx0 <- log(x0)
  l <- length(t)
  if(any(t == 0)){
    help <- matrix(rep(t, l-1), nrow = l-1, ncol = l-1)
    Th <- matrix(0, l-1, l-1)
    Th[lower.tri(Th, diag = TRUE)] <- t(help)[lower.tri(Th, diag = TRUE)]
    Th[upper.tri(Th)] <- help[upper.tri(help)]
    T_1 <- solve(Th)
    T_1 <- cbind(rep(0, l), rbind(rep(0, l-1), T_1))
  }else{
    help <- matrix(rep(t, l), nrow = l, ncol = l)
    Th <- matrix(0, l, l)
    Th[lower.tri(Th, diag = TRUE)] <- t(help)[lower.tri(Th, diag = TRUE)]
    Th[upper.tri(Th)] <- help[upper.tri(help)]
    T_1 <- solve(Th)
  }
  postPhi <- function(gamma2, thetaT, N){
    Vpost <- 1/( t%*%T_1%*%t/gamma2 + 1/prior$s_phi )
    mpost <- Vpost*( prior$mu_phi/prior$s_phi + 1/gamma2*t%*%T_1%*%(logX - logx0 + gamma2*t/2 - thetaT*N) )

    rnorm(1, mpost, sqrt(Vpost))
  }
  postThetaT <- function(gamma2, phi, N){
    Vpost <- 1/( N%*%T_1%*%N/gamma2 + 1/prior$s_th )
    mpost <- Vpost*( prior$mu_th/prior$s_th + 1/gamma2*N%*%T_1%*%(logX - logx0 - phi*t + gamma2*t/2) )

    rnorm(1, mpost, sqrt(Vpost))
  }
  proposalGamma <- function(gamma2, phi, thetaT, N){
    lt <- length(t)
    aPost <- prior$alpha+lt/2

    u_X <- logx0 + (phi-gamma2/2)*t + thetaT*N
    bPost <- prior$beta+(logX-u_X)%*%T_1%*%(logX-u_X)/2
    1/rgamma(1, aPost, bPost)
  }
  RatioGamma <- function(gamma2, gamma2_drawn, phi, thetaT, N){
    h1 <- t%*%T_1%*%t
    h2 <- (logX - logx0 - phi*t - thetaT*N)%*%T_1%*%t
    exp((1/gamma2_drawn - 1/gamma2)*( (gamma2-gamma2_drawn)*h2/2 + (gamma2-gamma2_drawn)*h1/8))
  }
  if(missing(N)){
    post_dN <- function(dN, i, phi, gamma2, thetaT, xi){
      Di <- dlX[i]-(phi-gamma2/2)*Delta[i]-thetaT*dN
      dLambda <- Lambda(t.l[i], xi)-Lambda(t.l[i-1], xi)
      exp(-dLambda)/prod(1:max(dN,1))*exp(-Di^2/(2*gamma2*Delta[i])+dN*log(dLambda))/sqrt(2*pi*gamma2*Delta[i])
    }
    drawN <- function(phi, gamma2, thetaT, xi, dN_old){
      dN_new <- dN_old
      cands <- lapply(dN_old, function(n) 0:(n*rangeN+5))
      prob <- lapply(1:l, function(i) sapply(cands[[i]], post_dN, i, phi, gamma2, thetaT, xi))
      diFu <- lapply(prob, cumsum)
      ind <- sapply(diFu, function(vec) any(is.na(vec)) | any(is.infinite(vec)) )
      u <- numeric(length(diFu))
      if(sum(ind) < length(diFu)){
        u[!ind] <- runif(length(diFu[!ind]), 0, sapply(diFu[!ind], max))
        for(j in (1:l)[!ind]){
          dN_new[j] <- cands[[j]][which(diFu[[j]]>=u[j])[1]]
        }
      }
      dN_new
    }
  }

  # starting values
  phi <- start$phi
  thetaT <- start$thetaT
  gamma2 <- start$gamma2
  xi <- start$xi

  if(missing(N)){
    if(is.null(start$dN)){
      dN <- diff(simN(t.l, xi, 1, start = c(t[1], 0), Lambda = Lambda)$N)
    }else{
      dN <- start$dN
    }
    N <- cumsum(dN)
    sample.N <- TRUE

    N_out <- matrix(0, l, n)

  }else{
    sample.N <- FALSE
  }

  # storage variables
  phi_out <- rep(0, n)
  thetaT_out <- rep(0, n)
  gamma2_out <- rep(0, n)
  xi_out <- matrix(0, length(xi), n)

  for(count in 1:n){

    if(sample.N){
      dN <- drawN(phi, gamma2, thetaT, xi, dN)
      N <- cumsum(dN)
      if(count %% 1000 == 0) message(paste(count, "iterations are calculated"))
    }
    xi <- est_NHPP(dNtoTimes(dN, t), t[l], xi, n = 1, Lambda = Lambda)

    phi <- postPhi(gamma2, thetaT, N)

    thetaT <- postThetaT(gamma2, phi, N)

    gamma2_drawn <- proposalGamma(gamma2, phi, thetaT, N)
    gamma2[runif(1) <= RatioGamma(gamma2, gamma2_drawn, phi, thetaT, N)] <- gamma2_drawn

    # storage
    phi_out[count] <- phi
    thetaT_out[count] <- thetaT
    gamma2_out[count] <- gamma2
    xi_out[, count] <- xi
    if(sample.N){
      N_out[, count] <- N
    }

  }
  if(sample.N){
    out <- list(phi = phi_out, gamma2 = gamma2_out, thetaT = thetaT_out, xi = xi_out, N = N_out)
  }else{
    out <- list(phi = phi_out, gamma2 = gamma2_out, thetaT = thetaT_out, xi = xi_out)
  }
  return(out)
}


#' Metropolis within Gibbs sampler
#'
#' @description Bayesian estimation of the parameter of the jump diffusion process define by SDE
#'   \eqn{dXt =  b(\phi,t,Xt) dt + s(\gamma,t,Xt) dWt + h(\theta,t,Xt) dNt}.
#' @param X vector of observed variables
#' @param N vector of Poisson process variables
#' @param t vector of time points
#' @param n length of Markov chain
#' @param start list of starting values, list(phi, theta, gamma2)
#' @param b drift function
#' @param s variance function
#' @param h jump high function
#' @param priorRatio list of functions for the prior ratio of MH step, if missing: non-informative
#' @param Lambda intensity rate function
#' @param int if Lambda is missing, one of "Weibull" or "Exp"
#' @param rangeN range for candidates for filtering N
#' @param propSd starting value for proposal standard deviation
#'
#'
#' @return
#' \item{phi}{estimator of \eqn{\phi}}
#' \item{gamma2}{estimator of \eqn{\gamma}2}
#' \item{theta}{estimator of \eqn{\theta}}
#' \item{xi}{estimator of \eqn{\xi}}
#' \item{N}{estimator of latent variable \eqn{N}, if not observed}
#' \item{prop}{storage of adapted proposal variances}

est_JD_Euler <- function(X, N, t, n = 1000, start, b, s, h, priorRatio, Lambda, int = c("Weibull","Exp"), rangeN = 2, propSd = 0.5){
  if(missing(b)) b <- function(phi, t, x) phi*x
  if(missing(s)) s <- function(gamma2, t, x) sqrt(gamma2)*x
  if(missing(h)) h <- function(theta, t, x) theta*x
  if(missing(priorRatio)){
    priorRatio <- list(
      phi = function(phi_drawn, phi_old) 1,
      theta = function(theta_drawn, theta_old) 1,
      gamma2 = function(gamma2_drawn, gamma2_old) 1 )
  }

  dX <- diff(X)
  dt <- diff(t)
  lt <- length(t)
  # starting values
  phi <- start$phi
  theta <- start$theta
  gamma2 <- start$gamma2
  xi <- start$xi


  if(missing(Lambda)){
    int <- match.arg(int)
    if(int == "Weibull"){
      Lambda <- function(t, xi){
        (t/xi[2])^xi[1]
      }
    }else{
      Lambda <- function(t, xi){
        exp(xi[1]*t+xi[2])-exp(xi[2])
      }
    }
  }

  if(missing(N)){

    post_dN <- function(dN, i, phi, gamma2, theta, xi){
      Di <- dX[i]-b(phi,t[i],X[i])*dt[i]-h(theta,t[i],X[i])*dN
      dLambda <- Lambda(t[i+1],xi)-Lambda(t[i],xi)
      exp(-dLambda)/prod(1:max(dN,1))*exp(-Di^2/(2*s(gamma2,t[i],X[i])^2*dt[i])+dN*log(dLambda))/sqrt(2*pi*s(gamma2,t[i],X[i])^2*dt[i])
    }
    drawN <- function(phi, gamma2, theta, xi, dN_old){
      dN_new <- dN_old
      cands <- lapply(dN_old, function(n) 0:(n*rangeN+5))
      prob <- lapply(1:(lt-1), function(i) sapply(cands[[i]], post_dN, i, phi, gamma2, theta, xi))
      diFu <- lapply(prob, cumsum)
      ind <- sapply(diFu, function(vec) any(is.na(vec)) | any(is.infinite(vec)) )
      u <- numeric(length(diFu))
      if(sum(ind) < length(diFu)){
        u[!ind] <- runif(length(diFu[!ind]), 0, sapply(diFu[!ind], max))
        for(j in (1:(lt-1))[!ind]){
          dN_new[j] <- cands[[j]][which(diFu[[j]] >= u[j])[1]]
        }
      }
      dN_new
    }
    if(is.null(start$N)){
      N <- simN(t, xi, 1, start = c(t[1], 0), Lambda)$N
    }else{
      N <- start$N
    }
    N_out <- matrix(0, lt, n)
    sample.N <- TRUE
  }else{
    sample.N <- FALSE
  }
  dN <- diff(N)

  likeli <- function(phi, gamma2, theta, dN){
    dnorm(dX, b(phi, t[-lt], X[-lt])*dt + h(theta, t[-lt], X[-lt])*dN, s(gamma2, t[-lt], X[-lt])*sqrt(dt))
  }

  # storage variables
  phi_out <- rep(0, n)
  theta_out <- rep(0, n)
  gamma2_out <- rep(0, n)
  xi_out <- matrix(0, length(xi), n)

  propSd_phi <- propSd*start$phi*5
  propSd_theta <- propSd/5*start$theta
  propSd_gamma2 <- propSd*start$gamma2

  prop <- matrix(0, 3, n)
  for(count in 1:n){
    if(sample.N){
      dN <- drawN(phi, gamma2, theta, xi, dN)
      N <- cumsum(c(0, dN))

    }
    xi <- est_NHPP(dNtoTimes(dN, t), t[lt], xi, n = 1, Lambda = Lambda)

    phi_drawn <- rnorm(1, phi, propSd_phi)
    ratio <- prod(likeli(phi_drawn, gamma2, theta, dN)/likeli(phi, gamma2, theta, dN))
    ratio <- ratio*priorRatio$phi(phi_drawn, phi)
    phi[runif(1) <= ratio] <- phi_drawn

    theta_drawn <- rnorm(1, theta, propSd_theta)
    ratio <- prod(likeli(phi, gamma2, theta_drawn, dN)/likeli(phi, gamma2, theta, dN))
    ratio <- ratio*priorRatio$theta(theta_drawn, theta)
    theta[runif(1) <= ratio] <- theta_drawn

    gamma2_drawn <- proposal(gamma2, propSd_gamma2)
    ratio <- prod(likeli(phi, gamma2_drawn, theta, dN)/likeli(phi, gamma2, theta, dN))
    ratio <- ratio*proposalRatio(gamma2, gamma2_drawn, propSd_gamma2)
    ratio <- ratio*priorRatio$gamma2(gamma2_drawn, gamma2)
    gamma2[runif(1) <= ratio] <- gamma2_drawn

    # storage
    phi_out[count] <- phi
    theta_out[count] <- theta
    gamma2_out[count] <- gamma2
    xi_out[, count] <- xi

    if(sample.N){
      N_out[, count] <- N
      if(count %% 1000 == 0) message(paste(count, "iterations are calculated"))
      
    }

    if (count%%50 == 0){
      propSd_phi <- ad.propSd(phi_out[(count-50+1):count], propSd_phi, count)
      propSd_theta <- ad.propSd(theta_out[(count-50+1):count], propSd_theta, count)
      propSd_gamma2 <- ad.propSd(gamma2_out[(count-50+1):count], propSd_gamma2, count)
    }
    prop[,count] <- c(propSd_phi, propSd_theta, propSd_gamma2)
  }
  if(sample.N){
    out <- list(phi = phi_out, gamma2 = gamma2_out, theta = theta_out, xi = xi_out, N = N_out, prop = prop)
  }else{
    out <- list(phi = phi_out, gamma2 = gamma2_out, theta = theta_out, xi = xi_out, prop = prop)
  }
  return(out)
}
