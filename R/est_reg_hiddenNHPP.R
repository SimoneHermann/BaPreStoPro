


#' Gibbs sampler
#'
#' @description Bayesian estimation of the parameter of the regression model
#'   \eqn{y_i = f(t_i, N_i, \theta) + \epsilon_i}.
#' @param Y vector of observation variables
#' @param N vector of Poisson process variables
#' @param t vector of time points
#' @param fun regression function
#' @param n length of Markov chain
#' @param start list of starting values
#' @param prior list of prior values
#' @param Lambda intensity rate function
#' @param int one out of "Weibull" or "Exp", if Lambda is missing
#' @param rangeN range for candidates of filtering N, if unobserved
#'
#' @return
#' \item{theta}{Markov chains of \eqn{\theta}}
#' \item{gamma2}{Markov chains of \eqn{\gamma^2}}
#' \item{N}{if hidden: Markov chains of \eqn{N}}
#' \item{xi}{Markov chains of \eqn{\xi}}

est_reg_hiddenNHPP <- function(Y, N, t, fun, n = 1000, start, prior, Lambda, int = c("Weibull","Exp"), rangeN = 2){
  lt <- length(t)
  # starting values
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
        xi[2]*exp(xi[1]*t)-xi[2]
      }
    }
  }

  if(missing(N)){

    Npart <- 10
    A.to.B = function(A){
      B <- matrix(0, Npart, lt)
      B[,lt] <- 1:Npart;
      for (t in seq(lt, 2, by = -1)){
        B[,t-1] <- A[B[,t], t-1]}
      return(B)
    }

    CSMC = function(theta, gamma2, xi, N.cond, B_fixed, conditional = TRUE){# conditional SMC
      # N.cond = the old samples

      x <- matrix(0, Npart, lt)
      w <- matrix(1, Npart, lt)
      W <- matrix(1, Npart, lt)
      parents <-  matrix(1, Npart, lt-1)

      # initialisation
      x[,1] <- 0  # poisson always starts in 0
      if(conditional){
        x[B_fixed[1], 1] <- N.cond[1]
      }else{
        N.cond <- numeric(lt)
      }
      w[,1] <- dnorm(Y[1], mean = fun(t[1], x[,1], theta), sd = sqrt(gamma2))
      W[,1] <- w[,1]/sum(w[,1])

      for (n in 2:lt){
        if(conditional){
          set.parents  <- (1:Npart)[-B_fixed[n]]
          
          On_1 <- rmultinom(1, Npart-1, W[,n-1])
          O <- On_1[B_fixed[n-1]] + 1
          he <- sample(set.parents, O-1)
          
          parents[B_fixed[n], n-1] <- B_fixed[n-1]
          parents[he, n-1] <- B_fixed[n-1]
          parents[-c(B_fixed[n],he), n-1] <- sample(set.parents, Npart-O, replace = TRUE, prob =  W[-B_fixed[n], n-1])
          
        }else{
          set.parents <- 1:Npart
          parents[, n-1] <- sample(1:Npart, Npart, replace = TRUE, prob = W[,n-1])
          
        }

        
        x[,1:(n-1)] <- x[parents[,n-1], 1:(n-1)]
        x.past <- x[,n-1]
        

        dr <- function(Ni, dNi_old){
          cands <- 0:(dNi_old*rangeN + 5)
          prob <- dnorm(Y[n], fun(t[n], Ni+cands, theta), sqrt(2*gamma2))*
            dpois(Ni+cands, Lambda(t[n], xi))
          diFu <- cumsum(prob)
          u <- runif(1, 0, max(diFu))
          Ni + cands[which(diFu >= u)[1]]
        }
        x.new <- sapply(x.past, dr, diff(N.cond)[n-1])

        x[set.parents, n] <- x.new[set.parents]
        x[-set.parents, n] <- N.cond[n]
        wn <- dpois(x[, n], Lambda(t[n], xi))*
          dnorm(Y[n], fun(t[n], x[, n], theta), sqrt(gamma2))/
          (dnorm(Y[n], fun(t[n], x[, n], theta), sqrt(2*gamma2))*
           dpois(x[, n], Lambda(t[n], xi)))
        if(sum(wn) == 0) wn <- rep(1, length(wn))
        w[, n] <- wn
        W[, n] <- w[, n]/sum(w[, n])
      }
      lign.B <- A.to.B(parents)
      indice <- sample(1:Npart, 1, prob = W[, lt])
      X <- x[indice, ]
      B_fixed <- lign.B[indice,]
      return(list(N = X, B_fixed = B_fixed) )
    }

    if(is.null(start$N)){
#      N <- simN(t, xi, 1, start = c(t[1], 0), Lambda = Lambda)$N
      he <- CSMC(theta, gamma2, xi, conditional = FALSE)
      N <- he$N
      B_fixed <- he$B_fixed
    }else{
      N <- start$N
    }
    N_out <- matrix(0, lt, n)
    sample.N <- TRUE
  }else{
    sample.N <- FALSE
  }
  propSd <- abs(prior$mu)/20
  ltheta <- length(start$theta)


  sVar <- function(t) 1   # generalize ?

  postTheta <- function(N, gamma2, lastPhi, propSd){  # bad acceptance rate?
    phi_old <- lastPhi
    phi_drawn <- rnorm(ltheta, phi_old, propSd)
    ratio <- prod(dnorm(phi_drawn, prior$mu, sqrt(prior$Omega)))/prod(dnorm(phi_old, prior$mu, sqrt(prior$Omega)))
    ratio <- ratio* prod(dnorm(Y, fun(t, N, phi_drawn), sqrt(gamma2*sVar(t)))/dnorm(Y, fun(t, N, phi_old), sqrt(gamma2*sVar(t))))
    if(is.na(ratio)) ratio <- 0
    if(runif(1) < ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }

  postGamma2 <- function(theta, N){
    alphaPost <- prior$alpha + lt/2
    betaPost <-  prior$beta + sum((Y-fun(t, N, theta))^2/sVar(t))/2
    1/rgamma(1, alphaPost, betaPost)
  }

  theta_out <- matrix(0, n, length(theta))
  gamma2_out <- numeric(n)
  xi_out <- matrix(0, n, length(xi))

  for(count in 1:n){
    if(sample.N){
      he <- CSMC(theta, gamma2, xi, N, B_fixed)
      N <- he$N
      B_fixed <- he$B_fixed
    }
    xi <- est_NHPP(dNtoTimes(diff(N), t[-1]), t[lt], xi, n = 1, Lambda = Lambda)

    theta <- postTheta(N, gamma2, theta, propSd)

    gamma2 <- postGamma2(theta, N)

    # storage
    theta_out[count, ] <- theta
    gamma2_out[count] <- gamma2
    xi_out[count, ] <- xi
    if(sample.N){
      N_out[, count] <- N
    }
    if (count%%50 == 0){
      propSd <- sapply(1:length(theta), function(i){
        ad.propSd(theta_out[(count-50+1):count, i], propSd[i], count/50) })
#      print(propSd)
    }

  }
  if(sample.N){
    out <- list(gamma2 = gamma2_out, theta = theta_out, N = N_out, xi = xi_out)
  }else{
    out <- list(gamma2 = gamma2_out, theta = theta_out, xi = xi_out)
  }
  return(out)
}

