


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
        exp(xi[1]*t+xi[2])-exp(xi[2])
      }
    }
  }

  if(missing(N)){
    post_N <- function(N, i, gamma2, theta, xi){
      Di <- Y[i]-fun(t[i], N, theta)
      Lambdai <- Lambda(t[i], xi)
      exp(-Lambdai)/prod(1:max(N,1))*exp(-Di^2/(2*gamma2)+N*log(Lambdai))/sqrt(2*pi*gamma2)
    }
    drawN_wrongbutworking <- function(gamma2, theta, xi, N_old){   # wrong but working...
      N_new <- N_old
      cands <- lapply(N_old, function(n) 0:(n*rangeN+5))
      prob <- lapply(1:lt, function(i) sapply(cands[[i]], post_N, i, gamma2, theta, xi))
      diFu <- lapply(prob, cumsum)
      ind <- sapply(diFu, function(vec) any(is.na(vec)) | any(is.infinite(vec)) )
      u <- numeric(length(diFu))
      if(sum(ind) < length(diFu)){
        u[!ind] <- runif(length(diFu[!ind]), 0, sapply(diFu[!ind], max))
        for(j in (1:lt)[!ind]){
                  N_new[j] <- cands[[j]][which(diFu[[j]]>=u[j])[1]]
        }
      }
        N_new
    }
    drawN2 <- function(gamma2, theta, xi, N_old){
      N_new <- N_old
      dN_old <- diff(N_old)
      for(i in 2:lt){
        cands <- 0:(dN_old[i-1]*rangeN+5) + N_new[i-1]
        prob <- post_N(cands, i, gamma2, theta, xi)
        diFu <- cumsum(prob)
        if(!(any(is.na(diFu)) | any(is.infinite(diFu)))){
          u <- runif(1, 0, max(diFu))
          N_new[i] <- cands[which(diFu>=u)[1]]
        }else{
          N_new[i] <- N_new[i-1] + dN_old[i]
        }
      }
      N_new
    }
    dY <- diff(Y)
    drawN <- function(gamma2, theta, xi, N_old){
      N_new <- N_old
      dN_old <- diff(N_old)
      for(i in 2:lt){
        cands <- 0:(dN_old[i-1]*rangeN+5)
        prob <- dnorm(dY[i-1], fun(t[i],N_new[i-1]+cands,theta)-fun(t[i-1],N_new[i-1],theta), sqrt(2*gamma2))*
          dpois(cands, Lambda(t[i],xi)-Lambda(t[i-1],xi))
        diFu <- cumsum(prob)
        if(!(any(is.na(diFu)) | any(is.infinite(diFu)))){
          u <- runif(1, 0, max(diFu))
          N_new[i] <- N_new[i-1] + cands[which(diFu>=u)[1]]
        }else{
          N_new[i] <- N_new[i-1] + dN_old[i]
        }
      }
      N_new
    }

    Npart <- 10
    B_fixed <- rep(1, lt)
    A.to.B = function(A){
      B <- matrix(0, Npart, lt)
      B[,lt] <- 1:Npart;
      for (t in seq(lt, 2, by = -1)){
        B[,t-1] <- A[B[,t], t-1]}
      return(B)
    }

    conditional_SMC = function(theta, gamma2, xi, N.cond, B_fixed){# conditional SMC
      # N.cond = the old samples

      x <- matrix(0, Npart, lt)
      w <- matrix(1, Npart, lt)
      W <- matrix(1, Npart, lt)
      parents <-  matrix(1, Npart, lt-1)

      # initialisation
      x[,1] <- 0  # poisson always starts in 0
      x[B_fixed[1],1] <- N.cond[1]
      w[,1] <- dnorm(Y[1], mean = fun(t[1], x[,1], theta), sd = sqrt(gamma2))
      W[,1] <- w[,1]/sum(w[,1])

      for (indtime in 2:lt){
        v  <- 1:Npart; v <- v[v != B_fixed[indtime]]
        parents[B_fixed[indtime], indtime-1] <- B_fixed[indtime-1]
        O <- sample(1:Npart, 1, prob = dbinom(1:Npart, Npart, prob = W[B_fixed[indtime], indtime-1]));
        r <- NULL
        if (O > 1){
          r <- sample(v, O-1, replace = FALSE);
          parents[r, indtime-1] <- B_fixed[indtime-1];
        }
        if(O != Npart){
          u <- 1:Npart; u <- u[-c(r, B_fixed[indtime])]
          parents[u, indtime-1] <-  sample(v, Npart-O, replace = TRUE, prob = W[-B_fixed[indtime], indtime-1])
        }

        x[, 1:(indtime-1)] <- x[parents[,indtime-1], 1:(indtime-1)]
        x.past <- x[, indtime-1]

        dr <- function(Ni, dNi_old){
          cands <- 0:(dNi_old*rangeN+5)
          prob <- dnorm(Y[indtime], fun(t[indtime], Ni+cands,theta), sqrt(2*gamma2))*
            dpois(Ni+cands, Lambda(t[indtime],xi))
          diFu <- cumsum(prob)
          u <- runif(1, 0, max(diFu))
          Ni + cands[which(diFu >= u)[1]]
        }
        x.new <- sapply(x.past, dr, diff(N.cond)[indtime-1])

        x[v,indtime] <- x.new[v]
        x[-v,indtime] <- N.cond[indtime]
        wn <- dpois(x[,indtime],Lambda(t[indtime],xi))*
          dnorm(Y[indtime], fun(t[indtime], x[,indtime], theta),sqrt(gamma2))/(dnorm(Y[indtime], fun(t[indtime], x[,indtime], theta), sqrt(2*gamma2))*
                                                                              dpois(x[,indtime], Lambda(t[indtime],xi)))

        w[,indtime] <- wn
        W[,indtime] <- w[,indtime]/sum(w[,indtime])
      }
      lign.B <- A.to.B(parents)
      indice <- sample(1:Npart, 1, prob = W[,lt])
      X <- x[indice,]
      B_fixed <- lign.B[indice,]
      return(list(N = X, B_fixed = B_fixed) )
    }

    if(is.null(start$N)){
      N <- simN(t, xi, 1, start = c(t[1], 0), Lambda = Lambda)$N
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
#      N <- drawN_wrongbutworking(gamma2, theta, xi, N)
      he <- conditional_SMC(theta, gamma2, xi, N, B_fixed)
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

