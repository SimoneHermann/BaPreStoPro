

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

#' SMC
#'
#' @description Sequential Monte Carlo step - conditional and unconditional.
#' @param phi parameter \eqn{\phi}
#' @param gamma2 parameter \eqn{\gamma^2}
#' @param sigma2 parameter \eqn{\sigma^2}
#' @param Npart number of particle
#' @param times vector of time points
#' @param y observation vector
#' @param b.fun drift function
#' @param y0.fun function for starting point dependent on \eqn{\phi}
#' @param sigmaTilde variance function
#' @param conditional logical(1), if TRUE conditional SMC
#' @param X.cond if conditional = TRUE, the series to be hold
#' @param B_fixed ancestral lineage of X.cond

#' @author
#' Simone Hermann and Adeline Leclercq-Samson

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
      parents[, n-1] <- sample(1:Npart, Npart, replace = TRUE, prob = Weights[,n-1])
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

