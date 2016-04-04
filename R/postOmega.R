
#' Posterior
#'
#' @description Posterior for parameters \eqn{\Omega}
#' @param alpha vector of prior variables
#' @param beta vector of prior variables
#' @param phi matrix of random effects
#' @param mu mean of random effects
#' @return one sample of posterior
#' @export


postOmega <- function(alpha, beta, phi, mu){  # length(alpha)=length(beta)=length(mu)
  p <- length(mu)
  Dia <- numeric(p)
  for(i in 1:p){
    Dia[i] <- 1/rgamma(1, alpha[i] + nrow(phi)/2, beta[i] + sum((phi[,i]-mu[i])^2)/2)
  }
  Dia
}

#' Posterior
#'
#' @description Posterior for parameters \eqn{\Omega}
#' @param R prior matrix of wishart distribution
#' @param phi matrix of random effects
#' @param mu mean of random effects
#' @return one sample of posterior

postOmega_matrix <- function(R, phi, mu){
  Rpost <- solve(R + (t(phi)-as.vector(mu))%*%t((t(phi)-as.vector(mu))))
  solve( rWishart(1,nrow(phi)+length(mu)+1,Rpost)[,,1])
}
