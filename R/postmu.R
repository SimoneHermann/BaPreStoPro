#' Posterior for mu
#'
#' @description Posterior for parameters \eqn{\mu}
#' @param phi matrix of random effects
#' @param m prior mean
#' @param v prior variance
#' @param Omega variance of the random effects
#' @return one sample of posterior
#' @export



# postmu <- function(phi, m, v, Omega){  # phi matrix
#   n <- nrow(phi)
#   V_ <- diag(1/v)
#   Vpost <- solve(V_ + n*solve(Omega))
#   mpost <- Vpost%*%( V_%*%m + apply((solve(Omega)%*%t(phi)),1,sum) )
#
#   rmvnorm(1,mpost,Vpost)
# }
# in the case of Omega = vector... possibly faster...
postmu <- function(phi, m, v, Omega){  # phi nxp-matrix, m mean of mu, v diagonal of variance matrix
  n <- nrow(phi)
  V_ <- 1/v
  Vpost <- 1/(V_ + n/Omega)
  mpost <- Vpost*( m/v + apply(t(phi)/Omega, 1, sum) )

  rnorm(length(mpost), mpost, sqrt(Vpost))
}
