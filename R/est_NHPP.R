
#' Metropolis-Hastings sampler
#'
#' @description Bayesian estimation of the parameter of the intensity rate.
#' @param jumpTimes vector of jump (or event) times
#' @param Tend last observation vector of the NHPP
#' @param start starting value
#' @param n length of Markov chain
#' @param int one out of "Weibull" or "Exp"
#' @param priorRatio function for prior ratio, if missing: non-informative
#' @param proposal "lognormal" (for positive parameters, default) or "normal"
#' @param Lambda intensity rate function
#'
#' @return p x n dimensional matrix of posterior samples
#' @examples
#' # exponential
#' Lambda <- function(t, xi){
#' xi[2]*exp(xi[1]*t)-xi[2]
#' }
#' jumpTimes <- simN(seq(0, 1, by = 0.001), c(4, 0.5), len = 1, Lambda = Lambda)$Times
#' chain <- est_NHPP(jumpTimes, 1, c(4, 0.5), Lambda = Lambda) 
#' # chain <- est_NHPP(jumpTimes, 1, c(4, 0.5), int = "Exp")  
#' plot(chain[1,], type="l"); abline(h = 4)
#' plot(chain[2,], type="l"); abline(h = 0.5)
#' 
#' # weibull
#' Lambda <- function(t, xi){
#' (t/xi[2])^xi[1]
#' }
#' jumpTimes <- simN(seq(0, 1, by = 0.001), c(3, 0.5), len = 1, Lambda = Lambda)$Times
#' chain <- est_NHPP(jumpTimes, 1, c(3, 0.5), Lambda = Lambda)   
#' # chain <- est_NHPP(jumpTimes, 1, c(3, 0.5), int = "Weibull")  
#' plot(chain[1,], type="l"); abline(h = 3)
#' plot(chain[2,], type="l"); abline(h = 0.5)
#' 
#' @export
est_NHPP <- function(jumpTimes, Tend, start, n = 5000, int = c("Weibull","Exp"), priorRatio, proposal = c("lognormal", "normal"), Lambda){
  proposal <- match.arg(proposal)

  if(missing(priorRatio)){
    priorRatio <- function(xi_drawn, xi_old) 1
  }
  if(missing(Lambda)){
    int <- match.arg(int)
    if(int=="Weibull"){
      Lambda <- function(t, xi){
        (t/xi[2])^xi[1]
      }
      lambda <- function(t, xi){
        xi[1]/xi[2]*(t/xi[2])^(xi[1]-1)
      }
      Lik <- function(xi){
        lambda_vec <- lambda(jumpTimes, xi)
        prod(lambda_vec)*exp(-Lambda(Tend, xi))
      }
    }else{ # int == "Exp
      Lambda <- function(t, xi){
        xi[2]*exp(xi[1]*t)-xi[2]
      }
      lambda <- function(t, xi){
        xi[1]*xi[2]*exp(xi[1]*t)
      }
      Lik <- function(xi){
        lambda_vec <- lambda(jumpTimes, xi)
        prod(lambda_vec)*exp(-Lambda(Tend, xi))
      }
    }
    proposal <- "lognormal"
    
  } else{
    lambda <- function(t, xi, h = 1e-05){
      (Lambda(t+h,xi)-Lambda(t,xi))/h
    }
    Lik <- function(xi){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
  }
  if(proposal == "lognormal"){
    proposals <- list()
    proposals$draw <- function(xi_old, propSd){
      proposal(xi_old, propSd)
    }
    proposals$ratio <- function(xi_drawn, xi_old, propSd){
      proposalRatio(xi_old, xi_drawn, propSd)
    }
  }else{
    proposals <- list()
    proposals$draw <- function(xi_old, propSd){ 
      rnorm(length(xi_old), xi_old, propSd)
    }
    proposals$ratio <- function(xi_drawn, xi_old, propSd) 1
  }
  propSd <- (abs(start)+0.1)/2
  xi_old <- start
  xi_out <- matrix(0, length(start), n)
  LikOld <- Lik(xi_old)

  for(count in 1:n){
    xi_drawn <- proposals$draw(xi_old, propSd)

    LikNew <- Lik(xi_drawn)
    ratio <- proposals$ratio(xi_drawn, xi_old, propSd)
    ratio <- ratio*priorRatio(xi_drawn, xi_old)
    ratio <- ratio*LikNew/LikOld
    if(is.na(ratio)) ratio <- 0

    if(runif(1) <= ratio){
      xi_old <- xi_drawn
      LikOld <- LikNew
    }
    xi_out[,count] <- xi_old

    if (count%%50 == 0){
      propSd <- sapply(1:length(start), function(i){
        ad.propSd(xi_out[i, (count-50+1):count], propSd[i], count/50) })
    }

  }
  xi_out
}
