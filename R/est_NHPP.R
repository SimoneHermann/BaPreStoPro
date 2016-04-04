
#' Metropolis-Hastings sampler
#'
#' @description Bayesian estimation of the parameter of the intensity rate.
#' @param jumpTimes vector of jump (or event) times
#' @param Tend last observation vector of the NHPP
#' @param start starting value
#' @param n length of Markov chain
#' @param int one out of "Weibull" or "Exp"
#' @param priorRatio function for prior ratio, if missing: non-informative
#' @param proposals list of functions "ratio" and "draw"
#' @param Lambda intensity rate function
#'
#' @return p x n dimensional matrix of posterior samples
#' @examples
#' jumpTimes <- simN(seq(0, 1, by = 0.01), c(4, -1), len = 1, int = "Exp")$Times
#' chain <- est_NHPP(jumpTimes, 1, c(4, -1), int = "Exp")
#' par(mfrow = c(2, 1))
#' plot(chain[1,], type="l")
#' abline(h = 4, col = 2)
#' plot(chain[2,], type = "l")
#' abline(h = -1, col = 2)
#'
#' print("acceptance rate:")
#' length(unique(chain[1,]))/length(chain[1,])
#' @export
est_NHPP <- function(jumpTimes, Tend, start, n = 5000, int = c("Weibull","Exp"), priorRatio, proposals, Lambda){

  if(!(missing(proposals) || all(c("ratio","draw")%in%names(proposals)))){
    warning("specify a ratio function and a draw function")
  }
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
    if(missing(proposals)){
      proposal <- function(parOld, propSd){
        help <- rnorm(length(parOld), sqrt(parOld), propSd)
        help^2
      }
      proposalRatio <- function(parOld, parNew, propSd){
        prod(sqrt(parNew)/sqrt(parOld))
      }

      propSd <- (abs(start)+0.1)/2

      proposals <- list()
      proposals$draw <- function(xi_old, propSd){
        proposal(xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd){
        proposalRatio(xi_old, xi_drawn, propSd)
      }
    }
    }else{ # int == "Exp

      Lambda <- function(t, xi){
        exp(xi[1]*t+xi[2])-exp(xi[2])
      }
      lambda <- function(t, xi){
        xi[1]*exp(xi[1]*t+xi[2])
      }
      Lik <- function(xi){
        lambda_vec <- lambda(jumpTimes, xi)
        prod(lambda_vec)*exp(-Lambda(Tend, xi))
      }
      if(missing(proposals)){
        proposal <- function(parOld, propSd){
          help <- rnorm(length(parOld), sqrt(parOld), propSd)
          help^2
        }
        proposalRatio <- function(parOld, parNew, propSd){
          prod(sqrt(parNew)/sqrt(parOld))
        }

        propSd <- (abs(start)+0.1)/10
        proposals <- list()
        proposals$draw <- function(xi_old, propSd){
          c(proposal(xi_old[1], propSd[1]), rnorm(1, xi_old[2], propSd[2]))
        }
        proposals$ratio <- function(xi_drawn, xi_old, propSd){
          proposalRatio(xi_old[1], xi_drawn[1], propSd[1])
        }
      }
    } # end else
  } else{
    lambda <- function(t, xi, h = 1e-05){
      (Lambda(t+h,xi)-Lambda(t,xi))/h
    }
    Lik <- function(xi){
      lambda_vec <- lambda(jumpTimes, xi)
      prod(lambda_vec)*exp(-Lambda(Tend, xi))
    }
    if(missing(proposals)){
      propSd <- (abs(start)+0.1)/10
      proposals <- list()
#       proposals$draw <- function(xi_old, propSd){  # change ? only for positive paramters...
#         proposal(xi_old, propSd)
#       }
#       proposals$ratio <- function(xi_drawn, xi_old, propSd) proposalRatio(xi_old, xi_drawn, propSd)
      proposals$draw <- function(xi_old, propSd){  # change ? only for positive paramters...
        rnorm(length(xi_old), xi_old, propSd)
      }
      proposals$ratio <- function(xi_drawn, xi_old, propSd) 1

    }
  }

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
