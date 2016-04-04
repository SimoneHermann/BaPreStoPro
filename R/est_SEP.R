#' Metropolis-Hastings sampler
#'
#' @description Bayesian estimation of the parameter of the self-exciting process.
#' @param jumpTimes vector of event times
#' @param Tend last observation time, if missing => last event time
#' @param start starting value
#' @param n length of Markov chain
#' @param s stress range
#' @param f_xi exp(-f_xi) intensity function
#' @param priorRatio function(old, new) for prior ratio, if missing => noninformative
#' @param proposals list of functions: draw(old) and ratio(new, old)
#' @param Lmax maximal number of tension wires that can break
#'
#' @return p x n dimensional matrix of posterior samples, p length of vector xi
#' @examples
#' wt <- rexp(20, exp(- 1 + 0.25*200/(35-(0:19))))
#' jumpTimes <- cumsum(wt)
#' chain <- est_SEP(jumpTimes, start = c(1, 0.25), s = 200)
#' par(mfrow = c(2,1))
#' plot(chain[1,], type = "l")
#' abline(h = 1, col = 2)
#' plot(chain[2,], type = "l")
#' abline(h = 0.25, col = 2)
#'
#' print("acceptance rate:")
#' length(unique(chain[1,]))/length(chain[1,])
#'
#' # or:
#' s <- seq(100, 300, length = 10)
#' jumpTimes <- lapply(1:10, function(a){
#' wt <- rexp(20, exp(- 1 + 0.25*s[a]/(35-(0:19))))
#' cumsum(wt)
#' })
#' chain <- est_SEP(jumpTimes, start = c(1, 0.25), s = s)
#' @export
est_SEP <- function(jumpTimes, Tend, start, n = 5000, s = 200, f_xi, priorRatio, proposals, Lmax = 35){

  if(is.list(jumpTimes)){
    N <- length(jumpTimes)
    m <- sapply(jumpTimes, length)
    if(missing(Tend)) Tend <- sapply(1:N, function(i) jumpTimes[[i]][m[i]])
    jumps <- lapply(1:N, function(i) c(0, jumpTimes[[i]], Tend[i]))
    aj <- lapply(1:N, function(i) s[i]/(Lmax-0:m[i]))
    Lambda <- function(xi){
      he <- sapply(1:N, function(i) sum(exp(-f_xi(aj[[i]], xi))*diff(jumps[[i]])))
      sum(he)
    }
    Lik <- function(xi){
      he <- sapply(1:N, function(i) exp(-f_xi(aj[[i]][-(m[i]+1)], xi)))
      prod(unlist(he))*exp(-Lambda(xi))
    }

  }else{
    m <- length(jumpTimes)
    if(missing(Tend)) Tend <- jumpTimes[m]
    jumps <- c(0, jumpTimes, Tend)
    aj <- s/(Lmax-0:m)
    Lambda <- function(xi){
      sum(exp(-f_xi(aj, xi))*diff(jumps))
    }
    lambda <- function(xi){
      exp(-f_xi(aj[-(m+1)], xi))
    }
    Lik <- function(xi){
      prod(lambda(xi))*exp(-Lambda(xi))
    }
  }


  if(missing(f_xi)) f_xi <- function(aj, xi) xi[1] - xi[2]*aj
  if(!(missing(proposals) || all(c("ratio","draw")%in%names(proposals)))){
    warning("specify a ratio function and a draw function")
  }
  if(missing(priorRatio)){
    priorRatio <- function(xi_drawn, xi_old) 1
  }


  if(missing(proposals)){
    proposal <- function(parOld, propSd){
      help <- rnorm(length(parOld), sqrt(parOld), propSd)
      help^2
    }
    proposalRatio <- function(parOld, parNew){
      prod(sqrt(parNew)/sqrt(parOld))
    }

    proposals <- list()
    proposals$draw <- function(xi_old){
      proposal(xi_old, (abs(start)+0.1)/10)
    }
    proposals$ratio <- function(xi_drawn, xi_old){
      proposalRatio(xi_old, xi_drawn)
    }
  }

  xi_old <- start
  xi_out <- matrix(0, length(start), n)
  LikOld <- Lik(xi_old)

  for(count in 1:n){
    xi_drawn <- proposals$draw(xi_old)

    LikNew <- Lik(xi_drawn)
    ratio <- proposals$ratio(xi_drawn, xi_old)
    ratio <- ratio*priorRatio(xi_drawn, xi_old)
    ratio <- ratio*LikNew/LikOld
    if(is.na(ratio)) ratio <- 0

    if(runif(1) <= ratio){
      xi_old <- xi_drawn
      LikOld <- LikNew
    }
    xi_out[,count] <- xi_old
  }


  xi_out
}

