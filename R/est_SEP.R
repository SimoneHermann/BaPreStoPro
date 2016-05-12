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
#' @param propSd proposal standard deviation
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
#' chain <- est_SEP(jumpTimes, start = c(1, 0.25), s = s, n = 1000)
#' @export
est_SEP <- function(jumpTimes, Tend, start, n = 5000, s = 200, f_xi, priorRatio, proposals, Lmax = 35, propSd){

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

  if(missing(propSd)) propSd <- abs(start)/10 + 0.001
  if(missing(proposals)){
    proposals <- list()
    proposals$draw <- function(xi_old, propSd){
      proposal(xi_old, propSd)
    }
    proposals$ratio <- function(xi_drawn, xi_old, propSd){
      proposalRatio(xi_old, xi_drawn, propSd)
    }
  }

  xi_old <- start
  xi_out <- matrix(0, n, length(start))
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
    xi_out[count,] <- xi_old
  }


  xi_out
}

#' Prediction of SEP
#'
#' @description Bayesian prediction of the self-exciting process.
#' @param xi Markov chains
#' @param s stress
#' @param J index
#' @param f_xi_mat matrix-vise function
#' @param len length of samples
#' @param end how many jumps
#' @param grid fineness degree
#' @param Lmax maximal number of tension wires that can break
#' 
#' @examples
#' s <- seq(100, 300, length = 10)
#' f_xi <- function(y, xi) xi[1] - xi[2]*y
#' jumpTimes <- lapply(1:10, function(a){
#' wt <- rexp(20, exp(- f_xi(s[a]/(35-(0:19)), c(1, 0.25))))
#' cumsum(wt)
#' })
#' \dontrun{
#' chain <- est_SEP(jumpTimes, start = c(1, 0.25), s = s, n = 11000)
#' pred <- pred_SEP(chain[seq(1001, 11000, by = 10),], s = s[1], J = 1, f_xi_mat = function(y, xi) xi[,1] - xi[,2]*y)
#' qu <- apply(pred, 2, quantile, c(0.025, 0.975))
#' plot(jumpTimes[[1]], ylim = range(c(jumpTimes, range(qu))))
#' lines(qu[1,])
#' lines(qu[2,])
#' }
#' @export

pred_SEP <- function(xi, s, J = 1, f_xi_mat, len = 100, end = 20, grid = 0.01, Lmax = 35){

  lambda <- function(j){
    exp(-f_xi_mat(s/(Lmax-j-1), xi)) 
  }
  Fun <- function(cand, j){
    1-mean(exp(-lambda(j)*cand))
  }
  
  drawT <- function(j){
    u <- runif(1)
    cand <- 1
    memory <- cand
    
    while(length(unique(memory)) == length(memory)){
      if(Fun(cand, j) < u){
        cand <- cand*2
        memory <- c(memory,cand)
      }else{
        cand <- cand/2
        memory <- c(memory,cand)
      }
    }  
    lower <- memory[length(memory)]
    upper <- memory[length(memory)-1]
    diff <- upper - lower
    while(diff >= grid){
      if(Fun(lower+diff/2,j) < u){
        lower <- lower+diff/2
      }else{
        upper <- lower+diff/2
      }
      diff <- upper - lower
    }  
    (lower+upper)/2
  }
  
  result <- sapply(J:end, function(k){
    replicate(len, drawT(k))
  })
  t(apply(result, 1, cumsum))
}



