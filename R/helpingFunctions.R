

#' Transformation of vector of event times to counting process
#'
#' @description Transformation of vector of event times to counting process.
#' @param times vector of event times
#' @param t times of counting process
#'
#' @return vector of counting process observations in t
#' @examples
#' t <- seq(0, 1, by = 0.01)
#' times <- simN(t, c(5, 0.5), len = 1)$Times
#' process <- TimestoN(times, t)
#' @export
TimestoN <- function(times, t){
  lt <- length(t)
  dN <- numeric(lt)
  ind <- unlist(sapply(times[times <= max(t)], function(a){which(abs(t-a) == min(abs(t-a)))}))
  for(a in 1:length(ind)) dN[ind[a]] <- dN[ind[a]] + 1
  cumsum(dN)
}


#' Transformation of counting process to vector of event times
#'
#' @description Transformation of vector of counting process to event times.
#' @param dN vector of differences of counting process
#' @param t times of counting process
#'
#' @return vector of event times
#' @examples
#' t <- seq(0, 1, by = 0.01)
#' process <- simN(t, c(5, 0.5), len = 1)$N
#' times <- dNtoTimes(diff(process), t)
#' @export
dNtoTimes <- function(dN, t){
  if(any(dN > 1)){
    m <- sum(dN > 0)
    res <- NULL
    he1 <- dN[dN > 0]
    he2 <- t[dN > 0]
    for(mi in 1:m){
      res <- c(res, rep(he2[mi], he1[mi]))
    }
  }else{
    res <- t[dN > 0]
  }
  res
}



#' Sampling from proposal density for strictly positive parameters
#'
#' @description Used in Metropolis Hastings algorithms.
#' @param parOld the parameter from the last iteration step
#' @param propSd proposal standard deviation
#'
#' @return candidate for MH ratio
#' @examples
#' plot(replicate(100, proposal(1, 0.1)), type = "l")
#' @export
proposal <- function(parOld, propSd){
  if(any(parOld < 1e-150)) parOld[parOld < 1e-150] <- 1e-150  # 1e-320 equal to zero ...
  mu <- log(parOld) - log( propSd^2/(parOld^2) + 1)/2
  sigma2 <- log( propSd^2/(parOld^2)+1)
  rlnorm(length(parOld), mu, sqrt(sigma2))
}

#' Sampling from proposal density for strictly positive parameters
#'
#' @description Used in Metropolis Hastings algorithms.
#' @param parOld the parameter from the last iteration step
#' @param parNew drawn candidate
#' @param propSd proposal standard deviation
#'
#' @return MH ratio
#' @examples
#' cand <- proposal(1, 0.01)
#' proposalRatio(1, cand, 0.01)
#' @export
proposalRatio <- function(parOld, parNew, propSd){
  muOld <- log(parOld) - log( propSd^2/exp(2*log(parOld)) + 1)/2
  sigma2Old <- log( propSd^2/exp(2*log(parOld))+1)
  muNew <- log(parNew) - log( propSd^2/exp(2*log(parNew)) + 1)/2
  sigma2New <- log( propSd^2/exp(2*log(parNew))+1)

  prod(dlnorm(parOld, muNew, sqrt(sigma2New))/dlnorm(parNew, muOld, sqrt(sigma2Old)))
}

#' Plot function for credibility or prediction intervals
#'
#' @description Plots intervals.
#' @param input list or matrix of samples from posterior or predictive distribution
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main an overall title for the plot
#' @param l lower bound
#' @param u upper bound
#' @param color color of the lines to be drawn
#'
plotQuantiles <- function(input, xlab = "", ylab = "", main = "", l = 0.025, u = 0.975, color = 1){
  if(is.list(input)){
    qu1 <- sapply(input, quantile, l)
    qu2 <- sapply(input, quantile, u)
    me <- sapply(input, median)
  }else{
    if(is.matrix(input)){
      qu1 <- apply(input, 2, quantile, l)
      qu2 <- apply(input, 2, quantile, u)
      me <- apply(input, 2, median)
    }else{print("input has to be a list or a matrix")}
  }
  len <- length(me)
  ra <- range(c(qu1, qu2))
  if(any(qu2[1:20] > (ra[2]-(ra[2]-ra[1])*0.1))){
    yplot <- ra + c(0, (ra[2]-ra[1])*0.1)
  }else{
    yplot <- ra
  }
  plot(me, pch = 20, ylim = yplot, xlab = xlab, ylab = ylab, main = main)
  segments(1:len, qu1, 1:len, qu2, col = color)
  segments(1:len-len/200, qu1, 1:len+len/200, qu1, col = color)
  segments(1:len-len/200, qu2, 1:len+len/200, qu2, col = color)
  legend("topleft", "median", pch = 20, inset = 0.01, box.lty = 0)
}

#' Calculation of interval score
#'
#' @description Scoring rule of Raftery and Gneiting (??).
#' @param l lower bound
#' @param u upper bound
#' @param x true value
#' @param alpha level
#'
#' @return interval score
#' @export
scoreRule <- function(l, u, x, alpha = 0.05){
  u-l + 2/alpha*(l-x)*(x < l) + 2/alpha*(x-u)*(x > u)
}


#' Binary Search Algorithm
#'
#' @description Binary Search Algorithm
#' @param Fun cumulative distribution function
#' @param len number of samples
#' @param candArea candidate area
#' @param grid fineness degree
#' @param method vectorial ("vector") or not ("free")
#'
#' @return vector of samples
#' @examples
#' test <- BinSearch(function(x) pnorm(x, 5, 1), 1000, candArea = c(0, 10), method = "free")
#' plot(density(test))
#' curve(dnorm(x, 5, 1), col = 2, add = TRUE)
#' @export
BinSearch <- function(Fun, len, candArea, grid = 1e-05, method = c("vector", "free")){
  method <- match.arg(method)
  if(missing(candArea)){
    candArea <- findCandidateArea(Fun)
    if(method == "vector") candArea <- seq(candArea[1], candArea[2], by = grid)
  }else{
    if(method == "vector" & length(candArea) == 2) candArea <- seq(candArea[1], candArea[2], by = grid)
    if(method == "free" & length(candArea) > 2) candArea <- c(min(candArea), max(candArea))
  }

  if(method == "vector"){
    diFu <- sapply(candArea, Fun)
    U <- runif(len, 0, max(diFu))
    res <- sapply(U, function(u) candArea[which(diFu >= u)[1]])
  }
  if(method == "free"){
    res <- numeric(len)
    U <- runif(len)
    for(i in 1:len){
      lower <- candArea[1]
      upper <- candArea[2]

      diff <- upper - lower
      while(diff >= grid){
        if(Fun(lower+diff/2) < U[i]){
          lower <- lower+diff/2
        }else{
          upper <- lower+diff/2
        }
        diff <- upper - lower
      }
      res[i] <- (lower+upper)/2
    }
  }
  res
}

#' Rejection Sampling Algorithm
#'
#' @description Rejection Sampling
#' @param Fun cumulative distribution function
#' @param dens density
#' @param len number of samples
#' @param cand candidate area
#' @param grid fineness degree
#' @param method vectorial ("vector") or not ("free")
#'
#' @return vector of samples
#' @examples
#' plot(density(RejSampling(Fun = function(x) pnorm(x, 5, 1), dens = function(x) dnorm(x, 5, 1), len = 500, cand = seq(2, 9, by = 0.001), method = "free")))
#' lines(density(RejSampling(function(x) pnorm(x, 5, 1), function(x) dnorm(x, 5, 1), 500, cand = seq(2, 9, by = 0.001), method = "vector")), col=2)
#' curve(dnorm(x, 5, 1), from = 2, to = 8, add = TRUE, col = 3)
#' @export
RejSampling <- function(Fun, dens, len, cand, grid = 1e-03, method = c("vector", "free")){  # for negative support?!?
  method <- match.arg(method)
  if(method == "free"){
    res <- numeric(len)
    for(i in 1:len){
      if(missing(cand)){
        ca <- findCandidateArea(function(t) Fun(t))
      }else{
        ca <- range(cand)
      }
      mp <- optimize(f = function(t) dens(t), ca, maximum = TRUE)$objective
      resi <- NULL
      while(is.null(resi)){
        u <- runif(1,0,mp)
        candi <- runif(1, ca[1], ca[2])
        prob <- dens(candi)
        if(u <= prob){
          resi <- candi
        }
      }
      res[i] <- resi
    }
  }
  if(method == "vector"){
    res <- numeric(len)
    if(missing(cand)){
      ca <- findCandidateArea(Fun)
      cand <- seq(ca[1], ca[2], by = grid)
    }
    prob <- vapply(cand, function(v) dens(v), FUN.VALUE = numeric(1))
    cand <- cand[prob != 0]
    prob <- prob[prob != 0]
    mp <- max(prob)
    count <- 1
    while(count <= len){
      u <- runif(1, 0, mp)
      ind <- sample(1:length(cand), 1)
      if(u <= prob[ind]){
        res[count] <- cand[ind]
        count <- count + 1
      }
    }
  }
  return(res)
}


#' Helping function
#'
#' @description Adaptive MCMC
#' @param chain Markov chain
#' @param propSd current proposal standard deviation
#' @param iteration current iteration (batch)
#' @param lower lower bound
#' @param upper upper bound
#' @param delta.n function of batch number
#'
#' @return adjusted proposal standard deviation
#' @export
ad.propSd <- function(chain, propSd, iteration, lower = 0.3, upper = 0.6, delta.n = function(n) min(0.05, 1/sqrt(n))){
  ar <- length(unique(chain))/length(chain)
  lsi <- log(propSd)

  lsi[ar < lower] <- lsi - delta.n(iteration)
  lsi[ar > upper] <- lsi + delta.n(iteration)
  exp(lsi)
}


#' Helping function
#'
#' @description Finding suitable candidate area
#' @param VFun cumulative distribution function
#' @param start starting point for search
#' @param pos.support if TRUE: only positive support
#' @param quasi.null size of values to be defined as zero, default: 10^{-5}
#'
#' @return adjusted proposal standard deviation
#' @export

findCandidateArea <- function(VFun, start = 1, pos.support = TRUE, quasi.null = 1e-05){
  if(pos.support){
    cand <- start
    while(1 - VFun(cand) > quasi.null){
      cand <- cand*2
    }
    upper <- cand

    cand <- start
    while(VFun(cand) > 0){
      cand <- cand/2
    }
    lower <- cand
  } else{
    cand <- start
    while(1 - VFun(cand) > quasi.null){
      cand <- cand*2
    }
    upper <- cand

    cand <- -start
    while(VFun(cand) > quasi.null){
      cand <- cand*2
    }
    lower <- cand
  }

  c(lower, upper)
}


#' Calcucation of burn-in phase and thin rate
#'
#' @description Proposal for burn-in and thin rate
#' @param chain vector of Markov chain samples
#' @param dependence allowed dependence for the chain
#' @param m number of blocks
#' @export
diagnostic <- function(chain, dependence = 0.8, m = 10) {
  lc <- length(chain)
  K <- floor(lc/m)
  thinning <- min(which(acf(chain[-(1:floor(lc/5))], plot=F)$acf <= dependence)[1], floor(K/10), na.rm = TRUE)
  he1 <- sapply(1:m, function(i) chain[((i - 1) * K + 1):(i * K)][seq(1, K, by = thinning)])

  he2 <- apply(he1, 2, quantile, c(0.025, 0.975))
  he.mean <- apply(he1, 2, mean)
  is.in <- (he.mean[-m] >= he2[1, -1] & he.mean[-m] <= he2[2, -1]) | (he.mean[-1] >= he2[1, -m] & he.mean[-1] <= he2[2, -m])
  #    burnIn <- 0
  burnIn <- K
  for (i in 1:(m-1)) {
    if (sum(is.in) < length(is.in)) {
      is.in <- is.in[-1]
      burnIn <- burnIn + K
    }
  }
  burnIn <- min((m-1)*K, burnIn)
  return(c(burnIn = burnIn, thinning = thinning))
}
