

########
#' Prediction for diffusion process
#'
#' @description Bayesian prediction of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param object class object of MCMC samples: "est.Diffusion"
#' @param t vector of time points to make predictions for
#' @param Euler.interval if TRUE: simple prediction intervals with Euler are made (in one step each)
#' @param level level of the prediction intervals
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param b.fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param y.start optional, if missing, first (which.series = "new") or last observation variable ("current") is taken
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param method vectorial ("vector") or not ("free")
#' @param sampling.alg sampling algorithm, inversion method ("BinSearch") or rejection sampling ("RejSamp")
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est_diff <- estimate(model, t, data, 1000) # better: 10000
#' plot(est_diff)
#' \dontrun{
#' pred_diff <- predict(est_diff, t = seq(0, 1, by = 0.1))
#' pred_diff <- predict(est_diff, b.fun.mat = function(phi, t, y) phi[,1])  # much faster
#' pred_diff2 <- predict(est_diff, which.series = "current", b.fun.mat = function(phi, t, y) phi[,1])
#' pred_diff3 <- predict(est_diff, which.series = "current", y.start = data[51], 
#'                t = t[seq(51, 100, by = 5)], b.fun.mat = function(phi, t, y) phi[,1])
#' }
#' pred_diff <- predict(est_diff, Euler.interval = TRUE, b.fun.mat = function(phi, t, y) phi[,1])  
#' # one step Euler approximation
#' pred_diff <- predict(est_diff, pred.alg = "simpleTrajectory", sample.length = 100)
#' for(i in 1:100) lines(t[-1], pred_diff[i,], col = "grey")
#' pred_diff <- predict(est_diff, pred.alg = "simpleBayesTrajectory")
#' @export
setMethod(f = "predict", signature = "est.Diffusion",
          definition = function(object,
                                t, Euler.interval = FALSE, level = 0.05, burnIn, thinning,
                                b.fun.mat, which.series = c("new", "current"), y.start, M2pred = 10,
                                cand.length = 1000, pred.alg = c("Distribution", "Trajectory", "simpleTrajectory", "simpleBayesTrajectory"),
                                method = c("vector", "free"),
                                sampling.alg = c("BinSearch", "RejSamp"), sample.length, grid, plot.prediction = TRUE) {

  pred.alg <- match.arg(pred.alg)
  method <- match.arg(method)
  sampling.alg <- match.arg(sampling.alg)
  which.series <- match.arg(which.series)

  if(missing(y.start)){
    if(which.series == "new") y.start <- object@Y[1]
    if(which.series == "current") y.start <- object@Y[length(object@Y)]
  }

  if(missing(burnIn)) burnIn <- object@burnIn
  if(missing(thinning)) thinning <- object@thinning

  if(missing(t)){
    if(which.series == "new") t <- object@t
    if(which.series == "current"){
      dt <- median( diff(object@t))
      t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
    }
  }

  ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
  samples <- list(phi = as.matrix(object@phi[ind,]), gamma2 = object@gamma2[ind])
  K <- length(samples$gamma2)
  if(missing(sample.length) | pred.alg == "Distribution") sample.length <- K
  b.fun <- object@model$b.fun
  sT.fun <- object@model$sT.fun
  n <- length(t)
  dt <- diff(t)

  if(pred.alg == "simpleTrajectory"){
    if(missing(sample.length)) sample.length <- 100
    result <- matrix(0, sample.length, n)
    phi.est <- apply(samples$phi, 2, mean); gamma2.est <- mean(samples$gamma2)
    cl <- set.to.class("Diffusion", parameter = list(phi = phi.est, gamma2 = gamma2.est), b.fun = b.fun, sT.fun = sT.fun)
    for(i in 1:sample.length){
      result[i,] <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
    }
    result <- result[,-1]
  }
  if(pred.alg == "simpleBayesTrajectory"){
    if(missing(sample.length)){
      sample.length <- K
    } else{
      if(sample.length > K) sample.length <- K
    }
    result <- matrix(0, sample.length, n)
    for(i in 1:sample.length){
      cl <- set.to.class("Diffusion", parameter = list(phi = samples$phi[i,], gamma2 = samples$gamma2[i]), b.fun = b.fun, sT.fun = sT.fun)
      result[i,] <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
    }
    result <- result[,-1]
      
  }
  if(pred.alg == "Distribution"){

      if(Euler.interval){
        result <- matrix(0, 2, n-1)
        for(i in 1:(n-1)){
          cand.Area <- c(y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])),
                         y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])) )
          if(!missing(b.fun.mat)){
            VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1]))))
          }else{
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1])))))
          }
          result[,i] <- prediction.intervals(samples, VFun, x0 = y.start, level = level, candArea = cand.Area)
        }
      }else{
  
        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
          dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
        }else{
          if(pred.alg == "Distribution"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
          }
          if(pred.alg == "Trajectory"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
          }
        }
        cand.Area <- c(y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]),
                       y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]) )
        if(missing(grid)) d <- diff(cand.Area)/cand.length
  
        result <- matrix(0, sample.length, n-1)
        result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
        for(i in 2:(n-1)){
  
          if(!missing(b.fun.mat)){
            VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
            dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
          }else{
            if(pred.alg == "Distribution"){
              VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
              dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
            }
            if(pred.alg == "Trajectory"){
              VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
              dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
            }
          }
  
          cand.Area <- c(min(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]),
                         max(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]) )
          if(missing(grid)) d <- diff(cand.Area)/cand.length
  
          result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
                                  pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)

          if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
        }
      }
    }
  
  if(pred.alg == "Trajectory"){
    
    
    if(Euler.interval){
      result <- matrix(0, 2, n-1)
      for(i in 1:(n-1)){
        cand.Area <- c(y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])),
                       y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])) )
        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1]))))
        }else{
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1])))))
        }
        result[,i] <- prediction.intervals(samples, VFun, x0 = y.start, level = level, candArea = cand.Area)
      }
    }else{
      
      if(!missing(b.fun.mat)){
        VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
        dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
      }else{
        if(pred.alg == "Distribution"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
        }
        if(pred.alg == "Trajectory"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
        }
      }
      cand.Area <- c(y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]),
                     y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]) )
      if(missing(grid)) d <- diff(cand.Area)/cand.length
      
      result <- matrix(0, sample.length, n-1)
      result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
                              pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
      
      for(i in 2:(n-1)){
        
        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
          dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
        }else{
          if(pred.alg == "Distribution"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
          }
          if(pred.alg == "Trajectory"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
          }
        }
        
        cand.Area <- c(min(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]),
                       max(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]) )
        if(missing(grid)) d <- diff(cand.Area)/cand.length
        
        result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
        if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
      }
      
    }
  }
  
#   if(pred.alg %in% c("Distribution", "Trajectory")){
#     
#     
#     if(Euler.interval){
#       result <- matrix(0, 2, n-1)
#       for(i in 1:(n-1)){
#         cand.Area <- c(y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])),
#                        y.start + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*(t[i+1]-t[1])) )
#         if(!missing(b.fun.mat)){
#           VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1]))))
#         }else{
#           VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1])))))
#         }
#         result[,i] <- prediction.intervals(samples, VFun, x0 = y.start, level = level, candArea = cand.Area)
#       }
#     }else{
#       
#       if(!missing(b.fun.mat)){
#         VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
#         dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
#       }else{
#         if(pred.alg == "Distribution"){
#           VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
#           dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
#         }
#         if(pred.alg == "Trajectory"){
#           VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
#           dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
#         }
#       }
#       cand.Area <- c(y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]),
#                      y.start + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], y.start) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], y.start)^2*dt[1]) )
#       if(missing(grid)) d <- diff(cand.Area)/cand.length
#       
#       result <- matrix(0, sample.length, n-1)
#       result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
#                               pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
#       
#       for(i in 2:(n-1)){
#         
#         if(!missing(b.fun.mat)){
#           VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
#           dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
#         }else{
#           if(pred.alg == "Distribution"){
#             VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
#             dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
#           }
#           if(pred.alg == "Trajectory"){
#             VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
#             dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
#           }
#         }
#         
#         cand.Area <- c(min(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]),
#                        max(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]) )
#         if(missing(grid)) d <- diff(cand.Area)/cand.length
#         
#         result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
#                                 pred.alg = pred.alg, sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
#         if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
#       }
#       
#     }
#   }
  
  if(plot.prediction){
    qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
    if(which.series == "new"){
      plot(object@t, object@Y, type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
      lines(t[-1], qu[1,], col = 2)
      lines(t[-1], qu[2,], col = 2)
    }
    if(which.series == "current"){
      plot(object@t, object@Y, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
      lines(t[-1], qu[1,], col = 2)
      lines(t[-1], qu[2,], col = 2)
    }
  }
  return(result)

})


########
#' Prediction for mixed diffusion process
#'
#' @description Bayesian prediction of a stochastic process
#'   \eqn{dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param object class object of MCMC samples: "est.mixedDiffusion"
#' @param t vector of time points to make predictions for
#' @param Euler.interval if TRUE: simple prediction intervals with Euler are made (in one step each)
#' @param level level of the prediction intervals
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param b.fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param y.start optional, if missing, first (which.series = "new") or last observation variable ("current") is taken
#' @param ind.pred index of series to be predicted, optional, if which.series = "current" and ind.pred missing, the last series is taken
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' mu <- 2; Omega <- 0.4; phi <- matrix(rnorm(21, mu, sqrt(Omega)))
#' model <- set.to.class("mixedDiffusion", 
#'          parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'          b.fun = function(phi, t, x) phi*x, sT.fun = function(t, x) x)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est_mixdiff <- estimate(model, t, data[1:20,], 2000)
#' plot(est_mixdiff)
#' \dontrun{
#' pred_mixdiff <- predict(est_mixdiff, b.fun.mat = function(phi, t, y) phi[,1]*y)
#' lines(t, data[21,], lwd = 2)
#' mean(apply(pred_mixdiff$Y, 2, quantile, 0.025) <= data[21, ] & 
#' apply(pred_mixdiff$Y, 2, quantile, 0.975) >= data[21, ])
#' mean(sapply(1:20, function(i){ 
#'    mean(apply(pred_mixdiff$Y, 2, quantile, 0.025) <= data[i, ] & 
#'    apply(pred_mixdiff$Y, 2, quantile, 0.975) >= data[i, ])}))
#' pred_mixdiff2 <- predict(est_mixdiff, b.fun.mat = function(phi, t, y) phi[,1]*y, 
#'      which.series = "current")
#' pred_mixdiff3 <- predict(est_mixdiff, b.fun.mat = function(phi, t, y) phi[,1]*y, 
#'      which.series = "current", y.start = data[20, 51], t = t[51:101])
#' }
#' pred_mixdiff <- predict(est_mixdiff, Euler.interval = TRUE, 
#'      b.fun.mat = function(phi, t, y) phi[,1]*y); lines(t, data[21,], lwd = 2)  
#'      # one step Euler approximation
#' pred_mixdiff <- predict(est_mixdiff, pred.alg = "simpleTrajectory", sample.length = 100)
#' for(i in 1:100) lines(t, pred_mixdiff$Y[i,], col = "grey")
#' pred_mixdiff <- predict(est_mixdiff, pred.alg = "simpleBayesTrajectory")
#'
#' # OU
#' b.fun <- function(phi, t, y) phi[1]-phi[2]*y; y0.fun <- function(phi, t) phi[3]
#' mu <- c(10, 1, 0.5); Omega <- c(0.9, 0.01, 0.01)
#' phi <- sapply(1:3, function(i) rnorm(21, mu[i], sqrt(Omega[i])))
#' cl <- set.to.class("mixedDiffusion", 
#'            parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'            y0.fun = y0.fun, b.fun = b.fun, sT.fun = function(t, x) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data[1:20,], 2000)
#' plot(est)
#' pred_mixdiff <- predict(est, b.fun.mat = function(phi, t, y) phi[,1]-phi[,2]*y)
#' lines(t, data[21,], lwd = 2)
#' mean(apply(pred_mixdiff$Y, 2, quantile, 0.025) <= data[21, ] & 
#'      apply(pred_mixdiff$Y, 2, quantile, 0.975) >= data[21, ])
#' mean(sapply(1:20, function(i){
#'     mean(apply(pred_mixdiff$Y, 2, quantile, 0.025) <= data[i, ] & 
#'     apply(pred_mixdiff$Y, 2, quantile, 0.975) >= data[i, ])}))
#'
#' @export
setMethod(f = "predict", signature = "est.mixedDiffusion",
          definition = function(object, t, Euler.interval = FALSE, level = 0.05, burnIn, thinning,
                                b.fun.mat, which.series = c("new", "current"), y.start, ind.pred, M2pred = 10,
                                cand.length = 1000, pred.alg = c("Distribution", "Trajectory", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, grid, plot.prediction = TRUE) {
# Y as list ?
    pred.alg <- match.arg(pred.alg)
    which.series <- match.arg(which.series)

#    if(is.list(object@Y)) J <- length(object@Y)
#    if(is.matrix(object@Y)) J <- nrow(object@Y)

    J <- nrow(object@Y)
    n <- ncol(object@Y)  # equal to length(object@t)

    if(missing(ind.pred) & which.series == "current") ind.pred <- J

    if(missing(burnIn)) burnIn <- object@burnIn
    if(missing(thinning)) thinning <- object@thinning

    if(missing(t)){
      if(which.series == "new") t <- object@t
      if(which.series == "current"){
        dt <- median( diff(object@t))
        t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
      }
    }
    ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
    K <- length(ind)
    if(which.series == "new"){
      samples <- list(mu = as.matrix(object@mu[ind,]), Omega = as.matrix(object@Omega[ind,]) )
      phi.pred <- predPhi(samples)

      samples <- list(phi = phi.pred, gamma2 = object@gamma2[ind])
      if(missing(y.start)){
        y.start <- sapply(1:K, function(i) object@model$y0.fun(samples$phi[i, ], t[1]))
      }
      
    }
    if(which.series == "current"){
      samples <- list(phi = sapply(1:ncol(object@mu) , function(i) phi_ij(object@phi[ind], ind.pred, i) ), gamma2 = object@gamma2[ind])
      if(missing(y.start)){
         y.start <- object@Y[ind.pred, n]
      }
    }

    if(missing(sample.length) | pred.alg == "Distribution") sample.length <- K
    b.fun <- object@model$b.fun
    sT.fun <- object@model$sT.fun
    n <- length(t)
    dt <- diff(t)

    if(pred.alg == "simpleTrajectory"){
      if(missing(sample.length)) sample.length <- 100
      result <- matrix(0, sample.length, n)
      if(which.series == "current"){
        phi.est <- apply(samples$phi, 2, mean); gamma2.est <- mean(samples$gamma2)
        cl <- set.to.class("Diffusion", parameter = list(phi = phi.est, gamma2 = gamma2.est), b.fun = b.fun, sT.fun = sT.fun)
        for(i in 1:sample.length){
          result[i,] <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
        }
      }
      if(which.series == "new"){
        mu.est <- apply(as.matrix(object@mu[ind,]), 2, mean)
        Omega.est <- apply(as.matrix(object@Omega[ind,]), 2, mean)
        gamma2.est <- mean(samples$gamma2)
        phi.pred <- matrix(0, sample.length, length(mu.est))
        for(i in 1:sample.length){
          phi.pred[i, ] <- rnorm(length(mu.est), mu.est, sqrt(Omega.est))
          cl <- set.to.class("Diffusion", parameter = list(phi = phi.pred[i, ], gamma2 = gamma2.est), b.fun = b.fun, sT.fun = sT.fun)
          result[i,] <- simulate(cl, t = t, y0 = object@model$y0.fun(phi.pred[i,], t[1]), plot.series = FALSE)
        }
      }
    }
    if(pred.alg == "simpleBayesTrajectory"){
      if(missing(sample.length)){
        sample.length <- K
      } else{
        if(sample.length > K) sample.length <- K
      }
      result <- matrix(0, sample.length, n)
      if(which.series == "current"){
        for(i in 1:sample.length){
          cl <- set.to.class("Diffusion", parameter = list(phi = samples$phi[i,], gamma2 = samples$gamma2[i]), b.fun = b.fun, sT.fun = sT.fun)
          result[i,] <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
        }
      }
      if(which.series == "new"){
        mu.est <- as.matrix(object@mu[ind,])
        Omega.est <- as.matrix(object@Omega[ind,])
        phi.pred <- matrix(0, sample.length, length(mu.est))
        for(i in 1:sample.length){
          phi.pred[i, ] <- rnorm(length(mu.est), mu.est[i,], sqrt(Omega.est[i,]))
          cl <- set.to.class("Diffusion", parameter = list(phi = phi.pred[i, ], gamma2 = samples$gamma2[i]), b.fun = b.fun, sT.fun = sT.fun)
          result[i,] <- simulate(cl, t = t, y0 = object@model$y0.fun(phi.pred[i,], t[1]), plot.series = FALSE)
        }
      }
    }
    if(pred.alg %in% c("Distribution", "Trajectory")){

      if(Euler.interval){
        result <- matrix(0, 2, n-1)
        for(i in 1:(n-1)){
          cand.Area <- c(min(y.start) + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*(t[i+1]-t[1])),
                         max(y.start) + (t[i+1]-t[1])*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*(t[i+1]-t[1])) )
          if(!missing(b.fun.mat)){
            VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1]))))
          }else{
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+(t[i+1]-t[1])*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*(t[i+1]-t[1])))))
          }
          result[,i] <- prediction.intervals(samples, VFun, x0 = y.start, level = level, candArea = cand.Area)
        }
        result <- cbind(quantile(y.start, c(level/2, 1-level/2)), result)
        
      }else{


        if(pred.alg == "Distribution") method <- "vector"
        if(pred.alg == "Trajectory") method <- "free"
        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
          dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
        }else{
          if(pred.alg == "Distribution"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
          }
          if(pred.alg == "Trajectory"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
          }
        }
        cand.Area <- c(min(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]),
                       max(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]) )
        if(missing(grid)) d <- diff(cand.Area)/cand.length

        result <- matrix(0, sample.length, n-1)
        result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)

        for(i in 2:(n-1)){

          if(!missing(b.fun.mat)){
            VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
            dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
          }else{
            if(pred.alg == "Distribution"){
              VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
              dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
            }
            if(pred.alg == "Trajectory"){
              VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
              dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
            }
          }

           cand.Area <- range(result[,i-1]) + c(-1, 1)*(which.series == "new") + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*c(-1,1)*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i])
          if(missing(grid)) d <- diff(cand.Area)/cand.length

          result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
                                  pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
          if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
        }
        result <- cbind(y.start, result)
        
      }

    }



    if(plot.prediction){
      qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
      if(which.series == "new"){
        plot(object@t, object@Y[1,], type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
        for(j in 2:J) lines(object@t, object@Y[j,])
        lines(t, qu[1,], col = 2, lwd = 2)
        lines(t, qu[2,], col = 2, lwd = 2)
      }
      if(which.series == "current"){
        plot(object@t, object@Y[ind.pred, ], type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y[ind.pred, ], range(qu))), xlab = "t", ylab = expression(Y[t]))
        lines(t, qu[1,], col = 2)
        lines(t, qu[2,], col = 2)
      }
    }
    if(which.series == "new") return(list(Y = result, phi = phi.pred))
    if(which.series == "current") return(result)

})



########
#' Prediction for noisy / hidden diffusion process
#'
#' @description Bayesian prediction of the model,
#'   \eqn{Z_i = Y_{t_i} + \epsilon_i, dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t}.
#' @param object class object of MCMC samples: "est.hiddenDiffusion"
#' @param t vector of time points to make predictions for
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param b.fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' model <- set.to.class("hiddenDiffusion", parameter = list(phi = 5, gamma2 = 1, sigma2 = 0.1))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est_hiddiff <- estimate(model, t, data$Z, 1000)
#' plot(est_hiddiff)
#' \dontrun{
#' pred_hiddiff <- predict(est_hiddiff, t = seq(0, 1, by = 0.1))
#' pred_hiddiff2 <- predict(est_hiddiff, which.series = "current")
#' }
#' pred_hiddiff <- predict(est_hiddiff, pred.alg = "simpleTrajectory", sample.length = 100)
#' pred_hiddiff <- predict(est_hiddiff, pred.alg = "simpleBayesTrajectory")
#' @export
setMethod(f = "predict", signature = "est.hiddenDiffusion",
          definition = function(object, t, burnIn, thinning,
                                b.fun.mat, which.series = c("new", "current"), M2pred = 10,
                                cand.length = 1000, pred.alg = c("Distribution", "Trajectory", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, grid, plot.prediction = TRUE) {

    pred.alg <- match.arg(pred.alg)
    which.series <- match.arg(which.series)

    if(missing(burnIn)) burnIn <- object@burnIn
    if(missing(thinning)) thinning <- object@thinning

    if(missing(t)){
      if(which.series == "new") t <- object@t
      if(which.series == "current"){
        dt <- median( diff(object@t))
        t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
      }
    }

    ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
    samples <- list(phi = as.matrix(object@phi[ind,]), gamma2 = object@gamma2[ind], sigma2 = object@sigma2[ind])
    K <- length(samples$gamma2)
    if(missing(sample.length) | pred.alg == "Distribution") sample.length <- K
    b.fun <- object@model$b.fun
    sT.fun <- object@model$sT.fun
    n <- length(t)
    dt <- diff(t)

    if(which.series == "new") y.start <- object@Y.est[ind, 1]
    if(which.series == "current") y.start <- object@Y.est[ind, ncol(object@Y.est)]

    if(pred.alg == "simpleTrajectory"){
      if(missing(sample.length)) sample.length <- 100
      result <- matrix(0, sample.length, n)
      phi.est <- apply(samples$phi, 2, mean); gamma2.est <- mean(samples$gamma2); sigma2.est <- mean(samples$sigma2)
      if(which.series == "new"){
        y0.fun <- object@model$y0.fun
      }else{
        y0.fun <- function(phi, t) y.start[sample(K, 1)]
      }
      
      cl <- set.to.class("hiddenDiffusion", parameter = list(phi = phi.est, gamma2 = gamma2.est, sigma2 = sigma2.est), y0.fun = y0.fun, b.fun = b.fun, sT.fun = sT.fun)
      Ypred <- matrix(0, sample.length, n)
     for(i in 1:sample.length){
        he <- simulate(cl, t = t, plot.series = FALSE)
        result[i,] <- he$Z
        Ypred[i,] <- he$Y
      }

    }
    if(pred.alg == "simpleBayesTrajectory"){
      if(missing(sample.length)){
        sample.length <- K
      } else{
        if(sample.length > K) sample.length <- K
      }
      result <- matrix(0, sample.length, n)
      if(which.series == "new"){
        y0.fun <- object@model$y0.fun
      }else{
        y0.fun <- function(phi, t) y.start[sample(K, 1)]
      }
      
      Ypred <- matrix(0, sample.length, n)
      for(i in 1:sample.length){
        cl <- set.to.class("hiddenDiffusion", parameter = list(phi = samples$phi[i,], gamma2 = samples$gamma2[i], sigma2 = samples$sigma2[i]), y0.fun = y0.fun, b.fun = b.fun, sT.fun = sT.fun)
        he <- simulate(cl, t = t, plot.series = FALSE)
        result[i,] <- he$Z
        Ypred[i,] <- he$Y
      }
      
      
    }  
    if(pred.alg %in% c("Distribution", "Trajectory")){

      if(pred.alg == "Distribution") method <- "vector"
      if(pred.alg == "Trajectory") method <- "free"
      if(!missing(b.fun.mat)){
        VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
        dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
      }else{
        if(pred.alg == "Distribution"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
        }
        if(pred.alg == "Trajectory"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
        }
      }
      cand.Area <- c(min(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]),
                     max(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]) )
      if(missing(grid)) d <- diff(cand.Area)/cand.length

      result <- matrix(0, sample.length, n-1)
      result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
                              pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)

      for(i in 2:(n-1)){

        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
          dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
        }else{
          if(pred.alg == "Distribution"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
          }
          if(pred.alg == "Trajectory"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
          }
        }

        cand.Area <- c(min(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]),
                       max(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]) )
        if(missing(grid)) d <- diff(cand.Area)/cand.length

        result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
        if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
      }

    Ypred <- cbind(y.start, result)
    N <- ncol(Ypred)
    result <- matrix(0, K, N)

    for(i in 1:N){
      cand.Area <- range(Ypred[, i]) + c(-1, 1)*sqrt(mean(samples$sigma2))*5
      VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1, sqrt(samples$sigma2)))
      result[,i] <- pred.base(samples = samples, VFun, dens, x0 = Ypred[,i], len = sample.length, method = method,
                              pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
    }


    }


    if(plot.prediction){
      qu <- apply(result, 2, quantile, c(0.05 / 2, 1 - 0.05 / 2))
      if(which.series == "new"){
        plot(object@t, object@Z, type = "l", xlab = "t", ylim = range(c(object@Z, range(qu))), ylab = "prediction")
        lines(t, qu[1,], col = 2)
        lines(t, qu[2,], col = 2)
      }
      if(which.series == "current"){
        plot(object@t, object@Z, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Z, range(qu))), xlab = "t", ylab = "prediction")
        lines(t, qu[1,], col = 2)
        lines(t, qu[2,], col = 2)
      }
    }
    return(list(Y = Ypred, Z = result))

  })





########
#' Prediction for noisy/hidden mixed diffusion process
#'
#' @description Bayesian prediction of a stochastic process
#'   \eqn{Z_{ij} = Y_{t_{ij}} + \epsilon_{ij}, dY_t = b(\phi_j,t,Y_t)dt + s(\gamma,t,Y_t)dW_t, \phi_j~N(\mu, \Omega)}.
#' @param object class object of MCMC samples: "est.hiddenmixedDiffusion"
#' @param t vector of time points to make predictions for
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param b.fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param ind.pred index of series to be predicted, optional, if which.series = "current" and ind.pred missing, the last series is taken
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' mu <- c(5, 1); Omega <- c(0.9, 0.04)
#' phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' y0.fun <- function(phi, t) phi[2]
#' model <- set.to.class("hiddenmixedDiffusion", y0.fun = y0.fun, 
#'      b.fun = function(phi, t, y) phi[1], 
#'      parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 1, sigma2 = 0.01))
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' \dontrun{
#' est_hidmixdiff <- estimate(model, t, data$Z[1:20,], 2000)
#' plot(est_hidmixdiff)
#' pred1 <- predict(est_hidmixdiff, b.fun.mat = function(phi, t, y) phi[,1])
#' pred2 <- predict(est_hidmixdiff, pred.alg = "Trajectory", b.fun.mat = function(phi, t, y) phi[,1])
#' pred3 <- predict(est_hidmixdiff, pred.alg = "simpleTrajectory", sample.length = nrow(pred1$Y))
#' lines(t, apply(pred1$Z, 2, quantile, 0.025), col = 3)
#' lines(t, apply(pred1$Z, 2, quantile, 0.975), col = 3)
#' lines(t, apply(pred2$Z, 2, quantile, 0.025), col = 4)
#' lines(t, apply(pred2$Z, 2, quantile, 0.975), col = 4)
#' pred4 <- predict(est_hidmixdiff, pred.alg = "simpleBayesTrajectory")
#' }
#' @export
setMethod(f = "predict", signature = "est.hiddenmixedDiffusion",
          definition = function(object, t, burnIn, thinning,
                                b.fun.mat, which.series = c("new", "current"), ind.pred, M2pred = 10,
                                cand.length = 1000, pred.alg = c("Distribution", "Trajectory", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, grid, plot.prediction = TRUE) {

    pred.alg <- match.arg(pred.alg)
    which.series <- match.arg(which.series)

    J <- nrow(object@Z)
    n <- ncol(object@Z)  # equal to length(object@t)

    if(missing(ind.pred) & which.series == "current") ind.pred <- J

    if(missing(burnIn)) burnIn <- object@burnIn
    if(missing(thinning)) thinning <- object@thinning

    if(missing(t)){
      if(which.series == "new") t <- object@t
      if(which.series == "current"){
        dt <- median( diff(object@t))
        t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
      }
    }
    ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
    K <- length(ind)


    if(which.series == "new"){
      samples <- list(mu = as.matrix(object@mu[ind,]), Omega = as.matrix(object@Omega[ind,]) )
      phi.pred <- predPhi(samples)

      samples <- list(phi = phi.pred, gamma2 = object@gamma2[ind], sigma2 = object@sigma2[ind])
    }
    if(which.series == "current"){
      samples <- list(phi = sapply(1:ncol(object@mu) , function(i) phi_ij(object@phi[ind], ind.pred, i) ), gamma2 = object@gamma2[ind], sigma2 = object@sigma2[ind])
    }
    if(which.series == "new") y.start <- sapply(1:K, function(i) object@model$y0.fun(samples$phi[i, ], t[1]))
    if(which.series == "current") y.start <- sapply(object@Y.est[ind], function(li) li[[ind.pred]][length(li[[ind.pred]])])

    b.fun <- object@model$b.fun
    sT.fun <- object@model$sT.fun
    n <- length(t)
    dt <- diff(t)

    if(pred.alg == "simpleTrajectory"){
      if(missing(sample.length)) sample.length <- 100

      result <- matrix(0, sample.length, n)
      if(which.series == "current"){
        phi.est <- apply(samples$phi, 2, mean); gamma2.est <- mean(samples$gamma2)

        cl <- set.to.class("Diffusion", parameter = list(phi = phi.est, gamma2 = gamma2.est), y0.fun = function(phi, t) y.start[sample(K, 1)], b.fun = b.fun, sT.fun = sT.fun)
        for(i in 1:sample.length){
          result[i,] <- simulate(cl, t = t, plot.series = FALSE)
        }
      }
      if(which.series == "new"){
        mu.est <- apply(as.matrix(object@mu[ind,]), 2, mean)
        Omega.est <- apply(as.matrix(object@Omega[ind,]), 2, mean)
        gamma2.est <- mean(samples$gamma2)
        phi.pred <- matrix(0, sample.length, length(mu.est))

        for(i in 1:sample.length){
          phi.pred[i, ] <- rnorm(length(mu.est), mu.est, sqrt(Omega.est))
          cl <- set.to.class("Diffusion", parameter = list(phi = phi.pred[i, ], gamma2 = gamma2.est), b.fun = b.fun, sT.fun = sT.fun)
          result[i,] <- simulate(cl, t = t, y0 = object@model$y0.fun(phi.pred[i,], t[1]), plot.series = FALSE)
        }
      }
      Ypred <- result
      result <- apply(Ypred, 2, function(vec) rnorm(length(vec), vec, sqrt(mean(samples$sigma2))))
    }
    if(pred.alg == "simpleBayesTrajectory"){
      if(missing(sample.length)){
        sample.length <- K
      } else{
        if(sample.length > K) sample.length <- K
      }
      
      result <- Ypred <- matrix(0, sample.length, n)
      if(which.series == "current"){

        for(i in 1:sample.length){
          cl <- set.to.class("hiddenDiffusion", parameter = list(phi = samples$phi[i,], gamma2 = samples$gamma2[i], sigma2 = samples$sigma2[i]), y0.fun = function(phi, t) y.start[sample(K, 1)], b.fun = b.fun, sT.fun = sT.fun)
          he <- simulate(cl, t = t, plot.series = FALSE)
          Ypred[i,] <- he$Y
          result[i,] <- he$Z
        }
      }
      if(which.series == "new"){
        mu.est <- as.matrix(object@mu[ind,])
        Omega.est <- as.matrix(object@Omega[ind,])
        phi.pred <- matrix(0, sample.length, length(mu.est))
        
        for(i in 1:sample.length){
          phi.pred[i, ] <- rnorm(length(mu.est), mu.est[i,], sqrt(Omega.est[i,]))
          cl <- set.to.class("hiddenDiffusion", parameter = list(phi = phi.pred[i, ], gamma2 = samples$gamma2[i], sigma2 = samples$sigma2[i]), 
                             b.fun = b.fun, sT.fun = sT.fun, y0.fun = object@model$y0.fun)
          he <- simulate(cl, t = t, plot.series = FALSE)
          result[i,] <- he$Z
          Ypred[i,] <- he$Y
        }
      }
    }
    
    if(pred.alg %in% c("Distribution", "Trajectory")){
      if(missing(sample.length)) sample.length <- K

      if(pred.alg == "Distribution") method <- "vector"
      if(pred.alg == "Trajectory") method <- "free"
      if(!missing(b.fun.mat)){
        VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
        dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[1]*b.fun.mat(samples$phi, t[1], yn_1), sqrt(samples$gamma2*sT.fun(t[1], yn_1)^2*dt[1])))
      }else{
        if(pred.alg == "Distribution"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[1]*b.fun(samples$phi[k,], t[1], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1[k])^2*dt[1]))))
        }
        if(pred.alg == "Trajectory"){
          VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
          dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[1]*b.fun(samples$phi[k,], t[1], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[1], yn_1)^2*dt[1]))))
        }
      }
      cand.Area <- c(min(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]),
                     max(y.start) + dt[1]*b.fun(apply(samples$phi, 2, mean), t[1], mean(y.start)) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[1], mean(y.start))^2*dt[1]) )
      if(missing(grid)) d <- diff(cand.Area)/cand.length

      result <- matrix(0, sample.length, n-1)
      result[,1] <- pred.base(samples = samples, VFun, dens, x0 = y.start, len = sample.length, method = method,
                              pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)

      for(i in 2:(n-1)){

        if(!missing(b.fun.mat)){
          VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
          dens <- function(yn, yn_1, samples)  mean(dnorm(yn, yn_1+dt[i]*b.fun.mat(samples$phi, t[i], yn_1), sqrt(samples$gamma2*sT.fun(t[i], yn_1)^2*dt[i])))
        }else{
          if(pred.alg == "Distribution"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1[k]+dt[i]*b.fun(samples$phi[k,], t[i], yn_1[k]), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1[k])^2*dt[i]))))
          }
          if(pred.alg == "Trajectory"){
            VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
            dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, yn_1+dt[i]*b.fun(samples$phi[k,], t[i], yn_1), sqrt(samples$gamma2[k]*sT.fun(t[i], yn_1)^2*dt[i]))))
          }
        }

        cand.Area <- c(min(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) - 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]),
                       max(result[,i-1]) + dt[i]*b.fun(apply(samples$phi, 2, mean), t[i], mean(result[,i-1])) + 5*sqrt(mean(samples$gamma2)*sT.fun(t[i], mean(result[,i-1]))^2*dt[i]) )
        if(missing(grid)) d <- diff(cand.Area)/cand.length

        result[,i] <- pred.base(samples = samples, VFun, dens, x0 = result[,i-1], len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
        if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
      }
      Ypred <- cbind(y.start, result)
      N <- ncol(Ypred)
      result <- matrix(0, K, N)  # different sample.length ?

      for(i in 1:N){
        cand.Area <- range(Ypred[, i]) + c(-1, 1)*sqrt(mean(samples$sigma2))*5
        VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1, sqrt(samples$sigma2)))
        result[,i] <- pred.base(samples = samples, VFun, dens, x0 = Ypred[,i], len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
      }

    }

    if(plot.prediction){
      qu <- apply(result, 2, quantile, c(0.05 / 2, 1 - 0.05 / 2))
      if(which.series == "new"){
        plot(object@t, object@Z[1,], type = "l", xlab = "t", ylim = range(c(object@Z, range(qu))), ylab = expression(Z[t]))
        for(j in 2:J) lines(object@t, object@Z[j,])
        lines(t, qu[1,], col = 2, lwd = 2)
        lines(t, qu[2,], col = 2, lwd = 2)
      }
      if(which.series == "current"){
        plot(object@t, object@Z[ind.pred, ], type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Z[ind.pred, ], range(qu))), xlab = "t", ylab = expression(Z[t]))
        lines(t, qu[1,], col = 2)
        lines(t, qu[2,], col = 2)
      }
    }

    if(which.series == "new") return(list(Z = result, Y = Ypred, phi = phi.pred))
    if(which.series == "current") return(list(Z = result, Y = Ypred))

})






########
#' Prediction for Poisson process
#'
#' @description Bayesian prediction of a nonhomogeneous Poisson process.
#' @param object class object of MCMC samples: "est.NHPP"
#' @param variable if prediction of event times ("eventTimes") or of Poisson process variables ("PoissonProcess")
#' @param t vector of time points to make predictions for (only for variable = "PoissonProcess")
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param Lambda.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param Tstart optional, if missing, first (which.series = "new") or last observation variable ("current") is taken
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param grid fineness degress of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)), 
#'                Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, plot.series = TRUE)
#' est_NHPP <- estimate(model, t, data$Times, 2000)
#' plot(est_NHPP)
#' \dontrun{
#' pred_NHPP <- predict(est_NHPP)
#' pred_NHPP <- predict(est_NHPP, variable = "PoissonProcess")
#' pred_NHPP2 <- predict(est_NHPP, which.series = "current")
#' pred_NHPP3 <- predict(est_NHPP, variable = "PoissonProcess", which.series = "current")
#' }
#' pred_NHPP <- predict(est_NHPP, pred.alg = "simpleTrajectory", M2pred = length(data$Times))
#' pred_NHPP <- predict(est_NHPP, variable = "PoissonProcess", pred.alg = "simpleTrajectory", 
#'                      M2pred = length(data$Times))
#' pred_NHPP <- predict(est_NHPP, variable = "PoissonProcess", pred.alg = "simpleBayesTrajectory", 
#'                      M2pred = length(data$Times), sample.length = 100)
#'
#' @export
setMethod(f = "predict", signature = "est.NHPP",
          definition = function(object, variable = c("eventTimes", "PoissonProcess"),
                                t, burnIn, thinning,
                                Lambda.mat, which.series = c("new", "current"), Tstart, M2pred = 10,
                                cand.length = 1000, pred.alg = c("Trajectory", "Distribution", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, grid = 1e-05, plot.prediction = TRUE) {

    pred.alg <- match.arg(pred.alg)
    which.series <- match.arg(which.series)
    variable <- match.arg(variable)

    if(missing(Tstart)){
      if(which.series == "new") Tstart <- 0
      if(which.series == "current"){
        Tstart <- max(object@t)
      }
    }

    if(missing(burnIn)) burnIn <- object@burnIn
    if(missing(thinning)) thinning <- object@thinning

    if(missing(t) & variable == "PoissonProcess"){
      if(which.series == "new"){
        t <- object@t
      }
      if(which.series == "current"){
        dt <- median( diff(object@t))
        t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
      }
    }
    K <- nrow(object@xi)
    ind <- seq(burnIn + 1, K, by = thinning)
    K <- length(ind)
    samples <- object@xi[ind, ]

    Lambda <- object@model$Lambda

    if(missing(Lambda.mat)){
      if(pred.alg == "Distribution"){
        Lambda.diff <- function(cand, Tn_1, samples){
          sapply(1:K, function(i) Lambda(cand, samples[i,]) - Lambda(Tn_1[i], samples[i,]) )
        }

        Fun <- function(cand, Tn_1){
          sol <- 1 - exp(- Lambda.diff(cand, Tn_1, samples))
          sum(sol[sol >= 0])/K
        }
      }else{
        Lambda.diff <- function(cand, Tn_1, samples){
          sapply(1:K, function(i) Lambda(cand + Tn_1, samples[i,]) - Lambda(Tn_1, samples[i,]) )
        }

        Fun <- function(cand, Tn_1) 1 - mean(exp(- Lambda.diff(cand, Tn_1, samples)))
      }
    }else{

      if(pred.alg == "Distribution"){
        Fun <- function(cand, Tn_1){
          sol <- 1 - exp(-(Lambda.mat(cand, samples) - Lambda.mat(Tn_1, samples)))
          sum(sol[sol >= 0])/K
        }
      }else{
        Fun <- function(cand, Tn_1){
          1-mean(exp(-(Lambda.mat(cand+Tn_1, samples)-Lambda.mat(Tn_1, samples))))
        }
      }

    }


    if( variable == "eventTimes"){

      if(pred.alg == "simpleTrajectory"){
        if(missing(sample.length)) sample.length <- 100
        xi <- apply(samples, 2, mean)
        drawTn <- function(Tn_1){
          cand <- seq(Tn_1, Tn_1 + 1, by = 1e-04)
          prob <- 1 - exp(-(Lambda(cand, xi)-Lambda(Tn_1, xi)))
          u <- runif(1)
          cand[which(abs(u-prob) == min(abs(u-prob)))]
        }
        result <- matrix(0, sample.length, M2pred + 1)
        result[, 1] <- Tstart
        for(i in 1:sample.length){
          for(j in 2:(M2pred + 1)){
            result[i,j] <- drawTn(result[i, j-1])
          }
        }
        result <- result[, -1]
      }
      
      if(pred.alg == "simpleBayesTrajectory"){
        if(missing(sample.length)){
          sample.length <- K
        } else{
          if(sample.length > K) sample.length <- K
        }
        drawTn <- function(Tn_1, xi){
          cand <- seq(Tn_1, Tn_1 + 1, by = 1e-04)
          prob <- 1 - exp(-(Lambda(cand, xi)-Lambda(Tn_1, xi)))
          u <- runif(1)
          cand[which(abs(u-prob) == min(abs(u-prob)))]
        }
        result <- matrix(0, sample.length, M2pred + 1)
        result[, 1] <- Tstart
        for(i in 1:sample.length){
          for(j in 2:(M2pred + 1)){
            result[i,j] <- drawTn(result[i, j-1], samples[i,])
          }
        }
        result <- result[, -1]
      }
      
      if(pred.alg == "Trajectory"){

        if(missing(sample.length)) sample.length <- 100

        drawTn <- function(Tn_1){
          u <- runif(1)
          cand <- 1
          memory <- cand

          while(length(unique(memory)) == length(memory)){
            if(Fun(cand, Tn_1) < u){
              cand <- cand*2
              memory <- c(memory, cand)
            }else{
              cand <- cand/2
              memory <- c(memory, cand)
            }
          }
          lower <- min(memory[length(memory)-1], memory[length(memory)])
          upper <- max(memory[length(memory)-1], memory[length(memory)])
          diff <- upper - lower
          while(diff >= grid){
            if(Fun(lower+diff/2, Tn_1) < u){
              lower <- lower+diff/2
            }else{
              upper <- lower+diff/2
            }
            diff <- upper - lower
          }
          Tn_1 + (lower+upper)/2
        }
        result <- matrix(0, sample.length, M2pred)
        result[,1] <- sapply(1:sample.length, function(i) drawTn(Tstart))

        if(M2pred > 1){
          for(i in 2:M2pred){
            result[,i] <- vapply(result[,i-1], drawTn, FUN.VALUE = numeric(1))
          }
        }

      }
      if(pred.alg == "Distribution"){

        if(missing(sample.length)) sample.length <- K

        result <- matrix(0, sample.length, M2pred)
        Tn_1 <- rep(Tstart, sample.length)
        for(l in 1:M2pred){
          cand <- seq(0, max(Tn_1)/2 + 1, length = cand.length)
          diFu <- vapply(cand, Fun, Tn_1, FUN.VALUE = numeric(1))
          U <- runif(sample.length, 0, max(diFu))
          Tn_1 <- sapply(U, function(u) cand[which(diFu >= u)[1]])
          result[,l] <- Tn_1
      }

    }

      if(plot.prediction){
        qu <- apply(result, 2, quantile, c(0.05 / 2, 1 - 0.05 / 2))
        if(which.series == "new"){
          plot(object@jumpTimes, type = "l", xlab = "event number", ylim = range(c(object@jumpTimes, range(qu))), ylab = "event times")
          lines(qu[1,], col = 2)
          lines(qu[2,], col = 2)
        }
        if(which.series == "current"){
          plot(object@jumpTimes, type = "l", xlim = c(1, length(object@jumpTimes) + M2pred), ylim = range(c(object@jumpTimes, range(qu))), xlab = "event number", ylab = "event times")
          lines(length(object@jumpTimes) + 1:M2pred, qu[1,], col = 2)
          lines(length(object@jumpTimes) + 1:M2pred, qu[2,], col = 2)
        }
      }
    }else{

      n <- length(t)

      if(pred.alg == "Distribution"){
        message("does not work")
        return(NULL)
      }
      if(pred.alg == "simpleTrajectory"){
        if(missing(sample.length)) sample.length <- 100
        xi <- apply(samples, 2, mean)

        if(which.series == "new"){
            cl <- set.to.class("NHPP", parameter = list(xi = xi), Lambda = Lambda)
            result <- matrix(0, sample.length, n)  # here: N_t
            for(i in 1:sample.length){
              result[i, ] <- simulate(cl, t = t, plot.series = FALSE)$N
            }

        }else{
          result <- matrix(0, sample.length, n)
          drawTn <- function(Tn_1){
            cand <- seq(Tn_1, Tn_1 + 1, by = 1e-04)
            prob <- 1 - exp(-(Lambda(cand, xi)-Lambda(Tn_1, xi)))
            u <- runif(1)
            cand[which(abs(u-prob) == min(abs(u-prob)))]
          }

          result[, 1] <- max(object@N)
          for(i in 1:sample.length){
            times <- drawTn(Tstart)
            while(times[length(times)] < t[n]){
              times <- c(times, drawTn(times[length(times)]))
            }
            result[i, ] <- result[i, 1] + TimestoN(times, t)
          }

        }

      }
      
      if(pred.alg == "simpleBayesTrajectory"){
        if(missing(sample.length)){
          sample.length <- K
        } else{
          if(sample.length > K) sample.length <- K
        }

        if(which.series == "new"){
          result <- matrix(0, sample.length, n)  # here: N_t
          for(i in 1:sample.length){
            cl <- set.to.class("NHPP", parameter = list(xi = samples[i,]), Lambda = Lambda)
            result[i, ] <- simulate(cl, t = t, plot.series = FALSE)$N
          }
          
        }else{
          result <- matrix(0, sample.length, n)
          drawTn <- function(Tn_1, xi){
            cand <- seq(Tn_1, Tn_1 + 1, by = 1e-04)
            prob <- 1 - exp(-(Lambda(cand, xi)-Lambda(Tn_1, xi)))
            u <- runif(1)
            cand[which(abs(u-prob) == min(abs(u-prob)))]
          }
          
          result[, 1] <- max(object@N)
          for(i in 1:sample.length){
            times <- drawTn(Tstart)
            while(times[length(times)] < t[n]){
              times <- c(times, drawTn(times[length(times)], samples[i,]))
            }
            result[i, ] <- result[i, 1] + TimestoN(times, t)
          }
          
        }
        
      }
      if(pred.alg == "Trajectory"){

        if(missing(sample.length)) sample.length <- 100

        drawTn <- function(Tn_1){
          u <- runif(1)
          cand <- 1
          memory <- cand

          while(length(unique(memory)) == length(memory)){
            if(Fun(cand, Tn_1) < u){
              cand <- cand*2
              memory <- c(memory, cand)
            }else{
              cand <- cand/2
              memory <- c(memory, cand)
            }
          }
          lower <- min(memory[length(memory)-1], memory[length(memory)])
          upper <- max(memory[length(memory)-1], memory[length(memory)])
          diff <- upper - lower
          while(diff >= grid){
            if(Fun(lower+diff/2, Tn_1) < u){
              lower <- lower+diff/2
            }else{
              upper <- lower+diff/2
            }
            diff <- upper - lower
          }
          Tn_1 + (lower+upper)/2
        }

        result <- matrix(0, sample.length, n)  # here: N_t
        if(which.series == "new") result[, 1] <- 0
        if(which.series == "current") result[, 1] <- max(object@N)
        for(i in 1:sample.length){

          times <- drawTn(Tstart)
          while(times[length(times)] < t[n]){
            times <- c(times, drawTn(times[length(times)]))
          }
          result[i, ] <- result[i, 1] + TimestoN(times, t)
        }

      }

      if(plot.prediction){
        qu <- apply(result, 2, quantile, c(0.05 / 2, 1 - 0.05 / 2))
        if(which.series == "new"){
          plot(object@t, object@N, type = "l", lwd = 2, xlab = "t", ylim = range(c(object@N, range(qu))), ylab = expression(N[t]))
          lines(t, qu[1,], col = 2)
          lines(t, qu[2,], col = 2)
        }
        if(which.series == "current"){
          plot(object@t, object@N, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@N, range(qu))), xlab = "t", ylab = expression(N[t]))
          lines(t, qu[1,], col = 2)
          lines(t, qu[2,], col = 2)
        }
      }
    }





    return(result)

})


########
#' Prediction for jump diffusion process
#'
#' @description Bayesian prediction of a stochastic process
#'   \eqn{dY_t = b(\phi,t,Y_t)dt + s(\gamma,t,Y_t)dW_t + h(\eta,t,Y_t)dN_t}.
#' @param object class object of MCMC samples: "est.jumpDiffusion"
#' @param t vector of time points to make predictions for
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param Lambda.mat matrix-wise definition of intensity rate function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param M2pred default 10, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector"), for jump diffusion
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' model <- set.to.class("jumpDiffusion", 
#'          parameter = list(theta = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)), 
#'          Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
#' est_jd <- estimate(model, t, data, 2000)
#' plot(est_jd)
#' \dontrun{
#' pred_jd <- predict(est_jd, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
#' est_jd2 <- estimate(model, t[1:81], data = list(N = data$N[1:81], Y = data$Y[1:81]), 2000)
#' pred_jd2 <- predict(est_jd2, t = t[81:101], which.series = "current", 
#'                      Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
#' lines(t, data$Y, type="l", lwd = 2)
#' pred_jd3 <- predict(est_jd, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1], 
#'                      pred.alg = "Distribution")
#' }
#' pred_jd4 <- predict(est_jd, pred.alg = "simpleTrajectory", sample.length = 100)
#' for(i in 1:100) lines(t[-1], pred_jd4$Y[i,], col = "grey")
#' @export
setMethod(f = "predict", signature = "est.jumpDiffusion",
          definition = function(object, t, burnIn, thinning, Lambda.mat,
                                which.series = c("new", "current"), M2pred = 10,
                                cand.length = 1000, pred.alg = c("Trajectory", "Distribution", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, plot.prediction = TRUE) {

    pred.alg <- match.arg(pred.alg)
    which.series <- match.arg(which.series)

    if(missing(burnIn)) burnIn <- object@burnIn
    if(missing(thinning)) thinning <- object@thinning
    Lambda <- object@model$Lambda
    b.fun <- object@model$b.fun
    s.fun <- object@model$s.fun
    h.fun <- object@model$h.fun

    ind <- seq(burnIn + 1, ncol(object@xi), by = thinning)
    K <- length(ind)

    if(which.series == "new"){
      Tstart <- 0
      y.start <- object@Y[1]
      if(missing(t)) t <- object@t
    }

    if(which.series == "current"){
      Tstart <- max(object@t)
      y.start <- object@Y[length(object@Y)]
      if(missing(t)){
        dt <- median( diff(object@t))
        t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
      }
      if(length(object@N.est) == 0){  # -> object@N is equal to observed Poisson process
        startN <- max(object@N)
      }else{  # Poisson variables are filtered
        startN <- object@N.est[nrow(object@N.est), ind]
      }
    }

    n <- length(t)

    if(pred.alg == "Trajectory" | pred.alg == "Distribution"){

      gridN <- median(diff(t))/10
      
      samples <- t(object@xi[, ind])

      if(missing(Lambda.mat)){
        Lambda.diff <- function(cand, Tn_1, samples){
          sapply(1:K, function(i) Lambda(cand + Tn_1, samples[i,]) - Lambda(Tn_1, samples[i,]) )
        }
        Fun <- function(cand, Tn_1) 1 - mean(exp(- Lambda.diff(cand, Tn_1, samples)))
      }else{
        Fun <- function(cand, Tn_1){
          1-mean(exp(-(Lambda.mat(cand + Tn_1, samples)-Lambda.mat(Tn_1, samples))))
        }
      }
      sample.lengthN <- K

      drawTn <- function(Tn_1){
        u <- runif(1)
        cand <- Tn_1 + 0.1
        memory <- cand

        while(length(unique(memory)) == length(memory)){
          if(Fun(cand, Tn_1) < u){
            cand <- cand*2
            memory <- c(memory, cand)
          }else{
            cand <- cand/2
            memory <- c(memory, cand)
          }
        }
        lower <- min(memory[length(memory)-1], memory[length(memory)])
        upper <- max(memory[length(memory)-1], memory[length(memory)])
        diff <- upper - lower
        while(diff >= gridN){
          if(Fun(lower+diff/2, Tn_1) < u){
            lower <- lower+diff/2
          }else{
            upper <- lower+diff/2
          }
          diff <- upper - lower
        }
        Tn_1 + (lower+upper)/2
      }

      result <- matrix(0, sample.lengthN, n)  # here: N_t
      if(which.series == "new") result[, 1] <- 0
      if(which.series == "current") result[, 1] <- startN
      for(i in 1:sample.lengthN){

        times <- drawTn(Tstart)
        while(times[length(times)] < t[n]){
          times <- c(times, drawTn(times[length(times)]))
        }
        result[i, ] <- result[i, 1] + TimestoN(times, t)
        if(i %% 100 == 0) message(paste(i, "of", sample.lengthN, " Poisson process predictions are calculated"))
      }
      Npred <- result

      dN <- t(apply(Npred, 1, diff))
## jump diffusion:

      dt <- diff(t)
      samples <- list(phi = object@phi[ind], theta = object@theta[ind], gamma2 = object@gamma2[ind])
      if(missing(sample.length) | pred.alg == "Distribution") sample.length <- K

      if(pred.alg == "Distribution") method <- "vector"
      if(pred.alg == "Trajectory") method <- "free"

      VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[1]*b.fun(samples$phi, t[1], yn_1) + h.fun(samples$theta, t[1], yn_1)*dN[,1], sqrt(s.fun(samples$gamma2, t[1], yn_1)^2*dt[1])))
      cand.Area <- y.start + dt[1]*b.fun(mean(samples$phi), t[1], y.start) + h.fun(mean(samples$theta), t[1], y.start)*range(dN[,1]) + 5* c(-1,1) *sqrt(s.fun(mean(samples$gamma2), t[1], y.start)^2*dt[1])
      d <- diff(cand.Area)/cand.length

      result <- matrix(0, sample.length, n-1)
      result[,1] <- pred.base(samples = samples, VFun, x0 = y.start, len = sample.length, method = method,
                              pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)

      for(i in 2:(n-1)){

        VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, yn_1+dt[i]*b.fun(samples$phi, t[i], yn_1) + h.fun(samples$theta, t[i], yn_1)*dN[,i], sqrt(s.fun(samples$gamma2, t[i], yn_1)^2*dt[i])))

        cand.Area <- range(result[,i-1]) + dt[i]*b.fun(mean(samples$phi), t[i], mean(result[,i-1])) +
                      h.fun(mean(samples$theta), t[i], mean(result[,i-1]))*range(dN[,i]) + c(-1, 1)*5*sqrt(s.fun(mean(samples$gamma2), t[i], mean(result[,i-1]))^2*dt[i])
        d <- diff(cand.Area)/cand.length

        result[,i] <- pred.base(samples = samples, VFun, x0 = result[,i-1], len = sample.length, method = method,
                                pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
        if(i %% 10 == 0) message(paste(i, "of", n-1, "jump diffusion predictions are calculated"))
      }

    }



    if(pred.alg == "simpleTrajectory"){
      if(missing(sample.length)) sample.length <- 100
      xi <- apply(object@xi[,ind], 1, mean)
      phi <- mean(object@phi[ind])
      gamma2 <- mean(object@gamma2[ind])
      theta <- mean(object@theta[ind])

      if(which.series == "new"){
        cl <- set.to.class("jumpDiffusion", parameter = list(phi = phi, theta = theta, gamma2 = gamma2, xi = xi),
                           Lambda = Lambda, b.fun = b.fun, s.fun = s.fun, h.fun = h.fun)
        result <- matrix(0, sample.length, n-1)
        Npred <- matrix(0, sample.length, n)  # here: N_t

        for(i in 1:sample.length){
          help <- simulate(cl, t = t, y0 = object@Y[1], plot.series = FALSE)
          result[i, ] <- help$Y[-1]
          Npred[i, ] <- help$N
        }

      }else{
        result <- matrix(0, sample.length, n-1)
        Npred <- matrix(0, sample.length, n)  # here: N_t
        
        for(i in 1:sample.length){
          if(length(startN) > 1){
            ind.i <- sample(K, 1)
            Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN[ind.i]), Lambda = Lambda)$N
          }else{
            ind.i <- 1
            Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
          }
          
          result[i, ] <- sim_JD_Euler(t-max(object@t), phi, theta, gamma2, b.fun, s.fun, h.fun, start = y.start, Npred[i,] - startN[ind.i])[-1]
        }
      }

    }
    
    if(pred.alg == "simpleBayesTrajectory"){
      if(missing(sample.length)){
        sample.length <- K
      } else{
        if(sample.length > K) sample.length <- K
      }
      
      xi <- object@xi[,ind]
      phi <- object@phi[ind]
      gamma2 <- object@gamma2[ind]
      theta <- object@theta[ind]
      
      if(which.series == "new"){
        result <- matrix(0, sample.length, n-1)
        Npred <- matrix(0, sample.length, n)  # here: N_t
        
        for(i in 1:sample.length){
          cl <- set.to.class("jumpDiffusion", parameter = list(phi = phi[i], theta = theta[i], gamma2 = gamma2[i], xi = xi[,i]),
                             Lambda = Lambda, b.fun = b.fun, s.fun = s.fun, h.fun = h.fun)
          help <- simulate(cl, t = t, y0 = object@Y[1], plot.series = FALSE)
          result[i, ] <- help$Y[-1]
          Npred[i, ] <- help$N
        }
        
      }else{
        result <- matrix(0, sample.length, n-1)
        Npred <- matrix(0, sample.length, n)  # here: N_t
        
        for(i in 1:sample.length){
          if(length(startN) > 1){
            ind.i <- sample(K, 1)
            Npred[i, ] <- simN(t, xi[,i], len = 1, start = c(max(object@t), startN[ind.i]), Lambda = Lambda)$N
          }else{
            ind.i <- 1
            Npred[i, ] <- simN(t, xi[,i], len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
          }
          
          result[i, ] <- sim_JD_Euler(t-max(object@t), phi[i], theta[i], gamma2[i], b.fun, s.fun, h.fun, start = y.start, Npred[i,] - startN[ind.i])[-1]
        }
      }
      
    }

    if(plot.prediction){
      qu <- apply(result, 2, quantile, c(0.05 / 2, 1 - 0.05 / 2))
      if(which.series == "new"){
        plot(object@t, object@Y, type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
        lines(t[-1], qu[1,], col = 2)
        lines(t[-1], qu[2,], col = 2)
      }
      if(which.series == "current"){
        plot(object@t, object@Y, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
        lines(t[-1], qu[1,], col = 2)
        lines(t[-1], qu[2,], col = 2)
      }
    }
    return(list(N = Npred[, -1], Y = result))

})



########
#' Prediction for jump diffusion process
#'
#' @description Bayesian prediction of a stochastic process
#'   \eqn{Y_t = y_0 \exp( \phi t - \gamma2/2 t+\gamma W_t + \log(1+\theta) N_t)}.
#' @param object class object of MCMC samples: "est.Merton"
#' @param t vector of time points to make predictions for
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param Lambda.mat matrix-wise definition of intensity rate function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param M2pred default 10, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param only.interval if TRUE: only calculation of prediction intervals (only for pred.alg = "Distribution")
#' @param level level of the prediction intervals
#' @param cand.length length of candidate samples (if method = "vector"), for jump diffusion
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"
#' @param sample.length length of samples to be drawn
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' cl <- set.to.class("Merton", 
#'                parameter = list(thetaT = 0.1, phi = 0.05, gamma2 = 0.1, xi = c(3, 1/4)), 
#'                Lambda = function(t, xi) (t/xi[2])^xi[1])
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, y0 = 0.5, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' \dontrun{
#' pred1 <- predict(est, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
#' pred2 <- predict(est, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1], pred.alg = "Trajectory")
#' pred3 <- predict(est, pred.alg = "simpleTrajectory")
#' pred4 <- predict(est, pred.alg = "simpleBayesTrajectory")
#' }
#' @export
setMethod(f = "predict", signature = "est.Merton",
          definition = function(object, t, burnIn, thinning, Lambda.mat,
                                which.series = c("new", "current"), M2pred = 10, only.interval = TRUE, level = 0.05,
                                cand.length = 1000, pred.alg = c("Distribution", "Trajectory", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, plot.prediction = TRUE) {
            
            pred.alg <- match.arg(pred.alg)
            which.series <- match.arg(which.series)
            
            if(missing(burnIn)) burnIn <- object@burnIn
            if(missing(thinning)) thinning <- object@thinning
            Lambda <- object@model$Lambda
           
            ind <- seq(burnIn + 1, ncol(object@xi), by = thinning)
            K <- length(ind)
            
            if(which.series == "new"){
              Tstart <- 0
              y.start <- object@Y[1]
              if(missing(t)) t <- object@t
            }
            
            if(which.series == "current"){
              Tstart <- max(object@t)
              y.start <- object@Y[length(object@Y)]
              if(missing(t)){
                dt <- median( diff(object@t))
                t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
              }
              if(length(object@N.est) == 0){  # -> object@N is equal to observed Poisson process
                startN <- max(object@N)
              }else{  # Poisson variables are filtered
                startN <- object@N.est[nrow(object@N.est), ind]
              }
            }
            
            n <- length(t)
            
            if(pred.alg == "Trajectory" | pred.alg == "Distribution"){
              
              samples <- t(object@xi[, ind])
              
              if(missing(Lambda.mat)){
                Lambda.diff <- function(cand, Tn_1, samples){
                  sapply(1:K, function(i) Lambda(cand + Tn_1, samples[i,]) - Lambda(Tn_1, samples[i,]) )
                }
                Fun <- function(cand, Tn_1) 1 - mean(exp(- Lambda.diff(cand, Tn_1, samples)))
              }else{
                Fun <- function(cand, Tn_1){
                  1-mean(exp(-(Lambda.mat(cand + Tn_1, samples)-Lambda.mat(Tn_1, samples))))
                }
              }
              sample.lengthN <- K
              gridN <- median(diff(t))/10
              
              drawTn <- function(Tn_1){
                u <- runif(1)
                cand <- Tn_1 + 0.1
                memory <- cand
                
                while(length(unique(memory)) == length(memory)){
                  if(Fun(cand, Tn_1) < u){
                    cand <- cand*2
                    memory <- c(memory, cand)
                  }else{
                    cand <- cand/2
                    memory <- c(memory, cand)
                  }
                }
                lower <- min(memory[length(memory)-1], memory[length(memory)])
                upper <- max(memory[length(memory)-1], memory[length(memory)])
                diff <- upper - lower
                while(diff >= gridN){
                  if(Fun(lower+diff/2, Tn_1) < u){
                    lower <- lower+diff/2
                  }else{
                    upper <- lower+diff/2
                  }
                  diff <- upper - lower
                }
                Tn_1 + (lower+upper)/2
              }
              
              result <- matrix(0, sample.lengthN, n)  # here: N_t
              if(which.series == "current") result[, 1] <- startN
              for(i in 1:sample.lengthN){
                
                times <- drawTn(Tstart)
                while(times[length(times)] < t[n]){
                  times <- c(times, drawTn(times[length(times)]))
                }
                result[i, ] <- result[i, 1] + TimestoN(times, t)
              }
              Npred <- result
              
              dN <- t(apply(Npred, 1, diff))
              ## jump diffusion:
              
              dt <- diff(t)
              samples <- list(phi = object@phi[ind], thetaT = object@thetaT[ind], gamma2 = object@gamma2[ind])
              
              if(missing(sample.length)) sample.length <- K
              
              if(pred.alg == "Distribution"){
                if(only.interval){
                  result <- matrix(0, 2, n-1)
                  
                  for(i in 1:(n-1)){
                    VFun <- function(yn, yn_1, samples)  mean(pnorm(log(yn), log(yn_1)+sum(dt[1:i])*(samples$phi+ samples$gamma2/2) + samples$thetaT*(Npred[,i+1] - Npred[,1]), sqrt(samples$gamma2*sum(dt[1:i]))))
                    cand.Area <- exp(log(y.start) + sum(dt[1:i])*(mean(samples$phi)-mean(samples$gamma2)/2) + mean(samples$thetaT)*range(Npred[,i+1] - Npred[,1]) + 6* c(-1,1) *sqrt(mean(samples$gamma2)*sum(dt[1:i])))
                    
                    d <- diff(cand.Area)/cand.length
                    
                    result[,i] <- prediction.intervals(samples, VFun, x0 = y.start, level = level, candArea = cand.Area, grid = d)
                  }
                  
                }else{
                  result <- matrix(0, sample.length, n-1)
                  
                  for(i in 1:(n-1)){
                    VFun <- function(yn, yn_1, samples)  mean(pnorm(log(yn), log(yn_1)+sum(dt[1:i])*(samples$phi+ samples$gamma2/2) + samples$thetaT*(Npred[,i+1] - Npred[,1]), sqrt(samples$gamma2*sum(dt[1:i]))))
                    cand.Area <- exp(log(y.start) + sum(dt[1:i])*(mean(samples$phi)-mean(samples$gamma2)/2) + mean(samples$thetaT)*range(Npred[,i+1] - Npred[,1]) + 5* c(-1,1) *sqrt(mean(samples$gamma2)*sum(dt[1:i])))
                    
                    d <- diff(cand.Area)/cand.length
                    
                    result[,i] <- pred.base(samples = samples, VFun, x0 = y.start, len = sample.length, method = "vector",
                                            pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
                    if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
                  }
                  
                }


              }else{  # "Trajectory"
                VFun <- function(yn, yn_1, samples)  mean(pnorm(log(yn), log(yn_1)+dt[1]*(samples$phi+ samples$gamma2/2) + samples$thetaT*dN[,1], sqrt(samples$gamma2*dt[1])))
                cand.Area <- exp(log(y.start) + dt[1]*(mean(samples$phi)-mean(samples$gamma2)/2) + mean(samples$thetaT)*range(dN[,1]) + 5* c(-1,1) *sqrt(mean(samples$gamma2)*dt[1]))
                d <- diff(cand.Area)/cand.length
                
                result <- matrix(0, sample.length, n-1)
                result[,1] <- pred.base(samples = samples, VFun, x0 = y.start, len = sample.length, method = "free",
                                        pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
                
                for(i in 2:(n-1)){
                  VFun <- function(yn, yn_1, samples)  mean(pnorm(log(yn), log(yn_1)+dt[i]*(samples$phi+ samples$gamma2/2) + samples$thetaT*dN[,i], sqrt(samples$gamma2*dt[i])))
                  cand.Area <- exp(range(log(result[,i-1])) + dt[i]*(mean(samples$phi)-mean(samples$gamma2)/2) + mean(samples$thetaT)*range(dN[,i]) + 5* c(-1,1) *sqrt(mean(samples$gamma2)*dt[i]))
                  
                  d <- diff(cand.Area)/cand.length
                  
                  result[,i] <- pred.base(samples = samples, VFun, x0 = result[,i-1], len = sample.length, method = "free",
                                          pred.alg = pred.alg, sampling.alg = "BinSearch", candArea = cand.Area, grid = d)
                  if(i %% 10 == 0) message(paste(i, "of", n-1, "predictions are calculated"))
                }
                
              }
            }
            
            
            if(pred.alg == "simpleTrajectory"){
              if(missing(sample.length)) sample.length <- 100
              xi <- apply(object@xi[,ind], 1, mean)
              phi <- mean(object@phi[ind])
              gamma2 <- mean(object@gamma2[ind])
              thetaT <- mean(object@thetaT[ind])
              
              if(which.series == "new"){
                cl <- set.to.class("Merton", parameter = list(phi = phi, thetaT = thetaT, gamma2 = gamma2, xi = xi),
                                   Lambda = Lambda)
                result <- matrix(0, sample.length, n-1)
                Npred <- matrix(0, sample.length, n)  # here: N_t
                
                for(i in 1:sample.length){
                  help <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
                  result[i, ] <- help$Y[-1]
                  Npred[i, ] <- help$N
                }
                
              }else{
               result <- matrix(0, sample.length, n-1)
                Npred <- matrix(0, sample.length, n)  # here: N_t
                
                for(i in 1:sample.length){
                  if(length(startN) > 1){
                    ind.i <- sample(K, 1)
                    Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN[ind.i]), Lambda = Lambda)$N
                  }else{
                    ind.i <- 1
                    Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
                  }
                  
                  result[i, ] <- simY(t-max(object@t), phi, thetaT, gamma2, start = y.start, Npred[i,] - startN[ind.i])[-1]
                }
                
              }
            }
            
            if(pred.alg == "simpleBayesTrajectory"){
              if(missing(sample.length)){
                sample.length <- K
              } else{
                if(sample.length > K) sample.length <- K
              }
              
              xi <- object@xi[,ind]
              phi <- object@phi[ind]
              gamma2 <- object@gamma2[ind]
              thetaT <- object@thetaT[ind]
              
              if(which.series == "new"){
                result <- matrix(0, sample.length, n-1)
                Npred <- matrix(0, sample.length, n)  # here: N_t
                
                for(i in 1:sample.length){
                  cl <- set.to.class("Merton", parameter = list(phi = phi[i], thetaT = thetaT[i], gamma2 = gamma2[i], xi = xi[,i]),
                                     Lambda = Lambda)
                  help <- simulate(cl, t = t, y0 = y.start, plot.series = FALSE)
                  result[i, ] <- help$Y[-1]
                  Npred[i, ] <- help$N
                }
                
              }else{
                result <- matrix(0, sample.length, n-1)
                Npred <- matrix(0, sample.length, n)  # here: N_t
                
                for(i in 1:sample.length){
                  if(length(startN) > 1){
                    ind.i <- sample(K, 1)
                    Npred[i, ] <- simN(t, xi[,i], len = 1, start = c(max(object@t), startN[ind.i]), Lambda = Lambda)$N
                  }else{
                    ind.i <- 1
                    Npred[i, ] <- simN(t, xi[,i], len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
                  }
                  
                  result[i, ] <- simY(t-max(object@t), phi[i], thetaT[i], gamma2[i], start = y.start, Npred[i,] - startN[ind.i])[-1]
                }
                
              }
            }
            
            if(plot.prediction){
              qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
              if(which.series == "new"){
                plot(object@t, object@Y, type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
                lines(t[-1], qu[1,], col = 2)
                lines(t[-1], qu[2,], col = 2)
              }
              if(which.series == "current"){
                plot(object@t, object@Y, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
                lines(t[-1], qu[1,], col = 2)
                lines(t[-1], qu[2,], col = 2)
              }
            }
            return(list(N = Npred[, -1], Y = result))
            
          })



########
#' Prediction for regression model dependent on Poisson process
#'
#' @description Bayesian prediction of a regression model
#'   \eqn{y_i = f(t_i, N_i, \theta) + \epsilon_i}.
#' @param object class object of MCMC samples: "est.reg_hiddenNHPP"
#' @param t vector of time points to make predictions for
#' @param only.interval if TRUE: only calculation of prediction intervals
#' @param level level of the prediction intervals
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param Lambda.mat matrix-wise definition of intensity rate function (makes it faster)
#' @param fun.mat matrix-wise definition of regression function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param M2pred default 10, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector"), for jump diffusion
#' @param pred.alg prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleTrajectory"
#' @param sample.length length of samples to be drawn
#' @param grid fineness degress of approximation, for Poisson process
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' t <- seq(0,1, by = 0.01)
#' cl <- set.to.class("reg_hiddenNHPP", fun = function(t, N, theta) theta[1]*t + theta[2]*N, 
#'              parameter = list(theta = c(1,2), gamma2 = 0.1, xi = c(3, 1/4)), 
#'              Lambda = function(t, xi) (t/xi[2])^xi[1])
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' \dontrun{
#' pred <- predict(est, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1], 
#'                  fun.mat = function(t, N, theta) theta[,1]*t + theta[,2]*N)
#' }
#' pred <- predict(est, pred.alg = "simpleTrajectory", sample.length = 100)
#' @export
setMethod(f = "predict", signature = "est.reg_hiddenNHPP",
          definition = function(object, t, only.interval = TRUE, level = 0.05, burnIn, thinning, Lambda.mat, fun.mat, 
                                which.series = c("new", "current"), M2pred = 10,
                                cand.length = 1000, pred.alg = c("Distribution", "simpleTrajectory", "simpleBayesTrajectory"),
                                sample.length, grid = 1e-05, plot.prediction = TRUE) {
            
            pred.alg <- match.arg(pred.alg)
            which.series <- match.arg(which.series)
            
            if(missing(burnIn)) burnIn <- object@burnIn
            if(missing(thinning)) thinning <- object@thinning
            Lambda <- object@model$Lambda
            fun <- object@model$fun
            
            ind <- seq(burnIn + 1, nrow(object@xi), by = thinning)
            K <- length(ind)
            
            if(which.series == "new"){
              Tstart <- 0
              if(missing(t)) t <- object@t
            }
            
            if(which.series == "current"){
              Tstart <- max(object@t)
              if(missing(t)){
                dt <- median( diff(object@t))
                t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
              }
              if(length(object@N.est) == 0){  # -> object@N is equal to observed Poisson process
                startN <- max(object@N)
              }else{  # Poisson variables are filtered
                startN <- object@N.est[nrow(object@N.est), ind]
              }
            }
            
            n <- length(t)
            
            if(pred.alg == "Distribution"){
              
              samples <- object@xi[ind, ]
              
              if(missing(Lambda.mat)){
                Lambda.diff <- function(cand, Tn_1, samples){
                  sapply(1:K, function(i) Lambda(cand + Tn_1, samples[i,]) - Lambda(Tn_1, samples[i,]) )
                }
                Fun <- function(cand, Tn_1) 1 - mean(exp(- Lambda.diff(cand, Tn_1, samples)))
              }else{
                Fun <- function(cand, Tn_1){
                  1-mean(exp(-(Lambda.mat(cand + Tn_1, samples)-Lambda.mat(Tn_1, samples))))
                }
              }
              sample.lengthN <- K
              
              drawTn <- function(Tn_1){
                u <- runif(1)
                cand <- Tn_1 + 0.1
                memory <- cand
                
                while(length(unique(memory)) == length(memory)){
                  if(Fun(cand, Tn_1) < u){
                    cand <- cand*2
                    memory <- c(memory, cand)
                  }else{
                    cand <- cand/2
                    memory <- c(memory, cand)
                  }
                }
                lower <- min(memory[length(memory)-1], memory[length(memory)])
                upper <- max(memory[length(memory)-1], memory[length(memory)])
                diff <- upper - lower
                while(diff >= grid){
                  if(Fun(lower+diff/2, Tn_1) < u){
                    lower <- lower+diff/2
                  }else{
                    upper <- lower+diff/2
                  }
                  diff <- upper - lower
                }
                Tn_1 + (lower+upper)/2
              }
              
              result <- matrix(0, sample.lengthN, n)  # here: N_t
              if(which.series == "current") result[, 1] <- startN
              for(i in 1:sample.lengthN){
                
                times <- drawTn(Tstart)
                while(times[length(times)] < t[n]){
                  times <- c(times, drawTn(times[length(times)]))
                }
                result[i, ] <- result[i, 1] + TimestoN(times, t)
              }
              Npred <- result
              
              samples <- list(theta = as.matrix(object@theta[ind,]), gamma2 = object@gamma2[ind])
              sample.length <- K
              method <- "vector"

           
              if(only.interval){
                result <- matrix(0, 2, n)
                for(i in 1:n){
                  cand.Area <- range(fun(t[i], Npred[, i], apply(samples$theta, 2, mean))) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))
                  if(!missing(fun.mat)){
                    VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(t[i], Npred[,i], samples$theta), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                  }else{
                    VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(t[i], Npred[k,i], samples$theta[k,]), sqrt(samples$gamma2[k]))))
                  }
                  result[,i] <- prediction.intervals(samples, VFun, x0 = 0, level = level, candArea = cand.Area)
                }
                if(plot.prediction){
                  if(which.series == "new"){
                    plot(object@t, object@Y, xlab = "t", ylim = range(c(object@Y, range(result))), ylab = expression(Y[t]))
                    lines(t, result[1,], col = 2)
                    lines(t, result[2,], col = 2)
                  }
                  if(which.series == "current"){
                    plot(object@t, object@Y, xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(result))), xlab = "t", ylab = expression(Y[t]))
                    lines(t, result[1,], col = 2)
                    lines(t, result[2,], col = 2)
                  }
                }
              }else{
                result <- matrix(0, sample.length, n)
                for(i in 1:n){
                  
                  if(!missing(fun.mat)){
                    VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(t[i], Npred[,i], samples$theta), sqrt(samples$gamma2)))
                    dens <- function(yn, yn_1, samples)  mean(dnorm(yn, fun.mat(t[i], Npred[,i], samples$theta), sqrt(samples$gamma2)))
                  }else{
                    VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(t[i], Npred[k,i], samples$theta[k,]), sqrt(samples$gamma2[k]))))
                    dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, fun(t[i], Npred[k,i], samples$theta[k,]), sqrt(samples$gamma2[k]))))
                  }
                  
                  cand.Area <- range(fun(t[i], Npred[, i], apply(samples$theta, 2, mean))) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))
                  if(missing(grid)) d <- diff(cand.Area)/cand.length
                  
                  result[,i] <- pred.base(samples = samples, VFun, dens, x0 = 0, len = sample.length, method = method,
                                          pred.alg = "Distribution", sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
                }
                if(plot.prediction){
                  qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
                  if(which.series == "new"){
                    plot(object@t, object@Y, xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
                    lines(t, qu[1,], col = 2)
                    lines(t, qu[2,], col = 2)
                  }
                  if(which.series == "current"){
                    plot(object@t, object@Y, xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
                    lines(t, qu[1,], col = 2)
                    lines(t, qu[2,], col = 2)
                  }
                }
              }
              
              
            }else{
              
              if(pred.alg == "simpleTrajectory"){
                if(missing(sample.length)) sample.length <- 100
                xi <- apply(object@xi[ind, ], 2, mean)
                gamma2 <- mean(object@gamma2[ind])
                theta <- apply(object@theta[ind,], 2, mean)
                
                if(which.series == "new"){
                  cl <- set.to.class("reg_hiddenNHPP", parameter = list(theta = theta, gamma2 = gamma2, xi = xi), Lambda = Lambda, fun = fun)
                  result <- matrix(0, sample.length, n)
                  Npred <- matrix(0, sample.length, n)  # here: N_t
                  
                  for(i in 1:sample.length){
                    help <- simulate(cl, t = t, plot.series = FALSE)
                    result[i, ] <- help$Y
                    Npred[i, ] <- help$N
                  }
                  
                }else{
                  result <- matrix(0, sample.length, n)
                  Npred <- matrix(0, sample.length, n)  # here: N_t
                  
                  for(i in 1:sample.length){
                    if(length(startN) > 1){
                      Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN[sample(K, 1)]), Lambda = Lambda)$N
                    }else{
                      ind.i <- 1
                      Npred[i, ] <- simN(t, xi, len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
                    }
                    
                    result[i, ] <- fun(t, Npred[i,], theta) + rnorm(length(t), 0, sqrt(gamma2))
                  }
                }

              }
              
              if(pred.alg == "simpleBayesTrajectory"){
                if(missing(sample.length)){
                  sample.length <- K
                } else{
                  if(sample.length > K) sample.length <- K
                }
                xi <- object@xi[ind, ]
                gamma2 <- object@gamma2[ind]
                theta <- object@theta[ind,]
                
                if(which.series == "new"){
                  result <- matrix(0, sample.length, n)
                  Npred <- matrix(0, sample.length, n)  # here: N_t
                  
                  for(i in 1:sample.length){
                    cl <- set.to.class("reg_hiddenNHPP", parameter = list(theta = theta[i,], gamma2 = gamma2[i], xi = xi[i,]), Lambda = Lambda, fun = fun)
                    help <- simulate(cl, t = t, plot.series = FALSE)
                    result[i, ] <- help$Y
                    Npred[i, ] <- help$N
                  }
                  
                }else{
                  result <- matrix(0, sample.length, n)
                  Npred <- matrix(0, sample.length, n)  # here: N_t
                  
                  for(i in 1:sample.length){
                    if(length(startN) > 1){
                      Npred[i, ] <- simN(t, xi[i,], len = 1, start = c(max(object@t), startN[sample(K, 1)]), Lambda = Lambda)$N
                    }else{
                      ind.i <- 1
                      Npred[i, ] <- simN(t, xi[i,], len = 1, start = c(max(object@t), startN), Lambda = Lambda)$N
                    }
                    
                    result[i, ] <- fun(t, Npred[i,], theta[i,]) + rnorm(length(t), 0, sqrt(gamma2[i]))
                  }
                }
                
              }
              
              if(plot.prediction){
                qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
                if(which.series == "new"){
                  plot(object@t, object@Y, type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
                  lines(t, qu[1,], col = 2)
                  lines(t, qu[2,], col = 2)
                }
                if(which.series == "current"){
                  plot(object@t, object@Y, type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
                  lines(t, qu[1,], col = 2)
                  lines(t, qu[2,], col = 2)
                }
              }
              
            }
            
           return(list(N = Npred, Y = result))
            
          })


########
#' Prediction for regression model
#'
#' @description Bayesian prediction of regression model
#'   \eqn{y_i = f(\phi, t_i) + \epsilon_i}.
#' @param object class object of MCMC samples: "est.Regression"
#' @param t vector of time points to make predictions for
#' @param only.interval if TRUE: only calculation of prediction intervals
#' @param level level of the prediction intervals
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param method vectorial ("vector") or not ("free")
#' @param sampling.alg sampling algorithm, inversion method ("BinSearch") or rejection sampling ("RejSamp")
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' t <- seq(0,1, by = 0.01)
#' cl <- set.to.class("Regression", fun = function(phi, t) phi[1]*t + phi[2], 
#'                    parameter = list(phi = c(1,2), gamma2 = 0.1))
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data, 1000)
#' plot(est)
#' pred <- predict(est, fun.mat = function(phi, t) phi[,1]*t + phi[,2])
#' pred2 <- predict(est, fun.mat = function(phi, t) phi[,1]*t + phi[,2], only.interval = FALSE)
#' plot(density(pred2[,10]))
#' @export
setMethod(f = "predict", signature = "est.Regression",
          definition = function(object,
                                t, only.interval = TRUE, level = 0.05, burnIn, thinning,
                                fun.mat, which.series = c("new", "current"), M2pred = 10,
                                cand.length = 1000, method = c("vector", "free"),
                                sampling.alg = c("BinSearch", "RejSamp"), sample.length, grid, plot.prediction = TRUE) {
            
            method <- match.arg(method)
            sampling.alg <- match.arg(sampling.alg)
            which.series <- match.arg(which.series)
            
            if(missing(burnIn)) burnIn <- object@burnIn
            if(missing(thinning)) thinning <- object@thinning
            
            if(missing(t)){
              if(which.series == "new") t <- object@t
              if(which.series == "current"){
                dt <- median( diff(object@t))
                t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
              }
            }
            
            ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
            samples <- list(phi = as.matrix(object@phi[ind,]), gamma2 = object@gamma2[ind])
            K <- length(samples$gamma2)
            if(missing(sample.length)) sample.length <- K
            fun <- object@model$fun
            sT.fun <- object@model$sT.fun
            n <- length(t)

            if(only.interval){
              result <- matrix(0, 2, n)
              for(i in 1:n){
                cand.Area <- fun(apply(samples$phi, 2, mean), t[i]) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))*sT.fun(t[i])
                if(!missing(fun.mat)){
                  VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                }else{
                  VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                }
                result[,i] <- prediction.intervals(samples, VFun, x0 = 0, level = level, candArea = cand.Area)
              }
              if(plot.prediction){
                if(which.series == "new"){
                  plot(object@t, object@Y, xlab = "t", ylim = range(c(object@Y, range(result))), ylab = expression(Y[t]))
                  lines(t, result[1,], col = 2)
                  lines(t, result[2,], col = 2)
                }
                if(which.series == "current"){
                  plot(object@t, object@Y, xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(result))), xlab = "t", ylab = expression(Y[t]))
                  lines(t, result[1,], col = 2)
                  lines(t, result[2,], col = 2)
                }
              }
            }else{
              result <- matrix(0, sample.length, n)
              for(i in 1:n){
                
                if(!missing(fun.mat)){
                  VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                  dens <- function(yn, yn_1, samples)  mean(dnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                }else{
                  VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                  dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                }
                
                cand.Area <- fun(apply(samples$phi, 2, mean), t[i]) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))*sT.fun(t[i])
                if(missing(grid)) d <- diff(cand.Area)/cand.length
                
                result[,i] <- pred.base(samples = samples, VFun, dens, x0 = 0, len = sample.length, method = method,
                                        pred.alg = "Distribution", sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
              }
              if(plot.prediction){
                qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
                if(which.series == "new"){
                  plot(object@t, object@Y, xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
                  lines(t, qu[1,], col = 2)
                  lines(t, qu[2,], col = 2)
                }
                if(which.series == "current"){
                  plot(object@t, object@Y, xlim = range(c(object@t, t)), ylim = range(c(object@Y, range(qu))), xlab = "t", ylab = expression(Y[t]))
                  lines(t, qu[1,], col = 2)
                  lines(t, qu[2,], col = 2)
                }
              }
            }
            

            return(result)
            
          })



########
#' Prediction for mixed regression model
#'
#' @description Bayesian prediction of the regression model
#'   \eqn{y_i = f(\phi_j, t_i) + \epsilon_i, \phi_j~N(\mu, \Omega)}.
#' @param object class object of MCMC samples: "est.mixedRegression"
#' @param t vector of time points to make predictions for
#' @param only.interval if TRUE: only calculation of prediction intervals
#' @param level level of the prediction intervals
#' @param burnIn burn-in period
#' @param thinning thinning rate
#' @param fun.mat matrix-wise definition of drift function (makes it faster)
#' @param which.series which series to be predicted, new one ("new") or further development of current one ("current")
#' @param ind.pred index of series to be predicted, optional, if which.series = "current" and ind.pred missing, the last series is taken
#' @param M2pred optional, if current series to predicted and t missing, M2pred variables will be predicted with dt of observation time points
#' @param cand.length length of candidate samples (if method = "vector")
#' @param sample.length length of samples to be drawn
#' @param grid fineness degree of approximation
#' @param plot.prediction if TRUE, result are plotted
#'
#' @examples
#' mu <- c(10, 5); Omega <- c(0.9, 0.01)
#' phi <- cbind(rnorm(21, mu[1], sqrt(Omega[1])), rnorm(21, mu[2], sqrt(Omega[2])))
#' cl <- set.to.class("mixedRegression", 
#'          parameter = list(phi = phi, mu = mu, Omega = Omega, gamma2 = 0.1), 
#'          fun = function(phi, t) phi[1]*t + phi[2], sT.fun = function(t) 1)
#' t <- seq(0, 1, by = 0.01)
#' data <- simulate(cl, t = t, plot.series = TRUE)
#' est <- estimate(cl, t, data[1:20,], 2000)
#' plot(est)
#' pred <- predict(est, fun.mat = function(phi, t) phi[,1]*t + phi[,2])
#' points(t, data[21,], pch = 20)
#'
#' @export
setMethod(f = "predict", signature = "est.mixedRegression",
          definition = function(object, t, only.interval = TRUE, level = 0.05, burnIn, thinning,
                                fun.mat, which.series = c("new", "current"), ind.pred, M2pred = 10,
                                cand.length = 1000, sample.length, grid, plot.prediction = TRUE) {
            which.series <- match.arg(which.series)
            
            J <- nrow(object@Y)
            n <- ncol(object@Y)  # equal to length(object@t)
            
            if(missing(ind.pred) & which.series == "current") ind.pred <- J

            if(missing(burnIn)) burnIn <- object@burnIn
            if(missing(thinning)) thinning <- object@thinning
            
            if(missing(t)){
              if(which.series == "new") t <- object@t
              if(which.series == "current"){
                dt <- median( diff(object@t))
                t <- object@t[length(object@t)] + cumsum(c(0, rep(dt, M2pred)))
              }
            }
            ind <- seq(burnIn + 1, length(object@gamma2), by = thinning)
            
            if(which.series == "new"){
              samples <- list(mu = as.matrix(object@mu[ind,]), Omega = as.matrix(object@Omega[ind,]) )
              phi.pred <- predPhi(samples)
              
              samples <- list(phi = phi.pred, gamma2 = object@gamma2[ind])
            }
            if(which.series == "current"){
              samples <- list(phi = sapply(1:ncol(object@mu) , function(i) phi_ij(object@phi[ind], ind.pred, i) ), gamma2 = object@gamma2[ind])
            }
            K <- length(samples$gamma2)
            if(missing(sample.length)) sample.length <- K
            fun <- object@model$fun
            sT.fun <- object@model$sT.fun
            n <- length(t)

            if(only.interval){
              
              result <- matrix(0, 2, n)
              for(i in 1:n){
                cand.Area <- fun(apply(samples$phi, 2, mean), t[i]) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))*sT.fun(t[i])
                if(!missing(fun.mat)){
                  VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                }else{
                  VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                }
                result[,i] <- prediction.intervals(samples, VFun, x0 = 0, level = level, candArea = cand.Area)
              }
              if(plot.prediction){
                if(which.series == "new"){
                  plot(object@t, object@Y[1,], type = "l", xlab = "t", ylim = range(c(object@Y, range(result))), ylab = expression(Y[t]))
                  for(j in 2:J) lines(object@t, object@Y[j,])
                  lines(t, result[1,], col = 2)
                  lines(t, result[2,], col = 2)
                }
                if(which.series == "current"){
                  plot(object@t, object@Y[ind.pred, ], type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y[ind.pred, ], range(result))), xlab = "t", ylab = expression(Y[t]))
                  lines(t, result[1,], col = 2)
                  lines(t, result[2,], col = 2)
                }
              }

            }else{
              result <- matrix(0, sample.length, n)
              for(i in 1:n){
                
                if(!missing(fun.mat)){
                  VFun <- function(yn, yn_1, samples)  mean(pnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                  dens <- function(yn, yn_1, samples)  mean(dnorm(yn, fun.mat(samples$phi, t[i]), sqrt(samples$gamma2*sT.fun(t[i])^2)))
                }else{
                  VFun <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) pnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                  dens <- function(yn, yn_1, samples)  mean(sapply(1:K, function(k) dnorm(yn, fun(samples$phi[k,], t[i]), sqrt(samples$gamma2[k]*sT.fun(t[i])^2))))
                }
                
                cand.Area <- fun(apply(samples$phi, 2, mean), t[i]) + 6*c(-1, 1)*sqrt(mean(samples$gamma2))*sT.fun(t[i])
                if(missing(grid)) d <- diff(cand.Area)/cand.length
                
                result[,i] <- pred.base(samples = samples, VFun, dens, x0 = 0, len = sample.length, method = method,
                                        pred.alg = "Distribution", sampling.alg = sampling.alg, candArea = cand.Area, grid = d)
              }
              if(plot.prediction){
                qu <- apply(result, 2, quantile, c(level / 2, 1 - level / 2))
                if(which.series == "new"){
                  plot(object@t, object@Y[1,], type = "l", xlab = "t", ylim = range(c(object@Y, range(qu))), ylab = expression(Y[t]))
                  for(j in 2:J) lines(object@t, object@Y[j,])
                  lines(t, qu[1,], col = 2, lwd = 2)
                  lines(t, qu[2,], col = 2, lwd = 2)
                }
                if(which.series == "current"){
                  plot(object@t, object@Y[ind.pred, ], type = "l", xlim = range(c(object@t, t)), ylim = range(c(object@Y[ind.pred, ], range(qu))), xlab = "t", ylab = expression(Y[t]))
                  lines(t, qu[1,], col = 2)
                  lines(t, qu[2,], col = 2)
                }
              }
              
            }
            
           if(which.series == "new") return(list(Y = result, phi = phi.pred))
            if(which.series == "current") return(result)
            
          })


#####################
#####################

#' Prediction function
#'
#' @description Bayesian prediction of the parameter of a stochastic process
#'   \eqn{dY_t = b(\phi, t, Y_t)dt + s(\gamma, t, Y_t)dW_t + h(\eta, t, Y_t)dN_t}.
#' @param samples MCMC samples
#' @param VFun cumulative distribution function
#' @param dens density function
#' @param len number of samples to be drawn
#' @param x0 starting point
#' @param method vectorial ("vector") or not ("free")
#' @param pred.alg prediction algorithm, "Distribution" or "Trajectory"
#' @param sampling.alg sampling algorithm, rejection sampling ("RejSamp") or inversion method ("BinSearch")
#' @param candArea candidate area
#' @param grid fineness degree
#'
#' @return
#' vector of samples from prediction
#' @export
pred.base <- function(samples, VFun, dens, len = 100, x0, method = c("vector", "free"),
                 pred.alg = c("Distribution", "Trajectory"), sampling.alg = c("RejSamp", "BinSearch"),
                 candArea, grid = 0.001){
  # without a candidate area: currently only sampling of positive values...
  method <- match.arg(method)
  pred.alg <- match.arg(pred.alg)
  sampling.alg <- match.arg(sampling.alg)
  if(missing(candArea)){
    candArea <- findCandidateArea(VFun)
    if(method == "vector") candArea <- seq(candArea[1], candArea[2], by = grid)
  }else{
    if(method == "vector" & length(candArea) == 2) candArea <- seq(candArea[1], candArea[2], by = grid)
    if(method == "free" & length(candArea) > 2) candArea <- c(min(candArea), max(candArea))
  }

  if(method == "vector" & pred.alg == "Trajectory") message("does not make much sense")


  if(pred.alg == "Distribution"){
    if(length(x0) == 1) x0 <- rep(x0, len)
    if(sampling.alg == "RejSamp"){
      sample.result <- RejSampling(Fun = function(t) VFun(t, x0, samples), dens = function(t) dens(t, x0, samples),
                  len = len, cand = candArea, grid = grid, method = method)
    }
    if(sampling.alg == "BinSearch"){
      sample.result <- BinSearch(Fun = function(t) VFun(t, x0, samples), len = len,
                                 candArea = candArea, grid = grid, method = method)
    }

  }
  if(pred.alg=="Trajectory"){
    if(length(x0) == 1) x0 <- rep(x0, len)
    if(sampling.alg == "RejSamp"){
      sample.result <- sapply(1:len, function(i) RejSampling(Fun = function(t) VFun(t, x0[i], samples), dens = function(t) dens(t, x0, samples),
                                                             len = 1, cand = candArea, grid = grid, method = method) )
    }
    if(sampling.alg == "BinSearch"){
      sample.result <- sapply(1:len, function(i) BinSearch(Fun = function(t) VFun(t, x0[i], samples), len = 1,
                                                           candArea = candArea, grid = grid, method = method))
    }

  }

  sample.result
}


#' Prediction interval function
#'
#' @description Bayesian prediction of the parameter of a stochastic process
#'   \eqn{dY_t = b(phi, t, Y_t)dt + s(\gamma, t, Y_t)dW_t + h(\eta, t, Y_t)dN_t}.
#' @param samples MCMC samples
#' @param Fun cumulative distribution function
#' @param x0 starting point
#' @param level level of prediction intervals
#' @param candArea candidate area
#' @param grid fineness degree
#'
#' @return
#' prediction interval
#' @export
prediction.intervals <- function(samples, Fun, x0, level = 0.05, candArea, grid = 0.001){
  # without a candidate area: currently only sampling of positive values...
  if(missing(candArea)) candArea <- findCandidateArea(Fun)

  if(length(candArea) > 2){
    method <- "vector"
  }else{
    method <- "free"
  }
  sampling <- function(Fun, qu){
    if(method == "vector"){
      diFu <- sapply(candArea, Fun)
      result <- sapply(qu, function(u) candArea[which(diFu >= u)[1]])
    }
    if(method == "free"){
      n <- length(qu)
      result <- numeric(n)
      for(i in 1:n){
        u <- qu[i]
        lower <- candArea[1]
        upper <- candArea[2]
        diff <- upper - lower
        while(diff >= grid){
          if(Fun(lower+diff/2) < u){
            lower <- lower+diff/2
          }else{
            upper <- lower+diff/2
          }
          diff <- upper - lower
        }
        result[i] <- (lower+upper)/2
      }
    }
    result
  }

  sampling(Fun = function(t) Fun(t, x0, samples), qu = c(level/2, 1-level/2))
}


#' Prediction Of The Random Effects In Mixed Stochastic Differential Equations
#'
#' @description Prediction of the random effects \eqn{\phi ~ N(\mu, \Omega)}
#' @param samples output of estSDE or estREg
#' @param cand candidates for phi (matrix with p columns)
#' @return matrix phi
#'
#'
predPhi <- function(samples, cand){

  K <- nrow(samples$mu)
  p <- ncol(samples$mu)
  if(missing(cand)){
    cand <- sapply(1:p, function(i){
      seq(mean(samples$mu[,i])-4*sqrt(mean(samples$Omega[,i])), mean(samples$mu[,i])+4*sqrt(mean(samples$Omega[,i])), length = 1000)
    })
  }
  phi.res <- matrix(0, K, p)
  for(i in 1:p){
    prob <- numeric(length(cand[,i]))
    for(a in 1:length(cand[,i])){
      prob[a] <- mean(dnorm(cand[a,i], samples$mu[,i], sqrt(samples$Omega[,i])))
    }
    cand2 <- cand[prob != 0, i]
    prob <- prob[prob != 0]
    mp <- max(prob)
    count <- 1
    while(count <= K){
      u <- runif(1, 0, mp)
      ind <- sample(1:length(cand2), 1)
      if(u <= prob[ind]){
        phi.res[count, i] <- cand2[ind]
        count <- count + 1
      }
    }
  }
  return(phi.res)
}
