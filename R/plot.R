#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.jumpDiffusion class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.jumpDiffusion", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  if(any(dim(x@N.est) == 0)){
    if(nrow(x@xi) == 1)
      op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
                cex.axis = 0.7)
    if(nrow(x@xi) > 1)
      op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
                cex.axis = 0.7)

    plot(x@phi, type = "l", ylab = expression(phi)); abline(h = x@model$phi, col = 2)
    plot(x@theta, type = "l", ylab = expression(theta)); abline(h = x@model$theta, col = 2)
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    for(i in 1:nrow(x@xi)){
      plot(x@xi[i,], type = "l", ylab = bquote(xi[.(i)])); abline(h = x@model$xi[i], col = 2)
    }
  } else {
    op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
              cex.axis = 0.7)
    plot(x@phi, type = "l", ylab = expression(phi)); abline(h = x@model$phi, col = 2)
    plot(x@xi[1,], type = "l", ylab = expression(xi[1])); abline(h = x@model$xi[1], col = 2)
    plot(x@theta, type = "l", ylab = expression(theta)); abline(h = x@model$theta, col = 2)
    if(nrow(x@xi) > 1) plot(x@xi[2,], type = "l", ylab = expression(xi[2])); abline(h = x@model$xi[2], col = 2)
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    plot(x@t, x@N.est[,length(x@phi)], type = "l", xlab = "t", ylab = "N"); for(i in 1:length(x@phi)) lines(x@t, x@N.est[,i])
  }

  par(old.settings)
})


#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.Merton class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.Merton", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  if(any(dim(x@N.est) == 0)){
    if(nrow(x@xi) == 1)
      op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
                cex.axis = 0.7)
    if(nrow(x@xi) > 1)
      op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
                cex.axis = 0.7)

    plot(x@phi, type = "l", ylab = expression(phi)); abline(h = x@model$phi, col = 2)
    plot(x@thetaT, type = "l", ylab = expression(tilde(theta))); abline(h = x@model$thetaT, col = 2)
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    for(i in 1:nrow(x@xi)){
      plot(x@xi[i,], type = "l", ylab = bquote(xi[.(i)])); abline(h = x@model$xi[i], col = 2)
    } 
  } else {
    op <- par(mfrow = c(3, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
              cex.axis = 0.7)
    plot(x@phi, type = "l", ylab = expression(phi)); abline(h = x@model$phi, col = 2)
    plot(x@xi[1,], type = "l", ylab = expression(xi[1])); abline(h = x@model$xi[1], col = 2)
    plot(x@thetaT, type = "l", ylab = expression(tilde(theta))); abline(h = x@model$thetaT, col = 2)
    if(nrow(x@xi)>1) plot(x@xi[2,], type = "l", ylab = expression(xi[2])); abline(h = x@model$xi[2], col = 2)
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    plot(x@t[-1], x@N.est[,length(x@phi)], type = "l", xlab = "t", ylab = "N"); for(i in 1:length(x@phi)) lines(x@t[-1], x@N.est[,i])
  }

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.Diffusion class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.Diffusion", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@phi); if(p == 1) mfr <- c(2,1) else mfr <- c(2,2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@phi[,i], type = "l", ylab = bquote(phi[.(i)])); abline(h = x@model$phi[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.mixedDiffusion class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.mixedDiffusion", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@mu); mfr <- c(p+1,2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@mu[,i], type = "l", ylab = bquote(mu[.(i)])); abline(h = x@model$mu[i], col = 2)
    plot(x@Omega[,i], type = "l", ylab = bquote(omega[.(i)])); abline(h = x@model$Omega[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.hiddenDiffusion class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.hiddenDiffusion", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@phi); mfr <- c(ceiling((p+3)/2),2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@phi[,i], type = "l", ylab = bquote(phi[.(i)])); abline(h = x@model$phi[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
  plot(x@sigma2, type = "l", ylab = expression(sigma^2)); abline(h = x@model$sigma2, col = 2)

  plot(x@t, x@Z, pch = 20, ylab = "data")
  for(i in 1:nrow(x@phi)) lines(x@t, x@Y.est[i,], col = 3)
  points(x@t, x@Z, pch = 20)

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.hiddenmixedDiffusion class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.hiddenmixedDiffusion", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@mu); mfr <- c(p+2,2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@mu[,i], type = "l", ylab = bquote(mu[.(i)])); abline(h = x@model$mu[i], col = 2)
    plot(x@Omega[,i], type = "l", ylab = bquote(omega[.(i)])); abline(h = x@model$Omega[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
  plot(x@sigma2, type = "l", ylab = expression(sigma^2)); abline(h = x@model$sigma2, col = 2)

  plot(x@t, x@Z[1,], pch = 20, ylab = "first data series")
  for(i in 1:nrow(x@mu)) lines(x@t, x@Y.est[[i]][[1]], col = 3)
  points(x@t, x@Z[1,], pch = 20)

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.reg_hiddenNHPP class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.reg_hiddenNHPP", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)
  if(any(dim(x@N.est) == 0)){

    p <- ncol(x@theta); q <- ncol(x@xi); mfr <- c(ceiling(p+q)/2 + 1, 2)
    op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
              cex.axis = 0.7)
    for(i in 1:p) {
      plot(x@theta[,i], type = "l", ylab = bquote(theta[.(i)]))
      abline(h = x@model$theta[i], col = 2)
    }
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    for(i in 1:q){
      plot(x@xi[,i], type = "l", ylab = bquote(xi[.(i)]))
      abline(h = x@model$xi[i], col = 2)
    } 
  } else {
    p <- ncol(x@theta); q <- ncol(x@xi); mfr <- c(ceiling((p+q)/2) + 1, 2)
    op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
              cex.axis = 0.7)

    for(i in 1:p) {
      plot(x@theta[,i], type = "l", ylab = bquote(theta[.(i)])); abline(h = x@model$theta[i], col = 2)
    }
    plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)
    for(i in 1:q){
      plot(x@xi[,i], type = "l", ylab = bquote(xi[.(i)]))
      abline(h = x@model$xi[i], col = 2)
    } 
    plot(x@t, x@N.est[,length(x@gamma2)], type = "l", xlab = "t", ylab = "N"); for(i in 1:length(x@gamma2)) lines(x@t, x@N.est[,i])
  }

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.NHPP class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.NHPP", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)
  q <- ncol(x@xi); mfr <- c(ceiling(q)/2, 2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)

  for(i in 1:q){
    plot(x@xi[,i], type = "l", ylab = bquote(xi[.(i)]))
    abline(h = x@model$xi[i], col = 2)
  }
  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.Regression class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.Regression", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@phi); if(p == 1) mfr <- c(2,1) else mfr <- c(2,2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@phi[,i], type = "l", ylab = bquote(phi[.(i)])); abline(h = x@model$phi[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)

  par(old.settings)
})

#' Plot method for the Bayesian estimation class object
#' 
#' @description Plot method for the S4 class Bayes.fit
#' @param x est.mixedRegression class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#' @export
setMethod(f = "plot", signature = "est.mixedRegression", definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  old.settings <- par(no.readonly = TRUE)

  p <- ncol(x@mu); mfr <- c(p+1,2)
  op <- par(mfrow = mfr, mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7,
            cex.axis = 0.7)
  for(i in 1:p) {
    plot(x@mu[,i], type = "l", ylab = bquote(mu[.(i)])); abline(h = x@model$mu[i], col = 2)
    plot(x@Omega[,i], type = "l", ylab = bquote(omega[.(i)])); abline(h = x@model$Omega[i], col = 2)
  }
  plot(x@gamma2, type = "l", ylab = expression(gamma^2)); abline(h = x@model$gamma2, col = 2)

  par(old.settings)
})

