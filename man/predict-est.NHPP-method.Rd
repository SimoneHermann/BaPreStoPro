% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/predict.R
\docType{methods}
\name{predict,est.NHPP-method}
\alias{predict,est.NHPP-method}
\title{Prediction for a non-homogeneous Poisson process}
\usage{
\S4method{predict}{est.NHPP}(object, variable = c("eventTimes",
  "PoissonProcess"), t, burnIn, thinning, Lambda.mat, which.series = c("new",
  "current"), Tstart, M2pred = 10, rangeN = c(0, 5), cand.length = 1000,
  pred.alg = c("Trajectory", "Distribution", "simpleTrajectory",
  "simpleBayesTrajectory"), sample.length, grid = 1e-05,
  plot.prediction = TRUE)
}
\arguments{
\item{object}{class object of MCMC samples: "est.NHPP", created with method \code{\link{estimate,NHPP-method}}}

\item{variable}{if prediction of event times ("eventTimes") or of Poisson process variables ("PoissonProcess")}

\item{t}{vector of time points to make predictions for (only for variable = "PoissonProcess")}

\item{burnIn}{burn-in period}

\item{thinning}{thinning rate}

\item{Lambda.mat}{matrix-wise definition of drift function (makes it faster)}

\item{which.series}{which series to be predicted, new one ("new") or further development of current one ("current")}

\item{Tstart}{optional, if missing, first (which.series = "new") or last observation variable ("current") is taken}

\item{M2pred}{optional, if current series to be predicted and t missing, \code{M2pred} variables will be predicted
with the observation time distances}

\item{rangeN}{vector of candidate area for differences of N, only if pred.alg = "Distribution" and variable = "PoissonProcess"}

\item{cand.length}{length of candidate samples (if method = "vector")}

\item{pred.alg}{prediction algorithm, "Distribution", "Trajectory", "simpleTrajectory" or "simpleBayesTrajectory"}

\item{sample.length}{number of samples to be drawn, default is the number of posterior samples}

\item{grid}{fineness degree of sampling approximation}

\item{plot.prediction}{if TRUE, prediction intervals are plotted}
}
\description{
Bayesian prediction of a non-homogeneous Poisson process with cumulative intensity function \eqn{\Lambda(t, \xi)}.
}
\examples{
model <- set.to.class("NHPP", parameter = list(xi = c(5, 1/2)),
               Lambda = function(t, xi) (t/xi[2])^xi[1])
t <- seq(0, 1, by = 0.01)
data <- simulate(model, t = t)
est <- estimate(model, t, data$Times, 1000)  # nMCMC should be much larger!
plot(est)
pred <- predict(est, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1],
   variable = "PoissonProcess", pred.alg = "Distribution")

\dontrun{
pred_NHPP <- predict(est, Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
pred_NHPP <- predict(est, variable = "PoissonProcess",
   Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
pred_NHPP2 <- predict(est, which.series = "current",
   Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
pred_NHPP3 <- predict(est, variable = "PoissonProcess", which.series = "current",
                      Lambda.mat = function(t, xi) (t/xi[,2])^xi[,1])
pred_NHPP4 <- predict(est, pred.alg = "simpleTrajectory", M2pred = length(data$Times))
}
pred_NHPP <- predict(est, variable = "PoissonProcess", pred.alg = "simpleTrajectory",
                     M2pred = length(data$Times))
pred_NHPP <- predict(est, variable = "PoissonProcess", pred.alg = "simpleBayesTrajectory",
                     M2pred = length(data$Times), sample.length = 100)
}
\references{
Hermann, S. (2016a). BaPreStoPro: an R Package for Bayesian Prediction of Stochastic Processes.
SFB 823 discussion paper 28/16.

Hermann, S. (2016b). Bayesian Prediction for Stochastic Processes based on the Euler Approximation Scheme.
SFB 823 discussion paper 27/16.
}

