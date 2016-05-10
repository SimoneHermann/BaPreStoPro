\name{BaPreStoPro-package}
\alias{BaPreStoPro-package}
\alias{BaPreStoPro}
\docType{package}
\title{
Bayesian Prediction of Stochastic Processes
}
\description{
This package contains simulate, estimate and predict methods for non-homogeneous Poisson processes (NHPP), jump diffusions, (mixed) diffusions, hidden (mixed) diffusion models, regression model including a NHPP, and (mixed) regression models for comparison.
}
\details{
\tabular{ll}{
Package: \tab BaPreStoPro\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-05-10\cr
License: \tab GLP-2, GLP-3\cr
}
Each of the models has its specific class, \code{"jumpDiffusion", "Merton", "Diffusion", "mixedDiffusion", "hiddenDiffusion", "hiddenmixedDiffusion", "jumpRegression", "NHPP", "Regression", "mixedRegression"}, created with the function \code{set.to.class}. For each of the model classes, a method \code{simulate} and \code{estimate}. The output of method \code{estimate} is a new class object. For the estimation classes, methods \code{plot} and \code{predict} are available. 
}
\author{
Simone Hermann
<hermann@statistik.tu-dortmund.de>
}
\references{
S. Hermann, K. Ickstadt and C. H.  Mueller (2016). Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model. \emph{Applied Stochastic Models in Business and Industry 2016}.
Hermann, S. and F. Ruggeri (2016). Modelling Wear in Cylinder Liners. Under review
in: Quality and Reliability Engineering International.
Hermann, S., K. Ickstadt, and C. H. Mueller (2016). Bayesian Prediction for a Jump Diffusion
Process with Application to Crack Growth in Fatigue Experiments. Submitted to Reliability Engineering & System Safety.
Heeke, G., S. Hermann, R. Maurer, K. Ickstadt, and C. H. Mueller (2016). Stochastic Modeling
and Statistical Analysis of Fatigue Tests on Prestressed Concrete Beams under Cyclic Loadings.
Submitted to Structural Safety.
}
\keyword{ package }

\examples{
model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
t <- seq(0, 1, by = 0.01)
data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
est <- estimate(model, t, data, 100) # better: 10000
plot(est)
pred <- predict(est, b.fun.mat = function(phi, t, y) phi[,1])

}
