\name{BaPreStoPro-package}
\alias{BaPreStoPro-package}
\alias{BaPreStoPro}
\docType{package}
\title{
Bayesian Prediction of Stochastic Processes
}
\description{
This package contains simulate, estimate and predict methods for jump diffusions, diffusions, mixed diffusions, hidden (mixed) diffusion models and regression models for comparison.
}
\details{
\tabular{ll}{
Package: \tab BaPreStoPro\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-04-25\cr
License: \tab GLP-2, GLP-3\cr
}
 An overview of how to use the package, including the most important functions 
}
\author{
Simone Hermann
<hermann@statistik.tu-dortmund.de>
}
\references{

Bayesian Prediction of Crack Growth Based on a Hierarchical Diffusion Model.  S. Hermann, K. Ickstadt and C. Mueller, \emph{appearing in: Applied Stochastic Models in Business and Industry 2016}.

}
\keyword{ package }

\examples{
model <- set.to.class("Diffusion", parameter = list(phi = 0.5, gamma2 = 0.01))
t <- seq(0, 1, by = 0.01)
data <- simulate(model, t = t, y0 = 0.5, plot.series = TRUE)
est_diff <- estimate(model, t, data, 100) # better: 10000
plot(est_diff)
pred_diff <- predict(est_diff, b.fun.mat = function(phi, t, y) phi[,1])  # much faster

}
