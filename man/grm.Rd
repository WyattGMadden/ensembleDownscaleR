% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grm.R
\name{grm}
\alias{grm}
\title{Fit the Geostatistical Regression Model (GRM)}
\usage{
grm(
  Y,
  X,
  L = NULL,
  M = NULL,
  coords,
  space.id,
  time.id,
  spacetime.id,
  include.additive.weekly.resid = T,
  include.additive.annual.resid = T,
  include.multiplicative.weekly.resid = T,
  include.multiplicative.annual.resid = T,
  num_neighbors = 10,
  n.iter = 25000,
  burn = 5000,
  thin = 4,
  covariance = "exponential",
  covariance.kernal = NULL,
  matern.nu = 1.5,
  tau.alpha.a = 0.5,
  tau.alpha.b = 0.005,
  tau.beta.a = 0.5,
  tau.beta.b = 0.005,
  omega.a = 0.5,
  omega.b = 0.005,
  theta.alpha.tune = 0.2,
  theta.alpha.a = 5,
  theta.alpha.b = 0.05,
  theta.alpha.init = 100,
  theta.beta.tune = 0.2,
  theta.beta.a = 5,
  theta.beta.b = 0.05,
  theta.beta.init = 100,
  rho.alpha.init = 0.9999,
  rho.beta.init = 0.9999,
  sigma.a = 0.001,
  sigma.b = 0.001,
  verbose = TRUE,
  verbose.iter = 1000
)
}
\arguments{
\item{Y}{Dependent variable vector (n)}

\item{X}{Unstandardized primary independent variable vector (n)}

\item{L}{Unstandardized spatial covariate matrix (n, p1)}

\item{M}{Unstandardized spatio-temporal covariate matrix (n, p2)}

\item{coords}{Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (n, 2)}

\item{space.id}{Spatial location ID vector (n)}

\item{time.id}{Temporal location ID vector (n)}

\item{spacetime.id}{ID vector of time points where spatial trends vary (n)}

\item{include.additive.weekly.resid}{Include additive temporal (weekly) residual bias not explained by other inputes}

\item{include.additive.annual.resid}{Include additive spatial (annual) residual bias not explained by other inputes}

\item{include.multiplicative.weekly.resid}{Include multiplicative temporal (weekly) residual bias not explained by other inputes}

\item{include.multiplicative.annual.resid}{Include multiplicative spatial (weekly) residual bias not explained by other inputes}

\item{num_neighbors}{Number of nearest neighbors to use in NNGP}

\item{n.iter}{Number of iterations used in predictions.}

\item{burn}{Number of pre-covergence simulations}

\item{thin}{Save every thin'th simulation}

\item{covariance}{Specify covariance function (from "exponential", "matern", "custom")}

\item{covariance.kernal}{Specify a custom covariance function if covariance = "custom". Must be a function with "distance" and "theta" parameters.}

\item{matern.nu}{Specify nu parameter for Matern covariance function if used (from 0.5, 1.5, and 2.5)}

\item{tau.alpha.a}{First tau alpha prior hyperparameter}

\item{tau.alpha.b}{Second tau alpha prior hyperparameter}

\item{tau.beta.a}{First tau beta prior hyperparameter}

\item{tau.beta.b}{Second tau beta prior hyperparameter}

\item{omega.a}{First omega prior hyperparameter}

\item{omega.b}{Second omega prior hyperparameter}

\item{theta.alpha.tune}{Theta alpha Metropolis-Hastings proposal tuning parameter}

\item{theta.alpha.a}{First theta alpha prior hyperparameter}

\item{theta.alpha.b}{Second theta alpha prior hyperparameter}

\item{theta.alpha.init}{Initial value for theta alpha}

\item{theta.beta.tune}{Theta beta Metropolis-Hastings proposal tuning parameter}

\item{theta.beta.a}{First theta beta prior hyperparameter}

\item{theta.beta.b}{Second theta beta prior hyperparameter}

\item{theta.beta.init}{Initial value for theta beta}

\item{rho.alpha.init}{Initial value for rho alpha}

\item{rho.beta.init}{Initial value for rho beta}

\item{sigma.a}{First sigma prior hyperparameter}

\item{sigma.b}{Second sigma prior hyperparameter}

\item{verbose}{Print MCMC output}

\item{verbose.iter}{print MCMC output step number each 'verbose.iter' iterations}
}
\value{
A list containing MCMC output
}
\description{
This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M
}
\examples{
# grm()


}
