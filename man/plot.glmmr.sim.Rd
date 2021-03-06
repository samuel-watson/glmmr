% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printfunctions.R
\name{plot.glmmr.sim}
\alias{plot.glmmr.sim}
\title{Plotting method for glmmr.sim}
\usage{
\method{plot}{glmmr.sim}(x, par, alpha = 0.05)
}
\arguments{
\item{x}{A `glmmr.sim` object}

\item{par}{Integer indicating the index of the parameter to plot if not a Bayesian analysis}

\item{alpha}{Numeric indictaing the type I error rate for non-Bayesian analysis}
}
\value{
A `ggplot2` plot
}
\description{
Plots a glmmr.sim object
}
\details{
For a Bayesian simulation analysis, this will plot the simulation based calibration ranks and the 
distribution of the posterior variance. Otherwise, it will plot the distribution of confidence interval 
half widths.
}
\examples{
...
}
