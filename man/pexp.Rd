% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extrafunctions.R
\name{pexp}
\alias{pexp}
\title{Power exponential covariance function}
\usage{
pexp(x)
}
\arguments{
\item{x}{A list with named elements `pars` and `data`. `pars` is a vector with one parameter value, and
`data` is the data `x`}
}
\value{
vector of values of the function
}
\description{
Power exponential covariance function
}
\details{
The function:

\dexp{f(x) = \theta_1^x}
}
\examples{
fexp(list(pars = c(0.8),data=runif(10)))
}
