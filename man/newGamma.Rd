% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{newGamma}
\alias{newGamma}
\title{Returns the matrix of paramters gamma}
\usage{
newGamma(object, ...)
}
\arguments{
\item{object}{an object that describes a matrix of negative binomial
distributions.}

\item{...}{Additional parameters.}
}
\value{
the matrix of gamma parameters
}
\description{
Given an object that describes a matrix of negative binomial distributions,
returns the matrix of parameters associated with V
}
\examples{
a <- newmodel(n=5, J=10)
newGamma(a)
}
