% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{newBeta}
\alias{newBeta}
\title{Returns the matrix of paramters beta}
\usage{
newBeta(object, ...)
}
\arguments{
\item{object}{an object that describes a matrix of negative binomial
distributions.}

\item{...}{Additional parameters.}
}
\value{
the matrix of beta parameters
}
\description{
Given an object that describes a matrix of negative binomial distributions,
returns the matrix of parameters associated with X
}
\examples{
a <- newmodel(n=5, J=10)
newBeta(a)
}
