% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/helper_methods.R
\docType{class}
\name{newmodel-class}
\alias{newmodel-class}
\alias{show,newmodel-method}
\alias{numberSamples,newmodel-method}
\alias{numberFeatures,newmodel-method}
\alias{numberFactors,newmodel-method}
\alias{newX,newmodel-method}
\alias{newV,newmodel-method}
\alias{newLogMu,newmodel-method}
\alias{newMu,newmodel-method}
\alias{newZeta,newmodel-method}
\alias{newPhi,newmodel-method}
\alias{newTheta,newmodel-method}
\alias{newEpsilon_beta,newmodel-method}
\alias{newEpsilon_gamma,newmodel-method}
\alias{newEpsilon_W,newmodel-method}
\alias{newEpsilon_alpha,newmodel-method}
\alias{newEpsilon_zeta,newmodel-method}
\alias{newW,newmodel-method}
\alias{newBeta,newmodel-method}
\alias{newGamma,newmodel-method}
\alias{newAlpha,newmodel-method}
\title{Class newmodel}
\usage{
\S4method{show}{newmodel}(object)

\S4method{numberSamples}{newmodel}(x)

\S4method{numberFeatures}{newmodel}(x)

\S4method{numberFactors}{newmodel}(x)

\S4method{newX}{newmodel}(object)

\S4method{newV}{newmodel}(object)

\S4method{newLogMu}{newmodel}(object)

\S4method{newMu}{newmodel}(object)

\S4method{newZeta}{newmodel}(object)

\S4method{newPhi}{newmodel}(object)

\S4method{newTheta}{newmodel}(object)

\S4method{newEpsilon_beta}{newmodel}(object)

\S4method{newEpsilon_gamma}{newmodel}(object)

\S4method{newEpsilon_W}{newmodel}(object)

\S4method{newEpsilon_alpha}{newmodel}(object)

\S4method{newEpsilon_zeta}{newmodel}(object)

\S4method{newW}{newmodel}(object)

\S4method{newBeta}{newmodel}(object)

\S4method{newGamma}{newmodel}(object)

\S4method{newAlpha}{newmodel}(object)
}
\arguments{
\item{object}{an object of class \code{newmodel}.}

\item{x}{an object of class \code{newmodel}.}
}
\value{
\code{numberSamples} returns the number of samples;
 \code{numberFeatures}returns the number of features; 
 \code{numberFactors} returns the number of latent factors.
}
\description{
Objects of this class store all the values needed to work with a
negative binomial model, as described in the vignette.
They contain all information to fit a model by penalized maximum likelihood
or simulate data from a model.
}
\details{
For the full description of the model see the model vignette.
  Internally, the slots are checked so that the matrices are of the
  appropriate dimensions: in particular, \code{X}, \code{O}
  and \code{W} need to have \code{n} rows, \code{V} needs to have \code{J}
  rows, \code{zeta} must be of length \code{J}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{show}: show useful info on the object.

\item \code{numberSamples}: returns the number of samples.

\item \code{numberFeatures}: returns the number of features.

\item \code{numberFactors}: returns the number of latent factors.

\item \code{newX}: returns the sample-level design matrix for mu.

\item \code{newV}: returns the gene-level design matrix for mu.

\item \code{newLogMu}: returns the logarithm of the mean of the non-zero
component.

\item \code{newMu}: returns the mean of the non-zero component.

\item \code{newZeta}: returns the log of the inverse of the dispersion
parameter.

\item \code{newPhi}: returns the dispersion parameter.

\item \code{newTheta}: returns the inverse of the dispersion parameter.

\item \code{newEpsilon_beta}: returns the regularization parameters for
\code{beta}.

\item \code{newEpsilon_gamma}: returns the regularization parameters for
\code{gamma}.

\item \code{newEpsilon_W}: returns the regularization parameters for
\code{W}.

\item \code{newEpsilon_alpha}: returns the regularization parameters for
\code{alpha}.

\item \code{newEpsilon_zeta}: returns the regularization parameters for
\code{zeta}.

\item \code{newW}: returns the matrix W of inferred sample-level
covariates.

\item \code{newBeta}: returns the matrix beta of inferred parameters.

\item \code{newGamma}: returns the matrix gamma of inferred parameters.

\item \code{newAlpha}: returns the matrix alpha of inferred parameters.
}}

\section{Slots}{

\describe{
\item{\code{X}}{matrix. The design matrix containing sample-level covariates, one
sample per row.}

\item{\code{V}}{matrix. The design matrix containing gene-level covariates, one gene
per row.}

\item{\code{X_intercept}}{logical. TRUE if X contains an intercept.}

\item{\code{V_intercept}}{logical. TRUE if V contains an intercept.}

\item{\code{W}}{matrix. The factors of sample-level latent factors.}

\item{\code{beta}}{matrix or NULL. The coefficients of X in the regression.}

\item{\code{gamma}}{matrix or NULL. The coefficients of V in the regression.}

\item{\code{alpha}}{matrix. The weight of sample-level latent factors.}

\item{\code{zeta}}{numeric. A vector of log of inverse dispersion parameters.}

\item{\code{epsilon_beta}}{nonnegative scalar. Regularization parameter for
beta}

\item{\code{epsilon_gamma}}{nonnegative scalar. Regularization parameter for
gamma}

\item{\code{epsilon_W}}{nonnegative scalar. Regularization parameter for W}

\item{\code{epsilon_alpha}}{nonnegative scalar. Regularization parameter for alpha}

\item{\code{epsilon_zeta}}{nonnegative scalar. Regularization parameter for zeta}
}}

