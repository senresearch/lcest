% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescaleVarCovMat_function.R
\name{rescaleVarCovMat}
\alias{rescaleVarCovMat}
\title{Rescale Varaince Covaraince Matrix}
\usage{
rescaleVarCovMat(varCovMat, transPar)
}
\arguments{
\item{varCovMat}{a 5x5 matrix which is the inverted hessian result of optim}

\item{transPar}{vector of length 5, the parameter results of optim still on transformed scale}
}
\value{
the varaiance covariance mat, now in orginal scale
}
\description{
The var cov mat is transformed when we get it from optim becuase optim works on the transformed scale. This function puts it back on the orginal or
inverse of transposed scale.
}

