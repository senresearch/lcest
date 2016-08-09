

#' Rescale Varaince Covaraince Matrix
#' 
#' The var cov mat is transformed when we get it from optim becuase optim works on the transformed scale. This function puts it back on the orginal or
#' inverse of transposed scale.
#' @param thetaT transformed vector of parameters to optimize such that c(xmu, ymu, log(xsig), log(ysig), tanh(r)).
#' @param varCovMat a 5x5 matrix which is the inverted hessian result of optim
#' @param transPar vector of length 5, the parameter results of optim still on transformed scale
#' @return the varaiance covariance mat, now in orginal scale
#' 
rescaleVarCovMat <- function(varCovMat, transPar){
  diagMat <- rbind( c(1,0,0,0,0), c(0,1,0,0,0), c(0,0,exp(transPar[3]),0,0),
                    c(0,0,0,exp(transPar[4]),0), c(1,0,0,0, 1/(cosh(transPar[5])^2) ) )
  diagMatTranspose <- t(diagMat)
  varCovMat <- (diagMat%*% varCovMat) %*% diagMatTranspose
}
