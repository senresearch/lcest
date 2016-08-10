#' Bivaraite Left Censored Estimates
#' 
#' This function computes the parameter estimates using maximum likelihood estimate (MLE) method.
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom with df=Inf implying normal. Must be greater than 3 and an intger.
#' @param thetaG vector of length 5, a guess of approximate values for (xmu,ymu,xsd,ysd,r), which
#'        by default uses the means, standard deviations and correlations not adjusting for censoring.
#' @param alpha confidence levels will be 1-alpha level.
#' @param control list of parameters to pass to control in optim. See detail and optim documentation for more information.
#' @details The maximum likelihood method is done with the optim function. Therefore thetaG is the initial 
#'          values for the parameters to be optimized over. 
#'          The control parameter allow the user to change certian parameters in optim. The most helpful is likely maxit.
#'          Increaing maxit may increase run time, but will decrease the change of no convergence 
#'          from convergence error code 1. For more information on convergence error codes see optim documentation.
#'          Please do not change fnscale from equaling -1. This is neccessay for MLE method.
#' @return Add Details 
#' @export
#' @example 
#' 
#' 
blcest <-function(cenData, df=Inf, thetaG = defaultGuess(cenData), alpha =.05, control=list(fnscale=-1, maxit = 1000)){
  # Warning for change of fnscale
  if( control$fnscale == -1 ) warning( " To properly run blcest fnscale must equal -1. ")
  #transform Starting Guess
  thetaGT <- transformTheta(thetaG)
  # control=list(fnscale=-1) maximizes instead of optims default of min
  # Run Optim Normal
  if ( df == Inf){
    results <- optim( par = thetaGT, likBiv.N, cenData = cenData,
                      control=control, hessian = TRUE )
  }
  # Run optim T
  if (df > 3 && df !=Inf){
    if ( df %% 1 != 0) stop( "df must be integer")
    results <- optim( par = thetaGT, likBiv.T, cenData = cenData, df=df,
                      control=list(fnscale=-1, maxit = maxit), hessian = TRUE )
  }
  if (df < 3 ) stop( "df must be greater than 3 or equal to 0 for normal")
  #Converts back to the bounds we interprete result on
  estimate <- transformThetaInv(results$par)
  # Changing Hessian to CovVarMat
  varCovMatT <- qr.solve(-results$hessian)
  varCovMat<-rescaleVarCovMat(varCovMatT, results$par)
  
  rownames(varCovMat) <- c("xMu", "yMu", "xSd", "ySd", "R")
  colnames(varCovMat) <- c("xMu", "yMu", "xSd", "ySd", "R")
  # finding Standard Errors
  stdError <-  sqrt(abs(diag(varCovMat)))
  # finding CI
  z <- abs(qnorm(alpha/2))
  lowerCI <- estimate - stdError*z
  upperCI <- estimate + stdError*z
  # compute T value
  tValue <- estimate/stdError
  # Creating coefficients Output (like lm)
  coefficients <- data.frame( estimate, stdError, tValue, upperCI, lowerCI,
                              row.names= c("xMu", "yMu", "xSd", "ySd", "R") )
  # Warns if does not convrege 
  if( results$convergence != 0){
    warning("An optim does not converge. Optim's Converge Output is ", results$convergence)}
  # returns 
  list( coefficients = round (coefficients, 4) , varCovMatrix=signif(varCovMat,4))
}