


#' Defualt Guess
#' 
#' This function picksa reasonable starting guess for optim, which is the mean, sd, and corr of the data ignoring censoring.
#' @param thetaT transformed vector of parameters to optimize such that c(xmu, ymu, log(xsig), log(ysig), tanh(r)).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @return guesstheta vector of length 5, guess for each xmu, ymu, xsd, ysd, rho
#'

defualtGuess <- function(cenData){
  gxMu = mean(cenData[,1] , na.rm = TRUE)
  gyMu = mean(cenData[,2] , na.rm = TRUE)
  gxSd = sd(cenData[,1] , na.rm = TRUE)
  gySd = sd(cenData[,2] , na.rm = TRUE)
  corMat = cor( cenData[,1:2])
  gr = corMat[1,2]
  guessTheta <- c(gxMu, gyMu, gxSd, gySd, gr)
  # return starting guess
  guessTheta
}