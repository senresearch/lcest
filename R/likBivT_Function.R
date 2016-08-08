
#' Likelihood Function for T Left Censored Data
#' 
#' This function is the maximum likelihood function of bivariate T distributed, left censored data.
#' @param thetaT transformed vector of parameters to optimize such that c(xmu, ymu, log(xsig), log(ysig), tanh(r)).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return calculated maximum likelihood function for optimazation using optim.
#' 
#Note: thetaIS denotes inverse of transformed theta and as scale matrix elements
#  ie theta = c(xmu, ymu, sig11, sig12, sig22)
likBivT <- function(thetaT, cenData, df){
  # Transform Theta
  thetaI <- transformThetaInv(thetaT)
  thetaIS <- convertThetaToScale(thetaI, df)
  #create vectors for each data point for any censoring scenario
  bothUncen <- bothUncen.T(thetaIS, cenData, df)
  bothCen <- bothCen.T(thetaIS, cenData, df)
  yCenOnly <- yCenOnly.T(thetaIS, cenData, df)
  xCenOnly <- xCenOnly.T(thetaIS, cenData, df)
  #creating a vector using indexing with flags to choose correct scenario
  sumVec = bothUncen*(1-cenData[,3])*(1-cenData[,4]) + bothCen*cenData[,3]*cenData[,4] + 
    xCenOnly*cenData[,3]*(1-cenData[,4]) +yCenOnly*cenData[,4]*(1-cenData[,3])
  #return sum
  sum(sumVec)
}

#' Log Likelihood Function Both Uncensored T
#' 
#' This is the log likelihood function when both x and y are uncensored for bivariate T distributed data.
#' @param thetaIS vector of parameters to optimize, as scale parameters and untransformed.
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return Vector whose sum would be log likelihood function if completely uncensored.
#'
bothUncen.T <- function(thetaIS, cenData, df){
  #Declare
  v <- vector( 'numeric', nrow(cenData))
  #Fill
  locVec <- c(thetaIS[1],thetaIS[2])
  scaleMat <-  matrix(c(thetaIS[3], thetaIS[5], thetaIS[5], thetaIS[4]), nrow=2, ncol=2 )
  dataMat <- cenData[,1:2]
  v <- dmt(dataMat, mean = locVec, S = scaleMat, df = df, log = TRUE)
  #return
  v
}


#' Log Likelihood Function Both Censored T
#' 
#' This is the portion of the log likelihood function when both x and y are censored for bivariate T distributed data.
#' @param thetaIS vector of parameters to optimize, as scale parameters and untransformed.
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return Vector whose to be summed for for this portion of the log liklihood function.
#'
bothCen.T <- function(thetaIS, cenData, df){
  #Declare
  v <- vector( 'numeric', nrow(cenData))
  #Fill
  locVec <- c(thetaIS[1],thetaIS[2])
  scaleMat <- matrix(c(thetaIS[3], thetaIS[5], thetaIS[5], thetaIS[4]), nrow=2, ncol=2 )
  v <- log(pmt(cenData[,1:2], mean = locVec, S = scaleMat, df = df) )
  #return
  v
}

#' Log Likelihood Function y Only Censored T
#' 
#' This is the portion of the log likelihood function when the y variable is cesnored and
#'    the x varaible isn't censored for bivariate t distributed data.
#' @param thetaIS vector of parameters to optimize, as scale parameters and untransformed.
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return Vector whose to be summed for this portion of the log liklihood function.
#'
yCenOnly.T <- function(thetaIS, cenData, df){
  #seperate data columns
  dataVecCol1 <- cenData[,1]
  dataVecCol2 <- cenData[,2]
  # Parameters for y|x
  p = 1 # because bivariate
  #  Mu
  mu <- thetaIS[2]+ thetaIS[5]/thetaIS[3]*(dataVecCol1 - thetaIS[1])
  #  Sigma
  d <- ((dataVecCol1-thetaIS[1])^2)/thetaIS[3]
  sigma22Given1  <- thetaIS[4] - (thetaIS[5]^2)/thetaIS[3]
  s  <- sigma22Given1*(df+d)/(df+p)
  # degrees of Freedom
  cdf <- df +p
  #Declare
  v <- vector( 'numeric', nrow(cenData))
  #Fill
  v <-dt((dataVecCol1 - thetaIS[1])/sqrt(thetaIS[3]), df, log = TRUE)-log(sqrt(thetaIS[3])) + #Pr(X=x) 
    pt((dataVecCol2 - mu)/sqrt(s), cdf, log.p = TRUE) #Pr(Y<=y|X=x)
  #return
  v
}

#' Log Likelihood Function x Only Censored T
#' 
#' This is the portion of the log likelihood function when the x variable is cesnored and
#'    the y varaible isn't censored for bivariate t distributed data.
#' @param thetaIS vector of parameters to optimize, as scale parameters and untransformed.
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return Vector whose to be summed for this portion of the log liklihood function.
#'
xCenOnly.T <- function(thetaIS, cenData, df){
  #seperate data columns
  dataVecCol1 <- cenData[,1]
  dataVecCol2 <- cenData[,2]
  # Parameters for x|y
  p = 1 # because bivariate
  #  Mu
  mu <- thetaIS[1]+ thetaIS[5]/thetaIS[4]*(dataVecCol2 - thetaIS[2])
  #  Sigma
  d <- ((dataVecCol1-thetaIS[2])^2)/thetaIS[4]
  sigma22Given1  <- thetaIS[3] - (thetaIS[5]^2)/thetaIS[4]
  s  <- sigma22Given1*(df+d)/(df+p)
  # degrees of Freedom
  cdf <- df +p
  #Declare
  v <- vector( 'numeric', nrow(cenData))
  #Fill
  v <-dt((dataVecCol2 - thetaIS[2])/sqrt(thetaIS[4]), df, log = TRUE)-log(sqrt(thetaIS[4])) + #Pr(Y=y) 
    pt((dataVecCol1 - mu)/sqrt(s), cdf, log.p = TRUE) #Pr(X<=y|Y=y)
  #return
  v
}
