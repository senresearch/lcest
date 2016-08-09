######################### likBivNorm Functions #######################
# Not public to Package User


#' Likelihood Function for Normal Left Censored Data
#' 
#' This function is the maximum likelihood function of bivariate, normally distributed, left censored data.
#' @param thetaT transformed vector of parameters to optimize such that c(xmu, ymu, log(xsig), log(ysig), tanh(r)).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @return calculated maximum likelihood function for optimazation using optim.
#' 
likBiv.N<- function(thetaT, cenData){
  # Transform Theta
  thetaI <- transformThetaInv(thetaT)
  #create vectors for each data point for any censoring scenario
  bothUncen <- bothUncen.N(thetaI, cenData[,1:2])
  bothCen <- bothCen.N(thetaI, cenData[,1:2])
  xOnlyCen <- xOnlyCen.N(thetaI, cenData[,1:2])
  yOnlyCen <- yOnlyCen.N(thetaI, cenData[,1:2])
  #creating a vector using binary cen vector to choose correct scenario
  sumVec = bothUncen*(1-cenData[,3])*(1-cenData[,4]) + bothCen*cenData[,3]*cenData[,4] +
    xOnlyCen*cenData[,3]*(1-cenData[,4]) + yOnlyCen*cenData[,4]*(1-cenData[,3])
  #return sum
  sum(sumVec)
}


#' Log Likelihood Function Both Uncensored Normal
#' 
#' This is the log likelihood function when both x and y are uncensored for bivariate, normally distributed data.
#' @param thetaI vector of parameters to optimize untransformed such that c( xmu, ymu, xsig, ysig, r).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @importFrom mnormt dmnorm
#' @return Vector whose sum would be log likelihood function if completely uncensored.
#' 
bothUncen.N <- function(thetaI, cenData){
  #initilaize
  v <- vector( "numeric", length(cenData[,1]))
  #declare
  ## make mu vec
  muVec = c(thetaI[1],thetaI[2]) 
  ## make covariance martix
  covMat <- matrix( c(thetaI[3]^2, thetaI[5]*thetaI[4]*thetaI[3],
                      thetaI[5]*thetaI[4]*thetaI[3], thetaI[4]^2),
                    nrow=2, ncol=2 )
  v <- dmnorm(cenData[,1:2], mean = muVec, varcov = covMat, log = TRUE)
  #return
  v
}

#' Log Likelihood Function Both Censored Normal
#' 
#' This is the portion of the log likelihood function when both x and y are censored for bivariate, normally distributed data.
#' @param thetaI vector of parameters to optimize untransformed such that c( xmu, ymu, xsig, ysig, r).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @return Vector whose to be summed for for this portion of the log liklihood function.
#'
bothCen.N <-function(thetaI, cenData){
  #Initialize
  v <- vector( "numeric", length(cenData[,1]))
  # make mean vec
  muVec <- c(thetaI[1], thetaI[2])
  # return log likelihood by using 2-dim distribution function
  #if no correlation                 
  #v <- pnorm(cenData[,1], muVec[1], sd = thetaI[3], log.p = TRUE) + 
  # pnorm(cenData[,2], muVec[2], sd = thetaI[4], log.p = TRUE) 
  # with correlation
  # designed by Dr Sen
  # no better results that ONE, just faster
  v <- pnorm2(cenData[,1:2], muVec, c(thetaI[3],thetaI[4]), thetaI[5] ,log.p=TRUE)
  #Return
  v
}

#' Log Likelihood Function y Only Censored Normal
#' 
#' This is the portion of the log likelihood function when the y variable is cesnored and
#'    the x varaible isn't censored for bivariate, normally distributed data.
#' @param thetaI vector of parameters to optimize untransformed such that c( xmu, ymu, xsig, ysig, r).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @return Vector whose to be summed for this portion of the log liklihood function.
#'
yOnlyCen.N <- function(thetaI, cenData){
  #Initialize
  v <- vector( "numeric", length(cenData[,1]))
  #Declare
  ## mean of y given x
  yGivenX <- thetaI[2] + thetaI[5]*thetaI[4]*(cenData[,1]-thetaI[1])/thetaI[3]
  ## conditional standard deviation of y
  sdYGivenX <- thetaI[4]*sqrt(1-thetaI[5]^2)
  #Return
  v <- pnorm(cenData[,2], yGivenX, sdYGivenX, log.p=TRUE) + dnorm(cenData[,1], thetaI[1], thetaI[3] ,log=TRUE)
}

#' Log Likelihood Function x Only Censored Normal
#' 
#' This is the portion of the log likelihood function when the x variable is cesnored and
#'    the y varaible isn't censored for bivariate, normally distributed data.
#' @param thetaI vector of parameters to optimize untransformed such that c( xmu, ymu, xsig, ysig, r).
#' @param cenData matrix of data in which first two columns hold data values for x and y and the second two
#'        columns hold flags for censoring such that 1 implies censored and 0 implies not censored.
#' @return Vector whose to be summed for this portion of the log liklihood function.
#'
xOnlyCen.N <- function(thetaI, cenData){
  #Initialize
  v <- vector( "numeric", length(cenData[,1]))
  #Declare
  ## mean of x given y
  xGivenY <- thetaI[1] + thetaI[5]*thetaI[3]*(cenData[,2]-thetaI[2])/thetaI[4]
  ## conditional standard deviation of x
  sdXGivenY <- thetaI[3]*sqrt(1-thetaI[5]^2)
  #Return
  v <- pnorm(cenData[,1], xGivenY, sdXGivenY, log.p=TRUE) + dnorm(cenData[,2], thetaI[2], thetaI[4] ,log=TRUE)
}