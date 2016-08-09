

#' Generate Bivariate T Data
#' 
#' This function will generate an uncensored dataset distributed with bivaraite t
#' @param n number of pairs of data to generate, number of rows in data set, or sample size
#' @param locVec a vector of length d, representing the location parameter (equal to the mean vector when df>1)
#' @param scaleMat a symmetric positive-definite matrix representing the scale matrix of the distribution, such that S*df/(df-2) is the variance-covariance matrix when df>2
#' @param df Degrees of Freedom must be greater than 3 and an integer.
#' @return data set matrix with nrow = n and ncol= d
#' @export
#' @examples 
#'
genData.T <- function(n, locVec, scaleMat, df){
  uncenData <- rmt( n, locVec, scaleMat, df)
}


#' Generate Bivariate Norm Data
#' 
#' This function will generate an uncensored dataset distributed with bivaraite normal
#' @param n number of pairs of data to generate, number of rows in data set, or sample size
#' @param xmu x variable's mean
#' @param ymu y variable's mean
#' @param xsd x variable's standard deviation
#' @param ysd y varaible's standard deviation
#' @param r correlation between the x and y variable
#' @return data frame of uncensored data with n rows and 2 columns
#' @export
#' @examples 
#'
genData.N <- function( n, xmu, ymu, xsd, ysd, r){
  #initialize
  x <- vector("numeric", n)
  y <- vector( "numeric", n)
  #declare
  x <- rnorm(n)
  y <- sqrt((1-r^2))*rnorm(n)+r*x
  x <- x*xsd + xmu
  y <- y*ysd + ymu
  #return
  uncenData <- data.frame(x,y)
}