#' Generate Bivariate Data
#' 
#' This function will generate an uncensored dataset distributed bivaraite normal or t.
#' @param n An integer, representing the number of pairs of data to generate, number of rows in data set, or sample size.
#' @param locVec A vector of length 2, representing the location parameter (equal to the mean vector when df>1).
#' @param scaleMat A symmetric positive-definite 2x2 matrix representing the scale matrix of the distribution, such that S*df/(df-2) is the variance-covariance matrix when df>2
#' @param df Integer greater than 3, representing degrees of freedom with df=Inf implying normal.
#' @return Returns a data set matrix with nrow = n and ncol= 2.
#' @details It is suggested to use buildScaleMat to generate the scaleMat parameter.
#' @importFrom mnormt rmt
#' @export
#' @examples 
#' xmu = 0
#' ymu = 0
#' xsd = 1
#' ysd = 1
#' r = 0
#' df = Inf #normal
#' scaleMat <- buildScaleMat( xsd, ysd, r, df)
#' genData(10, c(xmu, ymu), scaleMat, df)
#' 
genData<- function(n, locVec, scaleMat, df){
  if (df < 3 ) stop(" df must be greater than 3.")
  if ( df %% 1 != 0 && df != Inf) stop( "df must be integer.")
  uncenData <- rmt( n, locVec, scaleMat, df)
}

