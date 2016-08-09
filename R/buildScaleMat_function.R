#' Build Scale Matrix
#' 
#' This function builds a 2x2 scale matrix given xsd, ysd, corr, and df.
#' @param xsd standard deviation of x.
#' @param ysd standard deviation of y.
#' @param r correlation of x and y.
#' @param df Degrees of Freedom must be greater than 3 and an integer or Inf for normal.
#' @return 2x2 matrix
#' @export
#' @examples 
#' buildScaleMat( 1,1,0,Inf)
#'
# Purpose:  Builds the Var cov matrix from corr and sds requested and converts it using convertToScale
# Parameters: requested corr and sds to change to scale matrix
# Returns: a 2 x 2 scale matrix
buildScaleMat <- function(sdx, sdy, r, df){
  # get variance covaraince matrix
  varcovarMat <-matrix( c(sdx^2, sdx*sdy*r, sdx*sdy*r, sdy^2), nrow =2, ncol =2)
  # convert to scale matrix by calling convert to scale which multiplies by df based coeffeicent 
  scaleMat <- convertToScale( varcovarMat, df)
}