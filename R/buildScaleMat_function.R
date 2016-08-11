#' Build Scale Matrix
#' 
#' This function builds a 2x2 scale matrix given xsd, ysd, corr, and df.
#' @param xsd standard deviation of x. Must be positive.
#' @param ysd standard deviation of y. Must be positive.
#' @param r correlation of x and y. Cannot be less than -1 or more than 1.
#' @param df Degrees of Freedom must be greater than 3 and an integer or Inf for normal.
#' @return 2x2 matrix
#' @export
#' @examples 
#' buildScaleMat( 1,1,0,Inf) #normal 
#' buildScaleMat( 1,1,0,4) #t with df =4
#'

buildScaleMat <- function(xsd, ysd, r, df){
  # force bounds
  if( xsd <= 0 || ysd <= 0){ stop("xSd and ySd must be positive.")}
  if (r > 1 || r < -1) {stop("r must be between or equal to -1 and 1.")}
  if (df < 3 ) {stop(" df must be greater than 3.")}
  if ( df %% 1 != 0 && df != Inf) {stop( "df must be integer.")}
  # get variance covaraince matrix
  varcovarMat <-matrix( c(xsd^2, xsd*ysd*r, xsd*ysd*r, ysd^2), nrow =2, ncol =2)
  # convert to scale matrix by calling convert to scale which multiplies by df based coeffeicent 
  scaleMat <- convertToScale( varcovarMat, df)
  scaleMat
}
