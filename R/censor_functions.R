#' Left Censor Data
#' 
#' This function will take a data set and artificially left censor it to a request level of censoring.
#' @param uncenData Data matrix of any size.
#' @param cenLevelVec vector of length equal to rows of data, indicating censoring level requested in decimals.
#' @details The function sets the limit of detection at the value which the requested amount of data is below.
#'          All values below this limit of detection are replaced with the limit of detection.
#' @return Matrix such that the first half of the columns are the censored data and the 
#'         second haf are the isCensored flags with 1 indicating censoring, where censored data has a value
#'         equal to the limit of detection.
#' @export
#' @examples 
#' xmu = 0
#' ymu = 0
#' xsd = 1
#' ysd = 1
#' r = 0
#' df = Inf #normal
#' scaleMat <- buildScaleMat( xsd, ysd, r, df)
#' myData <- genData(10, c(xmu, ymu), scaleMat, df)
#' censorData(uncenData = myData, cenLevelVec =c(.2,.2))
#'
censorData <- function( uncenData, cenLevelVec){
  # find num of col and rows in data matrix
  numCol = length(uncenData[1,])
  numRow = length(uncenData[,1])
  #Declare flag matrix and fill with 0s
  flagMat <- matrix( rep(0, times = numCol*numRow), nrow = numRow, ncol = numCol)
  #loop thorugh data matrix to fill in flag matrix and replace "censored" values by LOD
  for ( c in 1:numCol){
    #Set limit of detection to catch level censoring requested
    LoD <- quantile(uncenData[,c], cenLevelVec[c])
    #Fill by row
    for ( r in 1:length(uncenData[,c])){
      if (uncenData[r,c] < LoD){
        flagMat[r,c] <- 1
        uncenData[r,c] <- LoD # now the uncenData is censored!!!
      }
    } 
  }
  # now the uncenData is censored!!!
  #bind and return altered data matrix and flag together
  cenData <- cbind(uncenData,flagMat)
}

#' Censor Impute
#' 
#' This function will take already censored data and replace LOD with LOD/sqrt2 or another value
#' @param cenData Data matrix such that the first half of the columns are the censored data and the 
#'         second haf are the isCensored flags with 1 indicating censoring, where censored data has a value
#'         equal to the limit of detection.
#' @param replacement The value request to fill censored values in the censored data matrix. 
#'        By Defualt equals LOD divided by sqrt2 as this is a common method.
#' @details To change the censored values to LOD/sqrt2 the data passed to the function must be censored
#'          such that censored values = LOD. If LOD is not availible, replacement should only equal numerical values. 
#' @return Matrix such that the first half of the columns are the censored data and the 
#'         second haf are the isCensored flags with 1 indicating censoring, where censored data has a value
#'         equal to replace parameter.
#' @export
#' @examples 
#'
censorImpute <- function( cenData, replacement =LOD/sqrt(2) ){
  # find num of col and rows in data matrix
  numCol = length(cenData[1,])
  numRow = length(cenData[,1])
  # loop through and change value if flag = 1
  for( dataCol in 1:(numCol/2)){
    flagCol = numCol/2 + dataCol
    for ( r in 1:numRow){
      if( cenData[r, flagCol] == 1){
        LOD <- cenData[r,dataCol]
        cenData[r,dataCol] <- replacement
      } else cenData[r,dataCol] <- cenData[r,dataCol]
    }
  }
  cenData
}


