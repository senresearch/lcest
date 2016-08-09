# Purpose: Converts a Var Cov Parameter (matrix, vec, or single parameter) into a scale parameter
# Parameter: varCovarPara = anything to multiply by df based constant to turn to scale parameter
# Returns: new parameter now in scale with df
# Note: based on the fact S*df/(df-2) = var cov matrix (see ?dmt)
# Note: does NOT convert theta
convertToScale <- function( varCovPara, df){
  if( df != Inf) scalePara <- varCovPara*(df -2)/df
  if( df == Inf) scalePara <- varCovPara
}

# Purpose: Converts a scale Parameter (matrix, vec, or single parameter) into a Var Cov parameter
# Parameter: scalePara = anything to multiply by df based constant to turn to varCov
# Returns: new parameter now as varCov
# Note: based on the fact S*df/(df-2) = var cov matrix (see ?dmt)
# Note: does NOT convert theta
convertFromScale <- function( scalePara, df){
  if( df != Inf) varCovPara <- scalePara*df/(df -2)
  if( df == Inf) varCovPara <- scalePara
}

# Purpose: Converts theta from c(xmu, ymu, xsd, ysd, r) to c( xmu, ymu, sigma11, sigma22, sigma21)
#     where sigma11 = xvar*(df -2)/ df , sigma22 = yvar*(df -2)/ df  , and sigma21 = r*ysd*xsd*(df-2)/df
# Parameters: theta with sds and corr
# Retruns: theta with sigma11, sigma21, and sigma22
convertThetaToScale<- function( oldTheta, df){
  #convert 
  sigma11 <- convertToScale(oldTheta[3]^2, df)
  sigma22 <- convertToScale(oldTheta[4]^2, df)
  sigma21 <- convertToScale( oldTheta[3]*oldTheta[4]*oldTheta[5], df)
  #Return
  c( oldTheta[1], oldTheta[2], sigma11, sigma22, sigma21)
}


# Purpose: Converts theta from c(xmu, ymu, sigma11, sigma22, sigma21) to c(xmu, ymu, xsd, ysd, r)
#     where sigma11 = xvar*(df -2)/ df , sigma22 = yvar*(df -2)/ df  , and sigma21 = r*ysd*xsd*(df-2)/df
# Parameters: theta with sigma11, sigma21, and sigma22
# Retruns: theta with sds and corr
convertThetaFromScale <- function( scaleTheta, df){
  #convert
  xsd <- sqrt(convertFromScale( scaleTheta[3], df))
  ysd <- sqrt(convertFromScale( scaleTheta[4], df))
  r <- scaleTheta[5]/sqrt(scaleTheta[3])/sqrt(scaleTheta[4])
  #return
  c(scaleTheta[1], scaleTheta[2], xsd, ysd, r)
}
