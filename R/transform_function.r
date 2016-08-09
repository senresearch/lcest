
#' Transform Theta
#' 
#' 
#' This function transforms theta to impose proper/logical bounds for std dev and rho
#' @param thetaI vector of length 5, such that c(mux, muy, sigx, sigy, r). I stands for inverted/orignal.
#' @return thetaT vector of length 5, such that c( mux, muy, log(sigx), log(sigy), atanh(r))
#' 
transformTheta <-function(thetaI){
  # assign
  thetaT <- thetaI
  # tarnsfrom sd
  thetaT[3:4] <- log(thetaI[3:4])
  # transform rho
  thetaT[5] <- atanh(thetaI[5])
  # return
  thetaT
}


#' Transform Theta Inverse
#' 
#' 
#' This function undoes the transformation to theta from the transfromTheta function.
#' @return thetaI vector of length 5, such that c(mux, muy, sigx, sigy, r). I stands for inverted/orignal.
#' @param thetaT vector of length 5, such that c( mux, muy, log(sigx), log(sigy), atanh(r))
#' 
transformThetaInv <-function(thetaT){
  # assign
  thetaI <- thetaT
  # untransfrom sd
  thetaI[3:4] <- exp(thetaT[3:4])
  # untransfrom rho
  thetaI[5] <- tanh(thetaT[5])
  # return
  thetaI
}

