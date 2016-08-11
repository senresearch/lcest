

################################################################
## function to return means, sds, and correlation of bivariate data
################################################################
bivariateSummary <- function(x)
{
  c(mean(x[,1]),mean(x[,2]),sd(x[,1]),sd(x[,2]),cor(x[,1],x[,2]))
}


################################################################
## function to run simulations and check confidence interval coverage
## nsim = number of simulations
## n = sample size
## theta = true values of bivariate t (xmu,ymu,xsd,ysd,r)
## df = degrees of freedom
## censorLevel = vector of length 2 indicating censoring level
## alpha = confidence levels will be 1-alpha level
################################################################


runTrial<- function(nsim,n,theta=c(0,0,1,1,0), df = Inf,
                         censorLevel=c(0,0),alpha=0.2, fileName ="trialOutput.csv")
{
  date <- date()
  ## make matrix for returning value
  outputData <- matrix(nrow=nsim,ncol=5*7)
  ## loop through simulations
  for( i in 1:nsim)
  {
    ## generate data
    locVec <- c(theta[1], theta[2])
    scaleMat <- buildScaleMat( theta[3], theta[4], theta[5], df)
    data1 <- rmt( n, locVec, scaleMat, df )
    ## censor data per specs
    cenData1 <- censorData( data1, censorLevel)
    cenData2 <- censorImpute( cenData1)
    ## get other censoring method parameters
    LODmethodPar <- defaultGuess(cenData1)
    LODSQRT2methodPar <- defaultGuess(cenData2)
    ## get optimResults
    result1 <- blcest(cenData1,alpha=alpha, df= df)
    ## get bivariate summary
    exact <- bivariateSummary(data1)
    
    ## assign values to matrix
    outputData[i,1:5] <- exact
    outputData[i,6:10] <- result1$coefficients$estimate
    outputData[i,11:15] <- result1$coefficients$stdError
    outputData[i,16:20] <- result1$coefficients$lowerCI
    outputData[i,21:25] <- result1$coefficients$upperCI
    outputData[i,26:30] <- LODmethodPar #LOD method
    outputData[i,31:35] <- LODSQRT2methodPar #LODsqrt2 method
  }
  
  ## names of the parameters
  parnames <- c("xmu","ymu","xsd","ysd","r")
  ## make names of all columns
  colnames(outputData) <- c(paste("exact",parnames,sep="."),
                            paste("estimate",parnames,sep="."),
                            paste("sd",parnames,sep="."),
                            paste("lowerCL",parnames,sep="."),
                            paste("upperCL",parnames,sep="."),
                            paste("LOD method", parnames, sep="."),
                            paste("LODSqrt2 method", parnames, sep="."))
  ## output to csv
  dataFile <- file(fileName, open="wt")
  on.exit(close(dataFile))
  writeLines(paste("# this csv file was created by trial_function.R on", date), con=dataFile)
  writeLines(paste("Parameters: nsim=", nsim, "; n=", n, "; theta= c(", theta[1], ";", theta[2], 
                  ";", theta[3], ";", theta[4], ";", theta[5], ")", "; df=", df,
                   "; censorLevel= c(", censorLevel[1], ";", censorLevel[2], ")", 
                   "; alpha=", alpha), con=dataFile)
  write.csv( outputData, dataFile)
}
#############################################################
## check coverage of confidence intervals
##
## outputData = output of runTrialNorm
## theta = true value of bivariate normal parameters
#############################################################

checkCICoverage <- function(outputData,theta=c(0,0,1,1,0))
{
  ## is lower limit less than true value
  l <- apply(outputData[,16:20],1,function(x) x<theta)
  ## is upper limit greater than true value
  u <- apply(outputData[,21:25],1,function(x) x>theta)
  ## is true value within confidence limits
  out <- apply((l)&(u),1,mean)
  ## names of the parameters
  names(out) <- c("xmu","ymu","xsd","ysd","r")
  ## return 
  out
}

#############################################################
## calculate bias of estimates
##
## outputData = output of runTrialNorm
## theta = true value of bivariate normal parameters
#############################################################

calcBiasMse <- function(outputData,theta=c(0,0,1,1,0))
{
  # for exact and estimate
  b0 <- apply(outputData[,1:5],1,function(x) (x-theta))
  b1 <- apply(outputData[,6:10],1,function(x) (x-theta))
  m0 <- apply(b0,1,mean)
  m1 <- apply(b1,1,mean)
  m2 <- apply(b0^2,1,mean)
  m3 <- apply(b1^2,1,mean)
  
  # for other methods
  b2 <- apply(outputData[,26:30],1,function(x) (x-theta))
  b3 <- apply(outputData[,31:35],1,function(x) (x-theta))
  m4 <- apply(b2,1,mean)
  m5 <- apply(b3,1,mean)
  m6 <- apply(b2^2,1,mean)
  m7 <- apply(b3^2,1,mean)
  
  out <- c(m0,m1,m2,m3, m4, m5, m6, m7)
              
  parnames <- c("xmu","ymu","xsd","ysd","r")
  ## make names of all columns
  names(out) <- c(paste("exact.bias",parnames,sep="."),
                  paste("estimate.bias",parnames,sep="."),
                  paste("exact.mse",parnames,sep="."),
                  paste("estimate.mse",parnames,sep="."),
                  paste("LOD.bias",parnames,sep="."),
                  paste("LODsqrt2.bias",parnames,sep="."),
                  paste("LOD.mse",parnames,sep="."),
                  paste("LODsqrt2.mse",parnames,sep="."))
  out
}

### usage:
### out <- runTrialNorm(400,700)
### checkCICoverage(out)
