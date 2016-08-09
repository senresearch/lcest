
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
## theta = true values of bivariate normal (xmu,ymu,xsd,ysd,r)
## censorLevel = vector of length 2 indicating censoring level
## alpha = confidence levels will be 1-alpha level
################################################################

runTrial.N <- function(nsim,n,theta=c(0,0,1,1,0),
                         censorLevel=c(0,0),alpha=0.2)
{
  ## make matrix for returning value
  outputData <- matrix(nrow=nsim,ncol=5*5)
  
  ## loop through simulations
  for( i in 1:nsim)
  {
    ## generate data
    data1 <- genData.N(n, theta[1],theta[2],theta[3],
                     theta[4],theta[5])
    ## censor data per specs
    cenData1 <- censorData( data1, censorLevel)
    ## get optimResults
    result1 <- optimResults(cenData1,alpha=alpha)
    ## get bivariate summary
    exact <- bivariateSummary(data1)
    
    ## assign values to matrix
    outputData[i,1:5] <- exact
    outputData[i,6:10] <- result1$coefficients$estimate
    outputData[i,11:15] <- result1$coefficients$stdError
    outputData[i,16:20] <- result1$coefficients$lowerCI
    outputData[i,21:25] <- result1$coefficients$upperCI
  }
  
  ## names of the parameters
  parnames <- c("xmu","ymu","xsd","ysd","r")
  ## make names of all columns
  colnames(outputData) <- c(paste("exact",parnames,sep="."),
                            paste("estimate",parnames,sep="."),
                            paste("sd",parnames,sep="."),
                            paste("lowerCL",parnames,sep="."),
                            paste("upperCL",parnames,sep="."))
  ## return
  outputData
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

runTrial.T <- function(nsim,n,theta=c(0,0,1,1,0), df = 4,
                         censorLevel=c(0,0),alpha=0.2)
{
  ## make matrix for returning value
  outputData <- matrix(nrow=nsim,ncol=5*5)
  
  ## loop through simulations
  for( i in 1:nsim)
  {
    ## generate data
    locVec <- c(theta[1], theta[2])
    scaleMat <- buildScaleMat( theta[3], theta[4], theta[5], df)
    data1 <- genData.T( n, locVec, scaleMat, df )
    ## censor data per specs
    cenData1 <- censorData( data1, censorLevel)
    ## get optimResults
    result1 <- optimResults(cenData1,alpha=alpha)
    ## get bivariate summary
    exact <- bivariateSummary(data1)
    
    ## assign values to matrix
    outputData[i,1:5] <- exact
    outputData[i,6:10] <- result1$coefficients$estimate
    outputData[i,11:15] <- result1$coefficients$stdError
    outputData[i,16:20] <- result1$coefficients$lowerCI
    outputData[i,21:25] <- result1$coefficients$upperCI
  }
  
  ## names of the parameters
  parnames <- c("xmu","ymu","xsd","ysd","r")
  ## make names of all columns
  colnames(outputData) <- c(paste("exact",parnames,sep="."),
                            paste("estimate",parnames,sep="."),
                            paste("sd",parnames,sep="."),
                            paste("lowerCL",parnames,sep="."),
                            paste("upperCL",parnames,sep="."))
  ## return
  outputData
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
  b0 <- apply(outputData[,1:5],1,function(x) (x-theta))
  b1 <- apply(outputData[,6:10],1,function(x) (x-theta))
  m0 <- apply(b0,1,mean)
  m1 <- apply(b1,1,mean)
  m2 <- apply(b0^2,1,mean)
  m3 <- apply(b1^2,1,mean)
  out <- c(m0,m1,m2,m3)
  
  parnames <- c("xmu","ymu","xsd","ysd","r")
  ## make names of all columns
  names(out) <- c(paste("exact.bias",parnames,sep="."),
                  paste("estimate.bias",parnames,sep="."),
                  paste("exact.mse",parnames,sep="."),
                  paste("estimate.mse",parnames,sep="."))
  out
}

### usage:
### out <- runTrialNorm(400,700)
### checkCICoverage(out)
