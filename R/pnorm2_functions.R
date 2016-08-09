#####################################################################
# Function to calculate bivariate normal distribution function
#
# Based on Genz(2004); tries to take advantage of the fact that we
# want the log of the probability.
#####################################################################

##################################################
# bivariate normal distribution function
# vectorized for first argument
#
# x = matrix with two columns
# m = mean (vector of length 2)
# s = standard deviation (vector of length 2)
# r = correlation
# log.p = if log probability is wanted
##################################################


pnorm2 <- function(x,m,s,r,log.p=TRUE)
{
  apply(x,1,pnorm2.uv,m=m,s=s,r=r,log.p=log.p)
}


##################################################
# bivariate standard normal distribution function
# 
# h,k = standard normal values
# rho = correlation
# log.p = if log probability is wanted
##################################################

pnorm2.std <- function(h,k,rho,log.p=TRUE)
{
  if(rho==0)
    ans <- pnorm(h) * pnorm(k)
  else 
    ans <- pnorm(h) * pnorm(k) + F(h,k,rho)
  
  if( log.p )
    ans <- log(ans)
  
  ans
}

##################################################
# bivariate normal distribution function
# (not vectorized)
#
# x = vector of length 2
# m = mean (vector of length 2)
# s = standard deviation (vector of length 2)
# r = correlation
# log.p = if log probability is wanted
##################################################

pnorm2.uv <- function(x,m,s,r,log.p=TRUE)
{
  h <- (x[1]-m[1])/s[1]
  k <- (x[2]-m[2])/s[2]
  
  pnorm2.std(h,k,r,log.p=log.p)
}

################################
# integrand in correction term
#
# r = correlation
# h,k = standard normal limits
################################

f <- function(r,h,k)
{
  fac <- (1-r^2)
  exp(-(h^2+k^2-2*r*h*k)/(2*fac))/sqrt(fac)/(2*pi)
}

################################
# correction term
#
# r = correlation
# h,k = standard normal limits
################################

F <- function(h,k,rho)
{
  integrate(f,0,rho,h=h,k=k)$value
}
