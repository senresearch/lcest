require(lcest)
set.seed(6724)
nsim = 500
n = 700

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(0,0),alpha=0.2, fileName ="trialOutput1T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(0,.1),alpha=0.2, fileName ="trialOutput2T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(.1,.1),alpha=0.2, fileName ="trialOutput3T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(.2,0),alpha=0.2, fileName ="trialOutput4T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(.2,.2),alpha=0.2, fileName ="trialOutput5T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(0,.4),alpha=0.2, fileName ="trialOutput6T.csv")

runTrial(nsim,n,theta=c(0,0,1,1,0), df = 4,
         censorLevel=c(.4,.4),alpha=0.2, fileName ="trialOutput7T.csv")

set.seed(6724)
runTrial(500,700,theta=c(0,0,1,1,0), df = Inf,
         censorLevel=c(.1,.1),alpha=0.2, fileName ="trialOutputMaxit.csv")

