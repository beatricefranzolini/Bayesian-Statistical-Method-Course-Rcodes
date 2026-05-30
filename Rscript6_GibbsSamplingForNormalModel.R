## IN THIS CODE WE PERFORM GIBBS SAMPLING ON A 
## UNIVARIATE NORMAL MODEL WITH UNKNOWN MEAN AND VARIANCE

##  DEFINE THE PRIOR HYPERPARAMETERS
mu0<-0
sigma0sq<-100
nu0<-1
xi0sq<-0.01

##  SIMULATE A VECTOR OF OBSERVED DATAPOINTS
n<-10
true_mean<-3; true_var<-2
x<-rnorm(n,mean = true_mean,sd = sqrt(true_var))
##  COMPUTE RELEVANT SUMMARY STATISTICS
mean.x<-mean(x)

## DEFINE EMPTY MATRIX OF POSTERIOR SAMPLES
s<-10000
theta.samples<-matrix(NA,nrow = s,ncol = 2)

## INITIALIZE THE MARKOV CHAIN (mu,sigmasq) AT SOME ARBITRARY VALUES
mu<-1
sigmasq<-1
## SET TIME TO 1 AND STORE THE FIRST SAMPLE
t<-1
theta.samples[t,]<-c(mu,sigmasq)

for (t in 2:s){# SIMULATE THE MARKOV CHAIN (mu,sigmasq) FROM t=1 to t=s
  
  # COMPUTE THE UPDATED PARAMETERS OF THE FULL CONDITIONAL OF mu
  mun<-(mu0/sigma0sq+n*mean.x/sigmasq)/(1/sigma0sq+n/sigmasq)
  sigmansq<-1/(1/sigma0sq+n/sigmasq)
  # SIMULATE A NEW mu FROM ITS FULL CONDITIONAL
  mu<-rnorm(1,mean = mun,sd = sqrt(sigmansq))
  
  # COMPUTE THE UPDATED PARAMETERS OF THE FULL CONDITIONAL OF sigmasq
  nun<-nu0+n
  xinsq<-(nu0*xi0sq+sum((x-mun)^2))/nun
  # SIMULATE A NEW sigmasq FROM ITS FULL CONDITIONAL
  sigmasq<-1/rgamma(n = 1,shape = nun/2,rate = nun*xinsq/2)
  
  # STORE THE NEW VALUE OF THE MARKOV CHAIN IN THE MATRIX OF POSTERIOR SAMPLES
  theta.samples[t,]<-c(mu,sigmasq)
}

### PLOT THE TRAJECTORIES OF (mu,sigmasq) AND THE CORRESPONDING 
### samples (removing the first few iterations as "burn in")
###
par(mfrow=c(2,1))
plot(theta.samples[,1],type="l",main="trajectory of mu^(t)")
plot(theta.samples[,2],type="l",main="trajectory of sigmasq^(t)")

plot(theta.samples[-c(1:10),1],type="p",main="posterior samples from mu^(t)")
plot(theta.samples[-c(1:10),2],type="p",main="posterior samples from sigmasq^(t)")
par(mfrow=c(1,1))

### NOW YOU CAN DO ANY INFERENCE OF INTEREST USING THE SAMPLES
### FOR EXAMPLE: A GOOD SANITY CHECK WOULD BE TO COMPARE THE 
### POSTERIOR DISTRIBUTION WITH THE "TRUE" VALUE WE USED TO GENERATE THE SAMPLES
###
par(mfrow=c(2,1))
hist(theta.samples[-c(1:10),1],
     xlab="values of mu",ylab="Estimated posterior density",
     main="Posterior histogram for mu")
points(x = true_mean,y = 0,col="red",pch=16)
hist(sqrt(theta.samples[-c(1:10),2]),
     xlab="values of sigma",
     ylab="Estimated posterior density",
     main="Posterior histogram for sigma")
points(x = sqrt(true_var),y = 0,col="red",pch=16)
par(mfrow=c(1,1))
## THIS IS USING KERNEL DENSITY ESTIMATORS TO GET SMOOTH CURVES
par(mfrow=c(2,1))
plot(density(theta.samples[-c(1:10),1]),xlab="values of mu",
     ylab="Estimated posterior density",
     main="Approximate posterior density for mu")
points(x = true_mean,y = 0,col="red",pch=16)

plot(density(sqrt(theta.samples[-c(1:10),2])),
     xlab="values of sigma",
     ylab="Estimated posterior density",
     main="Approximate posterior density for sigma")
points(x = sqrt(true_var),y = 0,col="red",pch=16)
par(mfrow=c(1,1))

### THIS IS HOW THE JOINT POSTERIOR DISTRIBUTION LOOKS LIKE (INCREASE s TO GET BETTER RESULTS)
###
plot(theta.samples[-c(1:10),],
     xlab="mu",
     ylab="sigma squared",
     main="joint posterior distribution of (mu,sigma^2)")
