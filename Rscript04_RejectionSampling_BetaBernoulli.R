rm(list=ls())
set.seed(0) #set seed for replicability

## LET'S ASSUME WE HAVE OBSERVED x_obs=1
x_obs<-1

## Let's use rejection sampling to simulate from the posterior

# STEP (1): SIMULATE theta FROM THE PRIOR AND X FROM THE LIKELIHOOD
theta<-runif(1)
x<-rbinom(1,size = 1,prob=theta)

# STEP (2): ACCEPT theta IF X=x_obs, OTHERWISE REPEAT STEP (1)
while(x!=x_obs){
  theta<-runif(1)
  x<-rbinom(1,size = 1,prob=theta)
}

## theta is a sample from the posterior p(theta|X=x_obs)
theta

## repeat the procedure s times
s<-10000
RS_output<-rep(NA,s)
for (i in 1:s){
  theta<-runif(1)
  x<-rbinom(1,size = 1,prob=theta)
  while(x!=x_obs){
    theta<-runif(1)
    x<-rbinom(1,size = 1,prob=theta)
  }
  RS_output[i]<-theta
}

## theta was initially simulated from a uniform distribution
## but the accept/reject mechanism "transforms" them into Beta (1,2)
## comparison with theta without accept/reject
plot(runif(s),main="Samples before accept/reject",xlab="Sample number",ylab="samples from the prior")
plot(RS_output,main="Samples after l'accept/reject",xlab="Sample number",ylab="samples from the posterior")

## plot the posterior
hist(RS_output,freq = F,main="Histogram of posterior samples")

