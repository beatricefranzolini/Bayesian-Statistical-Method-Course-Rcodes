## SIMULATE A DATASET
set.seed(0)
n<-30
y<-c(rpois(n/2,lambda = 10),rpois(n/2,lambda = 15))

## SET PRIOR HYPERPARAMETERS
alpha1<-1;beta1<-0.1
alpha2<-1;beta2<-0.1

# BUILD AN EMPTY MATRIX OF POSTERIOR SAMPLES
s<-10000
posterior.samples<-matrix(NA,nrow = s,ncol = 3)
# initialize unknown parameters (M,lambda1,lambda2)
M<-sample(x = c(1:(n-1)),size = 1)
lambda1<-mean(y[c(1:M)])
lambda2<-mean(y[-c(1:M)])
t<-1
posterior.samples[t,]<-c(M,lambda1,lambda2)

for (t in 2:s){# SIMULATE THE MARKOV CHAIN (M,lambda1,lambda2) FROM t=1 to t=s
  
  # COMPUTE FULL CONDITIONAL OF M
  log_weights<-cumsum(y[-n])*log(lambda1)+(sum(y)-cumsum(y[-n]))*log(lambda2)+(lambda2-lambda1)*c(1:(n-1))
  # SIMULATE A NEW M FROM ITS FULL CONDITIONAL
  M<-sample.int(n = n-1,size = 1,prob = exp(log_weights-max(log_weights)))  
  
  # SIMULATE A NEW lambda1 and lambda2 FROM ITS FULL CONDITIONAL
  lambda1<-rgamma(1,shape = alpha1+sum(y[c(1:M)]),rate = beta1+M) 
  lambda2<-rgamma(1,shape = alpha2+sum(y[-c(1:M)]),rate = beta2+n-M) 
  
  # STORE THE NEW VALUE OF THE MARKOV CHAIN IN THE MATRIX OF POSTERIOR SAMPLES
  posterior.samples[t,]<-c(M,lambda1,lambda2)
}

## PLOT THE TRAJECTORY OF THE MARKOV CHAIN
par(mfrow=c(3,1))
plot(posterior.samples[,1],type="l")
plot(posterior.samples[,2],type="l")
plot(posterior.samples[,3],type="l")
par(mfrow=c(1,1))

## PLOT THE DATASET AND THE POSTERIOR DISTRIBUTION OF M DEPENDING ON THE 
## LEVEL OF NOISE THE MODEL WILL BE ABLE TO RECOVER THE TRUE M OR 
## NOT (HERE THERE IS A TRUE M BECAUSE WE SIMULATED THE DATA FROM THE MODEL)
par(mfrow=c(2,1))
plot(y)
hist(posterior.samples[,1],breaks = c(0:(n-1))+0.5)
points(x = n/2,y = 0,col="red",pch=16)
par(mfrow=c(1,1))
## ONE INTERESTING ASPECT IS THAT THE BAYESIAN PROCEDURE PROVIDES
## UNCERTAINTY QUANTIFICATION FOR M
