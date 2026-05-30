## this code provides an illustration of the Metropolis-Hastings 
## algorithm

set.seed(0)
## DEFINE A TARGET DENSITY (GIVEN BY THE PROBLEM)
target_density<-function(x,mu,lambda){
  return(exp(-t(x-mu)%*%lambda%*%(x-mu)/2))
}

## CHOOSE A PROPOSAL DISTRIBUTION (CHOSEN BY THE USER)
## IN THIS CASE THE PROPOSAL IS A N(x,proposal_sd^2) distribution
proposal_sd<-0.5#proposal standard deviation of reasonable size
# IF YOU CHANGE THE VALUE OF proposal_sd THE EFFICIENCY OF
# THE ALGORITHM CHANGES
# proposal_sd<-0.05#proposal standard deviation too small
# proposal_sd<-3#proposal standard deviation too big

## PLOT THE TARGET DENSITY (EASY TO DO BECAUSE WE ARE IN 2D)
grid <- seq(-3,3,length.out=100)
z <- matrix(0,nrow=100,ncol=100)
mu <- c(0,0)
sigma <- matrix(c(1,0.5,0.5,0.5),nrow=2)
lambda<-solve(sigma)
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- target_density(x=c(grid[i],grid[j]),mu = mu,lambda = lambda)
  }
}
par(mfrow=c(1,1))
image(grid,grid,z,xlim=c(-3,3),ylim = c(-3,3))

## RUN THE METROPOLIS HASTINGS ALGORITHM STEP-BY-STEP
#readline(prompt="Press a key to propose a value ")
current_state<-c(-2,-1)
points(current_state[1],current_state[2],pch=1)
points(current_state[1],current_state[2],pch=20)

proposed_state<-current_state+rnorm(2,mean = 0,sd=proposal_sd)
points(proposed_state[1],proposed_state[2])
alpha<-min(1,target_density(proposed_state,mu,lambda)/
               target_density(current_state,mu,lambda))
print(paste("alpha=",alpha))
if(runif(1)<alpha){
  current_state<-proposed_state
  points(current_state[1],current_state[2],pch=20)
}else{
  current_state<-current_state
  points(proposed_state[1],proposed_state[2],col="red")
}


## RUN THE METROPOLIS HASTINGS ALGORITHM FOR s STEPS
## AND PLOT THE RESULTING SAMPLES
s<-10^4
x<-rep(NA,s)
y<-rep(NA,s)
x[1]<-current_state[1]
y[1]<-current_state[2]
for(t in 1:s){
  proposed_state<-c(x[t],y[t])+rnorm(2,mean = 0,sd=proposal_sd)
  alpha<-min(1,target_density(proposed_state,mu,lambda)/
               target_density(c(x[t],y[t]),mu,lambda))
  if(runif(1)<alpha){
    x[t+1]<-proposed_state[1]
    y[t+1]<-proposed_state[2]
  }else{
    x[t+1]<-x[t]
    y[t+1]<-y[t]
  }
}
points(x,y,pch=20,cex=0.1)
title(paste("Samples after",s,"steps"))




## THE METROPOLIS-HASTINGS ALGORITHM CAN BE
## APPLIED TO ANY TARGET DENSITY (PROVIDED WE CAN 
## EVALUATE IT POINT-WISE UP TO NORMALIZING CONSTANT)
target_density<-function(x,B){
  return(
    exp((
      -x[1]^2/200-0.5*(x[2]+B*x[1]^2-100*B)^2
      +0.5*(100*B)^2
    ))
  )
}
proposal_sd<-1.5
grid <- seq(-15,15,length.out=100)
z <- matrix(0,nrow=100,ncol=100)
B<-0.05
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- target_density(x=c(grid[i],grid[j]),B=B)
  }
}
image(grid,grid,z, xlim = c(-15,15),ylim = c(-15,15))
current_state<-c(-2,-1)
points(current_state[1],current_state[2],pch=1)
points(current_state[1],current_state[2],pch=20)

## RUN THE METROPOLIS HASTINGS ALGORITHM FOR 10^2 STEPS
s<-10^2
x<-rep(NA,s)
y<-rep(NA,s)
x[1]<-current_state[1]
y[1]<-current_state[2]
for(t in 1:s){
  proposed_state<-c(x[t],y[t])+rnorm(2,mean = 0,sd=proposal_sd)
  alpha<-min(1,target_density(proposed_state,B=B)/
               target_density(c(x[t],y[t]),B=B))
  if(runif(1)<alpha){
    x[t+1]<-proposed_state[1]
    y[t+1]<-proposed_state[2]
  }else{
    x[t+1]<-x[t]
    y[t+1]<-y[t]
  }
}
points(x,y,pch=20,cex=0.1)
title(paste("Samples after",s,"steps"))
## RUN THE METROPOLIS HASTINGS ALGORITHM FOR 10^4 STEPS
last_state<-c(x[t],y[t])
s<-10^4
x<-rep(NA,s)
y<-rep(NA,s)
x[1]<-last_state[1]
y[1]<-last_state[2]
for(t in 1:s){
  proposed_state<-c(x[t],y[t])+rnorm(2,mean = 0,sd=proposal_sd)
  alpha<-min(1,target_density(proposed_state,B=B)/
               target_density(c(x[t],y[t]),B=B))
  if(runif(1)<alpha){
    x[t+1]<-proposed_state[1]
    y[t+1]<-proposed_state[2]
  }else{
    x[t+1]<-x[t]
    y[t+1]<-y[t]
  }
}
image(grid,grid,z,xlim = c(-15,15),ylim = c(-15,15))
points(x,y,pch=20,cex=0.1)
title(paste("Samples after",s,"steps"))
