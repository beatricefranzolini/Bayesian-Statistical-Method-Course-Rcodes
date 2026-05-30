rm(list=ls())
set.seed(0) #set seed for replicability

### Set number of samples
s<-1000

### Simulate points
U_matrix<-matrix(runif(2*s,min = -1,max = 1),nrow = s,ncol = 2)

## Plot of simulated points and circle
plot(U_matrix,xlim=c(-1,1),ylim=c(-1,1),asp=1)
angle<-seq(0,2*pi,0.01)
lines(cos(angle),sin(angle),col="red",lwd=3)

### Compute Y sample
Y_samples<-4*(rowSums(U_matrix^2)<1)

## Apply law of large numbers to get the approximation
print(paste("approximation of pi=",mean(Y_samples)))
