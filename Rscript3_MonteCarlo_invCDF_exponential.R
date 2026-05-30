rm(list=ls())
set.seed(0) #set seed for replicability

# Monte Carlo samples from an Exponential distribution
# set some value for lambda
lambda<-3
#sample pseudo-random number uniform
U_vec <- runif(10000)
#Cdf inverse transformation
X_vec <- -log(1-U_vec)/lambda

hist(X_vec, prob = TRUE)

#plot the density of an exponential random variable to check the algorithm 
grid = seq(0,3,0.1)
dens <- function(grid){lambda*exp(-lambda*grid)}

lines(grid, dens(grid), col = "red")

