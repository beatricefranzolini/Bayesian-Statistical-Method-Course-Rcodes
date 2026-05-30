#Gaussian processes regression 
rm(list=ls())

library(MASS) #to simulate from multivariate normals

#set the prior
#define a non informative mean function
m <- function(x){
  0*x
}
#define a Exponential square kernel with periodic component
K <- function(x, xprime, sigma2 = 2, ell = 10, 
              sigma2_2 = 0.1, ell2 = 1, p = 1) {
  
  # Force it to column matrix
  x = matrix(x)
  xprime = matrix(xprime)
  p= dim(x)[1];    
  n= dim(xprime)[1];
  
  # Broadcast the matrix
  x = matrix(x, nrow = p, ncol = n)
  xprime = matrix(xprime, nrow = n, ncol = p)
  
  # Compute the kernel
  xprime = t(xprime)
  return(sigma2 * exp(-(x - xprime)^2 / ell) +
          + sigma2_2 * exp(- 2/ell2 * sin(pi /p *abs(x - xprime)) ) )
}

#Simulate data 
true_link <- function(x){
  1.5*x + 5 * sin(x)
}
x_obs = rnorm(500, 0, 4)
y_obs = true_link(x_obs) + rnorm(100, 0, 5)
par(mfrow = c(1, 1))
plot(x_obs,y_obs, ylim =c(-30,30))
grid = seq(min(x_obs),max(x_obs),0.2)
lines(grid, true_link(grid), type = 'l', lwd = 3, col = 'red')

#estimate the model 
m_post <- function(x, y_obs, x_obs, sigmay = 5){
  K(x,x_obs)%*%solve(K(x_obs,x_obs) + sigmay*diag(length(x_obs)))%*%(y_obs - m(y_obs)) 
}

K_post  <- function(x, xprime, y_obs, x_obs, sigmay = 5){
  return(K(x,xprime) - 
           K(x,x_obs)%*%solve(K(x_obs,x_obs) + sigmay*diag(length(x_obs)))%*%
           K(x_obs,xprime) )
}


#plot the posterior
set.seed(0)
x = seq(min(x_obs)-2, max(x_obs)+2, 0.5)

# Compute Mean
mu = m_post(x, y_obs, x_obs)

# Compute Covariance
Sigma = K_post(x, x, y_obs, x_obs)

nsim = 10000
out = array(NA, dim = c(nsim, length(x)))
for (i in 1:nsim){
  out[i, ] = mvrnorm(1, mu, Sigma)
}

par(mfrow = c(1, 1))
matplot(x, t(out), type = 'l', col = 'dodgerblue', main = 'Posterior', ylim =c(-30,30))
lines(x, mu, type = 'l', lwd = 3, col = 'black')
lines(x, true_link(x), type = 'l', lwd = 3, col = 'red')
points(x_obs,y_obs)
