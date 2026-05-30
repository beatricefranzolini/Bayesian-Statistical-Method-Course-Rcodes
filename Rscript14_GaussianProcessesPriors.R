#Gaussian processes simulation of paths

library(MASS) #to simulate from multivariate normals

#define an informative mean function
m <- function(x){
  10 * x * sin(0.6 * x)
}

x = seq(0, 10, 0.1)  #create a grid
y = m(x)
plot(x, y, type = 'l', lwd = 2, cex.axis = 1, cex.lab = 1)

#define a non informative mean function
m <- function(x){
  0*x
}

x = seq(0, 10, 0.1)  #create a grid
y = m(x)
plot(x, y, type = 'l', lwd = 2, cex.axis = 1, cex.lab = 1)

#Exponential square kernel ###########################################
K <- function(x, xprime, sigma2 = 2, ell = 10) {
  
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
  return(sigma2 * exp(-(x - xprime)^2 / ell))
}

#plot a the Gaussian process
set.seed(0)
x = seq(0, 10, 0.1)

# Compute Mean
mu = m(x)

# Compute Covariance
Sigma = K(x, x, sigma2 = 2, ell = 10)

nsim = 100
out = array(NA, dim = c(nsim, length(x)))
for (i in 1:nsim){
  out[i, ] = mvrnorm(1, mu, Sigma)
}

par(mfrow = c(1, 2))
# Plot 1: Plotting all the simulated paths
matplot(x, out[1,], type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')
matplot(x, t(out), type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')


# Compute Covariance
Sigma = K(x, x, sigma2 = 0.1, ell = 0.1)

nsim = 100
out = array(NA, dim = c(nsim, length(x)))
for (i in 1:nsim){
  out[i, ] = mvrnorm(1, mu, Sigma)
}
# Plot 2: Plotting all the simulated paths
matplot(x, out[1,], type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')
matplot(x, t(out), type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')

# Compute Covariance
Sigma = K(x, x, sigma2 = 0.1, ell = 1)

nsim = 100
out = array(NA, dim = c(nsim, length(x)))
for (i in 1:nsim){
  out[i, ] = mvrnorm(1, mu, Sigma)
}
# Plot 2: Plotting all the simulated paths
matplot(x, out[1,], type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')
matplot(x, t(out), type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')

#White noise kernel ###########################################
K_WN <- function(x, xprime, sigma2 = 2) {
  
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
  return(sigma2 * (x == xprime))
}

#plot a the Gaussian process
set.seed(0)
x = seq(0, 10, 0.1)

# Compute Mean
mu = m(x)

# Compute Covariance
Sigma = K_WN(x, x, sigma2 = 2)

nsim = 100
out = array(NA, dim = c(nsim, length(x)))
for (i in 1:nsim){
  out[i, ] = mvrnorm(1, mu, Sigma)
}

par(mfrow = c(1, 1))
matplot(x, out[1,], type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')

par(mfrow = c(1, 2))
matplot(x, out[2,], type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')

matplot(x, t(out), type = 'l', col = 'dodgerblue', main = 'Simulated Paths')
lines(x, mu, type = 'l', lwd = 3, col = 'black')
