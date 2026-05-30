library(MCMCprecision) #used to sample from a Dirichlet distribution via rdirichlet
library(salso) #to derive a point estimate of the partition
library(fossil) #to compute rand index
#Finite mixture model with Normal kernel and NormalInverseGaussian prior

set.seed(0)
#simulate data
n = 500
H_true = 2 
mean_true = c(-2, 2)
sigma2_true = c(1, 1)
lambda_true = c(0.4, 0.6)
z_true = sample(c(1,2), n, replace = TRUE, prob = lambda_true)

y = rep(NA,n)
for (h in 1: H_true){
  y[z_true==h]=rnorm(sum(z_true==h), mean_true[h], sqrt(sigma2_true[h]))
}
y = (y - mean(y))/sd(y)
hist(y, breaks = 30)

#prior 
H = 5
alpha = rep(1/H, H)
mu0 = 0 
k0 = 1 
v0 = 4 
sigma20 = 2

#number of iterations 
S = 10000


#objects where to save the posterior samples 
z = matrix(NA, nrow = S+1, ncol = n)
lambda = matrix(NA, nrow = S, ncol = H)
mu = matrix(NA, nrow = S, ncol = H)
sigma2 = matrix(NA, nrow = S, ncol = H)

#Gibbs Sampler 
#initialize 
z[1,] = sample(H, n, replace = TRUE) #random initialization

#the following lines create a progress bar 
pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = S, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


for (s in 1:S){
  nh_vector = as.vector(table(factor(z[s,], levels = 1:H)))
  lambda[s, ] = rdirichlet(1, alpha + nh_vector)
  for (h in 1:H){
    vn = v0 + nh_vector[h]
    kn = k0 + nh_vector[h]
    if(nh_vector[h]>0){
      meanh = mean(y[z[s,]==h])
      sigma2n = (v0*sigma20 + sum((y[z[s,]==h] - meanh)**2) + 
        k0*nh_vector[h]/kn*(meanh - mu0)**2)/vn
    }else{
      meanh = mu0; sigma2n = sigma20
    }
    sigma2[s,h] =  1/ rgamma(1, vn / 2, rate = vn / 2 * sigma2n)
    mu[s,h] = rnorm(1, (k0*mu0 + nh_vector[h]*meanh)/kn, sqrt(sigma2[s,h]/kn))
  }
  for (i in 1:n){
    prob = lambda[s,]*dnorm(y[i], mu[s,], sqrt(sigma2[s,]))
    z[s+1,i] = sample(H, 1, prob = prob)
  }
  setTxtProgressBar(pb, s)
}
close(pb)

library(T4cluster) #library to compute the posterior coclustering matrix
library(plot.matrix) #library to plot a matrix
PSM_estimate = psm(z[7000:10000,1:10])
plot(PSM_estimate[1:10,1:10])
PSM_true = psm(matrix(z_true[1:10], nrow = 1, ncol = 10))
plot(PSM_true)

z_est = salso(z[5000:10000,], loss="binder")
binder(z_true, z_est)  
rand.index(z_true, z_est)  

hist(y, breaks = 30)
stripchart(y,
           method = "jitter",
           pch = 23,
           bg = z_est,
           add = TRUE)

hist(y, breaks = 30)
stripchart(y,
           method = "jitter",
           pch = 23,
           bg = z_true,
           add = TRUE)


# === Posterior predictive density over histogram ===

# Use a subset of iterations (post burn-in)
idx <- 5000:10000  
# Grid for evaluation
yy <- seq(min(y), max(y), length.out = 500)

# Posterior density:
#   f̂(y) = (1/|idx|) * Σ_s Σ_h λ_{s,h} N(y | μ_{s,h}, σ²_{s,h})
dens_pred <- sapply(yy, function(yy0) {
  mean(
    rowSums(
      lambda[idx, ] *
        dnorm(yy0, mean = mu[idx, ], sd = sqrt(sigma2[idx, ]))
    )
  )
})

# Plot histogram with posterior density
hist(y, freq = FALSE, breaks = 30, col = "grey90",
     main = "Estimated Density", xlab = "y")
lines(yy, dens_pred, lwd = 2)
