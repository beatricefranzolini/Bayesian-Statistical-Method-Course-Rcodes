set.seed(0)
#Normal hierarchical partial exchangeable model

#num of groups and num of samples
J = 100 
n_j = rep(20, J)

#true data generating process
theta_true = runif(J, -1, 1)

#here we simulate n_dat datasets 
#to compare the MSE of freq and Bayes estimators
n_dat = 100
#PREPARE matrix TO SAVE POINT ESTIMATES 
theta_bayes = matrix(NA, nrow = n_dat, ncol = J)
theta_MLE = matrix(NA, nrow = n_dat, ncol = J)

#number of MCMC iterations for Bayes
S = 10000
# BUILD EMPTY MATRIX OF POSTERIOR SAMPLES
posterior.samples.theta = matrix(NA, nrow = S, ncol = J)

#the following lines create a progress bar 
pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = n_dat, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

#assume sigma2 known 
sigma2 = 10
for (num_dataset in 1:n_dat){
  #SIMULATE DATA 
  y = list()
  for (j in 1:J){
    y[[j]] = rnorm(n_j[j], theta_true[j], sqrt(sigma2))
  }
  
  #Compute the MLE  
  for(j in 1:J){
    theta_MLE[num_dataset, j] = mean(y[[j]])
  }
  
  #compute summary statistics
  S_j = rep(NA, J)
  for(j in 1:J){
    S_j[j] = sum(y[[j]])
  }
  
  #initialize the MCMC
  theta = theta_MLE[num_dataset, ]
  mu = 0
  
  t = 1
  posterior.samples.theta[t,] = theta
  
  for (t in 2: S){
    mu_post_vec = sigma2/n_j / (1 + sigma2/n_j) * mu + 
      1 / (1 + sigma2/n_j) * theta_MLE[num_dataset,]
    sigma_post_vec = 1 / (1 + n_j/sigma2)
    theta = rnorm(J, mu_post_vec, sqrt(sigma_post_vec))
    
    mu = rnorm(1, 0.1/(0.1+1/J)*mean(theta), sqrt(1 / (1/0.1 + J)) )
    posterior.samples.theta[t,] = theta
  }
  theta_bayes[num_dataset, ] = colMeans(posterior.samples.theta[5001:10000,])
  setTxtProgressBar(pb, num_dataset)
}

MSE_Bayes = mean(rowSums(t(apply(theta_bayes, 1, function(x) x-theta_true))^2))
MSE_Freq = mean(rowSums(t(apply(theta_MLE, 1, function(x) x-theta_true))^2))

#let's see what happens in the first group
hist(theta_MLE[,1], xlim = c(-3.5,3.5))
abline(v = theta_true[1], col = 'red', lwd = 2, lty = 'dashed')
abline(v = mean(theta_MLE[,1]), col = 'blue', lwd = 2, lty = 'dashed')

hist(theta_bayes[,1], xlim = c(-3.5,3.5))
abline(v = theta_true[1], col = 'red', lwd = 2, lty = 'dashed')
abline(v = mean(theta_bayes[,1]), col = 'blue', lwd = 2, lty = 'dashed')
