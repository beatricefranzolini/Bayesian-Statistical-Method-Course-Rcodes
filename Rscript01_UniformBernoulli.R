#Uniform-Bernoulli model, coin tossing 
rm(list = ls())

####uniform prior for theta
#plot prior distribution
#define x-axis
x = seq(0, 1, length=500)

#calculate uniform distribution probabilities
p = dunif(x, min = 0, max = 1)

#plot uniform distribution
plot(x, p, type = 'l', lwd = 3,
     xlab='theta', ylab='p(theta)', 
     main='Prior on theta')

####simulate n=3 data
theta_true = 0.3

set.seed(1)
n = 3
y = rbinom(n = n, size = 1, prob = theta_true)


####posterior distribution 
####(I know the expression because I have compute it analytically)
#plot the posterior distribution 
#calculate uniform distribution probabilities
p = dbeta(x, shape1 = 1 + sum(y), shape2 = 1 + n - sum(y))

#plot uniform distribution
plot(x, p, type = 'l',  lwd = 3,
     xlab='theta', ylab='p(theta|y1,y2,y3)', 
     main='Posterior on theta n = 3')
abline(v = c((1 + sum(y))/(2+n)), col = c("green"), lwd = c(3))
p_MLE = c(mean(y))
asym_MLE_lb_CI = p_MLE-1.96*sqrt(p_MLE*(1-p_MLE)/n)
asym_MLE_ub_CI = p_MLE+1.96*sqrt(p_MLE*(1-p_MLE)/n)
abline(v = c(asym_MLE_lb_CI, p_MLE, asym_MLE_ub_CI), 
       col=c("blue", "red", "blue"), lty=c(1, 2, 1), lwd=c(1, 3, 1))

###simulate n=7 additional data
set.seed(3)
n = 10
y = c(y, rbinom(n = 7, size = 1, prob = theta_true))

####posterior distribution 
####(I know the expression because I have compute it analytically)
#plot the posterior distribution 
#calculate uniform distribution probabilities
p = dbeta(x, shape1 = 1 + sum(y), shape2 = 1 + n - sum(y))

#plot uniform distribution
plot(x, p, type = 'l',  lwd = 3,
     xlab='theta', ylab='p(theta|y1,...,y10)', 
     main='Posterior on theta n = 10')
abline(v = c((1 + sum(y))/(2+n)), col = c("green"), lwd = c(3))
p_MLE = c(mean(y))
asym_MLE_lb_CI = p_MLE-1.96*sqrt(p_MLE*(1-p_MLE)/n)
asym_MLE_ub_CI = p_MLE+1.96*sqrt(p_MLE*(1-p_MLE)/n)
abline(v = c(asym_MLE_lb_CI, p_MLE, asym_MLE_ub_CI), 
       col=c("blue", "red", "blue"), lty=c(1, 2, 1), lwd=c(1, 3, 1))

#simulate 90 additional data 
n = 100
y = c(y, rbinom(n = 90, size = 1, prob = theta_true))

####posterior distribution 
####(I know the expression because I have compute it analytically)
#plot the posterior distribution 
#calculate uniform distribution probabilities
p = dbeta(x, shape1 = 1 + sum(y), shape2 = 1 + n - sum(y))

#plot uniform distribution
plot(x, p, type = 'l',  lwd = 3,
     xlab='theta', ylab='p(theta|y1,...,y100)', 
     main='Posterior on theta n = 100')
abline(v = c((1 + sum(y))/(2+n)), col = c("green"), lwd = c(3))
p_MLE = c(mean(y))
asym_MLE_lb_CI = p_MLE-1.96*sqrt(p_MLE*(1-p_MLE)/n)
asym_MLE_ub_CI = p_MLE+1.96*sqrt(p_MLE*(1-p_MLE)/n)
abline(v = c(asym_MLE_lb_CI, p_MLE, asym_MLE_ub_CI), 
       col=c("blue", "red", "blue"), lty=c(1, 2, 1), lwd=c(1, 3, 1))

#simulate 9000 additional data
n = 10000
y = c(y, rbinom(n = 9900, size = 1, prob = theta_true))

####posterior distribution 
####(I know the expression because I have compute it analytically)
#plot the posterior distribution 
#calculate uniform distribution probabilities
p = dbeta(x, shape1 = 1 + sum(y), shape2 = 1 + n - sum(y))

#plot uniform distribution
plot(x, p, type = 'l',  lwd = 3,
     xlab='theta', ylab='p(theta|y1,...,y10000)', 
     main='Posterior on theta n = 10000')
abline(v = c((1 + sum(y))/(2+n)), col = c("green"), lwd = c(3))
p_MLE = c(mean(y))
asym_MLE_lb_CI = p_MLE-1.96*sqrt(p_MLE*(1-p_MLE)/n)
asym_MLE_ub_CI = p_MLE+1.96*sqrt(p_MLE*(1-p_MLE)/n)
abline(v = c(asym_MLE_lb_CI, p_MLE, asym_MLE_ub_CI), 
       col=c("blue", "red", "blue"), lty=c(1, 2, 1), lwd=c(1, 3, 1))
