library(ggplot2)
library(gridExtra)

#A. Dirichlet-categorical model (one sample)
#simulate data for n = 20 subjects
set.seed(1234)
m = 10 #number of categories
n = 20 #number of data points
truth = data.frame( names = factor(1:m), 
                    values = rep(1/m,m)) #true frequencies
y = sample( m, n, replace = TRUE, prob = truth$values) #data from a uniform distribution

#set Dirichlet prior 
alpha = rep(2,m) #Dirichlet-prior parameter 
alpha0 = sum(alpha)

#summary statistics: frequencies / MLE estimate
obs_freq = as.vector(table(factor(y, levels = 1:m)))
c = data.frame(names = factor(1:m), 
               values = obs_freq / n,
               sd = obs_freq/n*(1 - obs_freq/n) / n)

#Bayesian posterior summaries 
alpha_n = obs_freq+alpha #posterior parameter
alpha0_n = sum(alpha_n) #=alpha0+n
Bayes = data.frame(names=factor(1:m), 
          postparam = alpha_n,
          values = alpha_n / alpha0_n)

#plot results 
#par(mfrow=c(3,1)) 
truth_plot = ggplot(truth) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.4) + theme_classic()+ ggtitle("TRUTH")

MLE_plot = ggplot(c) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.4) +
  geom_errorbar( aes(x=names, ymin=pmax(rep(0,m),values-1.96*sd), ymax=values+1.96*sd), width=0.4, 
                 colour="red", alpha=0.9, linewidth=1.3) + theme_classic()+ ggtitle("MLE")+
  geom_abline(slope=0, intercept=0.1,  col = "blue",lty=2)

Bay_plot = ggplot(Bayes) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.4) +
  geom_errorbar( aes(x=names, ymin=qbeta(0.05,postparam,sum(postparam)-postparam),
                     ymax=qbeta(0.95,postparam,sum(postparam)-postparam)), width=0.4, 
                 colour="orange", alpha=0.9, linewidth=1.3) + theme_classic()+ ggtitle("Bayes")+
  geom_abline(slope=0, intercept=0.1,  col = "blue",lty=2)


grid.arrange(truth_plot, MLE_plot, Bay_plot)


#simulate data for n = 1000 subjects
set.seed(1234)
m = 10 #number of categories
n = 1000 #number of data points
truth = data.frame( names = factor(1:m), 
                    values = rep(1/m,m)) #true frequencies
y = sample( m, n, replace = TRUE, prob = truth$values) #data from a uniform distribution

#set Dirichlet prior 
alpha = rep(2,m) #Dirichlet-prior parameter 
alpha0 = sum(alpha)

#summary statistics: frequencies / MLE estimate
obs_freq = as.vector(table(factor(y, levels = 1:m)))
c = data.frame(names = factor(1:m), 
               values = obs_freq / n,
               sd = obs_freq/n*(1 - obs_freq/n) / n)

#Bayesian posterior summaries 
alpha_n = obs_freq+alpha #posterior parameter
alpha0_n = sum(alpha_n) #=alpha0+n
Bayes = data.frame(names=factor(1:m), 
                   postparam = alpha_n,
                   values = alpha_n / alpha0_n)

#plot results 
#par(mfrow=c(3,1)) 
truth_plot = ggplot(truth) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.15) + theme_classic()+ ggtitle("TRUTH")

MLE_plot = ggplot(c) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.15) +
  geom_errorbar( aes(x=names, ymin=pmax(rep(0,m),values-1.96*sd), ymax=values+1.96*sd), width=0.4, 
                 colour="red", alpha=0.9, linewidth=1.3) + theme_classic()+ ggtitle("MLE")+
  geom_abline(slope=0, intercept=0.1,  col = "blue",lty=2)

Bay_plot = ggplot(Bayes) +
  geom_bar( aes(x=names, y=values), stat="identity", fill="skyblue", alpha=0.7 ) +
  ylim(0, 0.15) +
  geom_errorbar( aes(x=names, ymin=qbeta(0.05,postparam,sum(postparam)-postparam),
                     ymax=qbeta(0.95,postparam,sum(postparam)-postparam)), width=0.4, 
                 colour="orange", alpha=0.9, linewidth=1.3) + theme_classic()+ ggtitle("Bayes")+
  geom_abline(slope=0, intercept=0.1,  col = "blue",lty=2)


grid.arrange(truth_plot, MLE_plot, Bay_plot)
