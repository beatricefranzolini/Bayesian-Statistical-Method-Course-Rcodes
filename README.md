# Bayesian Statistical Methods Course R Codes

This repository contains the R scripts used in a Bayesian statistical methods course. 
The examples are organized as standalone scripts that progress from elementary conjugate Bayesian models and Monte Carlo simulation to Gibbs sampling, Metropolis-Hastings, finite mixture models, Gaussian processes, Bayesian classification, and Bayesian network models.

The code is intended for learning, experimentation, and classroom demonstration. 
Most scripts simulate their own data, fit or sample from a Bayesian model, and then visualize the prior, posterior, Markov chain output, or predictive behavior.

**Author**: [Beatrice Franzolini](https://beatricefranzolini.github.io/)

## Repository contents

| Script | Main topic | What it demonstrates |
| --- | --- | --- |
| `Rscript01_UniformBernoulli.R` | Uniform-Bernoulli model | Bayesian updating for a coin-tossing model with a uniform prior, beta posterior, posterior mean, MLE, and asymptotic confidence intervals as sample size grows. |
| `Rscript02_MonteCarloforPi.R` | Monte Carlo estimation | Approximation of $\pi$ by simulating points in a square and estimating the area of the unit circle. |
| `Rscript03_MonteCarlo_invCDF_exponential.R` | Inverse-CDF sampling | Simulation from an exponential distribution using inverse transform sampling and comparison with the analytical density. |
| `Rscript04_RejectionSampling_BetaBernoulli.R` | Rejection sampling | Sampling from a beta posterior in a Bernoulli model by accepting prior-likelihood simulations matching an observed Bernoulli outcome. |
| `Rscript05_RejectionSampling_BetaBinomial.R` | Rejection sampling | Extension of the previous accept/reject idea to a beta-binomial posterior with observed binomial data. |
| `Rscript06_GibbsSamplingForNormalModel.R` | Gibbs sampling | Gibbs sampler for a univariate normal model with unknown mean and variance, including trace plots, burn-in removal, posterior summaries, and density plots. |
| `Rscript07_GibbsSamplingForMultivariateNormal.R` | Gibbs sampling | Gibbs sampler for a bivariate normal target distribution, including conditional updates and convergence diagnostics. |
| `Rscript08_GibbsSamplingChangepointPoisson.R` | Changepoint model | Gibbs sampling for a Poisson changepoint model with two rates and an unknown changepoint. |
| `Rscript09_GibbsSamplingToyclustering.R` | Bayesian clustering | Gibbs sampling for a toy clustering example and visualization of posterior co-clustering probabilities. |
| `Rscript10_BiasVarianceTradeOff.R` | Bias-variance trade-off | Simulation study comparing polynomial regression fits and illustrating model complexity, bias, and variance. |
| `Rscript11_MetropolisHastings.R` | Metropolis-Hastings | Random-walk Metropolis-Hastings examples for posterior simulation, including proposal tuning and acceptance behavior. |
| `Rscript12_CategoricalDirichlet.R` | Categorical-Dirichlet model | Conjugate Bayesian analysis for categorical outcomes using a Dirichlet prior and posterior visualization. |
| `Rscript12b_BayesianNaiveBayes.R` | Bayesian Naive Bayes | Dirichlet-multinomial Bayesian Naive Bayes implementation with posterior predictive probabilities and an iris-data example. |
| `Rscript13_FiniteMixturemodel.R` | Finite mixture model | Gibbs sampler for a normal finite mixture model with latent allocations, Dirichlet mixture weights, posterior clustering summaries, and posterior predictive density. |
| `Rscript14_GaussianProcessesPriors.R` | Gaussian process priors | Simulation of Gaussian process prior paths under several mean and covariance kernel choices. |
| `Rscript15_GaussianProcessesRegression.R` | Gaussian process regression | Gaussian process regression with a squared-exponential plus periodic kernel, posterior mean, posterior covariance, and posterior sample paths. |
| `Rscript16_NetworkData.R` | Network data simulation | Simulation and visualization of Erdős-Rényi, preferential attachment, and stochastic block model networks. |
| `Rscript17_BayesianStochasticBlockModel.R` | Bayesian stochastic block models | Gibbs samplers for undirected and directed stochastic block models, posterior partition summaries, and applications to example network datasets. |


## Questions or bug reports

For questions or bug reports, please contact:

- [Beatrice Franzolini](https://beatricefranzolini.github.io/) — franzolini@pm.me

## Learning path

A suggested order is to follow the script numbering:

1. **Conjugate Bayesian models and basic simulation**: scripts 01-05 introduce posterior updating, Monte Carlo integration, inverse-CDF sampling, and rejection sampling.
2. **Markov chain Monte Carlo**: scripts 06-11 introduce Gibbs sampling, changepoint models, clustering with latent variables, and Metropolis-Hastings.
3. **Bayesian classification and mixture modeling**: scripts 12, 12b, and 13 cover categorical-Dirichlet models, Bayesian Naive Bayes, and finite mixtures.
4. **Complex and structured models**: scripts 14-17 cover Gaussian processes and Bayesian models for network data.

## Requirements

### R

The scripts are written in base R plus a small set of CRAN packages. Use a recent R 4.x installation when possible.

### R packages

Several scripts rely only on base R. The following packages are used by at least one script:

- `coda`
- `fossil`
- `ggplot2`
- `gridExtra`
- `igraph`
- `latentnet`
- `latex2exp`
- `MASS`
- `MCMCprecision`
- `mvtnorm`
- `network`
- `plot.matrix`
- `salso`
- `T4cluster`

You can install the package dependencies with:

```r
install.packages(c(
  "coda",
  "fossil",
  "ggplot2",
  "gridExtra",
  "igraph",
  "latentnet",
  "latex2exp",
  "MASS",
  "MCMCprecision",
  "mvtnorm",
  "network",
  "plot.matrix",
  "salso",
  "T4cluster"
))
```

## How to run the scripts

Clone the repository and open it as your working directory in RStudio.

Many scripts create plots. 

- Most simulation scripts set a random seed with `set.seed(...)`, so repeated runs should produce the same simulated data or MCMC trajectories within the same R/package environment.
- The scripts are standalone teaching examples rather than a single R package. They clear or create objects in the global environment and are best run one at a time.
- Some MCMC and Gaussian process scripts can take longer than the introductory Monte Carlo examples because they draw many posterior samples or manipulate larger covariance matrices.
- The code emphasizes transparency for instruction. In several places, loops are used deliberately so that the algorithmic steps are easy to follow.

## Script groups and key ideas

### 1. Conjugacy and basic Monte Carlo

The first scripts introduce simple probabilistic simulation and closed-form posterior updating:

- Bernoulli likelihood with uniform/beta posterior updates.
- Monte Carlo approximation by empirical averages.
- Inverse-transform sampling from a known CDF.
- Rejection sampling by simulating from a prior predictive mechanism and retaining parameter values that match observed data.

### 2. Gibbs sampling and Metropolis-Hastings

The MCMC examples show how posterior distributions can be explored when direct sampling is difficult:

- Conditional posterior updates for normal models.
- Latent-variable updates for changepoint and clustering models.
- Markov chain trace plots, burn-in handling, posterior summaries, and density estimates.
- Proposal distributions and acceptance probabilities in Metropolis-Hastings.

### 3. Classification and mixture models

The classification and mixture scripts introduce Bayesian models with categorical latent structure:

- Dirichlet priors for categorical probabilities.
- Bayesian Naive Bayes with posterior predictive class probabilities.
- Finite mixture modeling with latent cluster assignments and posterior co-clustering summaries.

### 4. Gaussian processes and networks

The final scripts demonstrate more structured Bayesian models:

- Gaussian process priors as distributions over functions.
- Gaussian process regression with posterior mean and posterior uncertainty.
- Random graph simulation under Erdős-Rényi, preferential attachment, and stochastic block models.
- Bayesian stochastic block model inference for undirected and directed network data.

## Troubleshooting

### `there is no package called ...`

Install the missing package with:

```r
install.packages("packageName")
```

Then rerun the script.
3. Set a random seed for examples that depend on simulation.
4. Update this README with the new script and any new package dependency.
