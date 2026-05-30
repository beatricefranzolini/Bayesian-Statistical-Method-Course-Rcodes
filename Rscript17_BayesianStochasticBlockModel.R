##############################
##  Libraries
##############################

# install.packages("igraph")
# install.packages("latentnet")
# install.packages("network")
# install.packages("salso")

library(igraph)
library(latentnet)
library(network)
library(salso)

##############################
##  Basic helpers
##############################

rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

##############################
##  Block statistics
##############################

## Undirected SBM: symmetric, no self-loops
## Given A (n x n), labels z (length n), K:
## returns m[h,j] = #edges between blocks h,j
##         M[h,j] = #possible edges between blocks h,j
sbm_block_stats_undirected <- function(A, z, K) {
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  
  m_mat <- matrix(0, K, K)
  M_mat <- matrix(0, K, K)
  
  for (h in 1:K) {
    idx_h <- which(z == h)
    n_h <- length(idx_h)
    if (n_h == 0) next
    
    for (j in h:K) {
      idx_j <- which(z == j)
      n_j <- length(idx_j)
      if (n_j == 0) next
      
      if (h == j) {
        subA <- A[idx_h, idx_h, drop = FALSE]
        m_hj <- sum(subA[upper.tri(subA)])
        M_hj <- n_h * (n_h - 1) / 2
      } else {
        subA <- A[idx_h, idx_j, drop = FALSE]
        m_hj <- sum(subA)
        M_hj <- n_h * n_j
      }
      
      m_mat[h, j] <- m_hj
      m_mat[j, h] <- m_hj
      M_mat[h, j] <- M_hj
      M_mat[j, h] <- M_hj
    }
  }
  
  list(m = m_mat, M = M_mat)
}

## Directed SBM: no symmetry, no self-loops
## m[h,j] = #edges from block h to block j
## M[h,j] = #possible directed dyads from block h to block j
sbm_block_stats_directed <- function(A, z, K) {
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  
  m_mat <- matrix(0, K, K)
  M_mat <- matrix(0, K, K)
  
  for (h in 1:K) {
    idx_h <- which(z == h)
    if (length(idx_h) == 0) next
    for (j in 1:K) {
      idx_j <- which(z == j)
      if (length(idx_j) == 0) next
      
      subA <- A[idx_h, idx_j, drop = FALSE]
      if (h == j) {
        # remove self-loops
        diag(subA) <- NA
      }
      m_hj <- sum(subA, na.rm = TRUE)
      M_hj <- sum(!is.na(subA))  # number of possible dyads
      
      m_mat[h, j] <- m_hj
      M_mat[h, j] <- M_hj
    }
  }
  
  list(m = m_mat, M = M_mat)
}

##############################
##  Gibbs sampler: undirected SBM
##############################

sbm_gibbs_undirected <- function(A, K,
                                 alpha = rep(1, K),
                                 a = 1, b = 1,
                                 n_iter = 5000,
                                 burn_in = 1000,
                                 thin = 10,
                                 verbose = TRUE) {
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  
  # initial labels
  z <- sample(1:K, n, replace = TRUE)
  
  n_save <- floor((n_iter - burn_in) / thin)
  Z_save <- matrix(NA_integer_, nrow = n, ncol = n_save)
  W_save <- array(NA_real_, dim = c(K, K, n_save))
  p_save <- matrix(NA_real_, nrow = K, ncol = n_save)
  iter_save <- 0
  
  for (it in 1:n_iter) {
    ## p | Z
    n_k <- tabulate(z, nbins = K)
    p <- rdirichlet(alpha + n_k)
    
    ## W | A,Z (use block stats)
    stats <- sbm_block_stats_undirected(A, z, K)
    m_mat <- stats$m
    M_mat <- stats$M
    W <- matrix(0, K, K)
    for (h in 1:K) {
      for (j in h:K) {
        m_hj <- m_mat[h, j]
        M_hj <- M_mat[h, j]
        if (M_hj == 0) {
          shape1 <- a
          shape2 <- b
        } else {
          shape1 <- a + m_hj
          shape2 <- b + M_hj - m_hj
        }
        W[h, j] <- rbeta(1, shape1, shape2)
        W[j, h] <- W[h, j]
      }
    }
    
    ## Z | rest
    for (u in 1:n) {
      log_q <- log(p)
      for (v in 1:n) {
        if (v == u) next
        a_uv <- A[u, v]
        z_v <- z[v]
        log_q <- log_q +
          a_uv * log(W[, z_v]) +
          (1 - a_uv) * log(1 - W[, z_v])
      }
      log_q <- log_q - max(log_q)
      q <- exp(log_q); q <- q / sum(q)
      z[u] <- sample(1:K, 1, prob = q)
    }
    
    if (it > burn_in && ((it - burn_in) %% thin == 0)) {
      iter_save <- iter_save + 1
      Z_save[, iter_save] <- z
      W_save[, , iter_save] <- W
      p_save[, iter_save] <- p
    }
    
    if (verbose && it %% 500 == 0) {
      cat("Undirected iteration", it, "of", n_iter, "\n")
    }
  }
  
  list(
    Z = Z_save,
    W = W_save,
    p = p_save,
    A = A,
    K = K,
    alpha = alpha,
    a = a,
    b = b,
    directed = FALSE
  )
}

##############################
##  Gibbs sampler: directed SBM
##############################

sbm_gibbs_directed <- function(A, K,
                               alpha = rep(1, K),
                               a = 1, b = 1,
                               n_iter = 5000,
                               burn_in = 1000,
                               thin = 10,
                               verbose = TRUE) {
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  
  z <- sample(1:K, n, replace = TRUE)
  
  n_save <- floor((n_iter - burn_in) / thin)
  Z_save <- matrix(NA_integer_, nrow = n, ncol = n_save)
  W_save <- array(NA_real_, dim = c(K, K, n_save))
  p_save <- matrix(NA_real_, nrow = K, ncol = n_save)
  iter_save <- 0
  
  for (it in 1:n_iter) {
    ## p | Z
    n_k <- tabulate(z, nbins = K)
    p <- rdirichlet(alpha + n_k)
    
    ## W | A,Z  (use directed block stats)
    stats <- sbm_block_stats_directed(A, z, K)
    m_mat <- stats$m
    M_mat <- stats$M
    W <- matrix(0, K, K)
    for (h in 1:K) {
      for (j in 1:K) {
        m_hj <- m_mat[h, j]
        M_hj <- M_mat[h, j]
        if (M_hj == 0) {
          shape1 <- a
          shape2 <- b
        } else {
          shape1 <- a + m_hj
          shape2 <- b + M_hj - m_hj
        }
        W[h, j] <- rbeta(1, shape1, shape2)
      }
    }
    
    ## Z | rest
    for (u in 1:n) {
      log_q <- log(p)
      for (v in 1:n) {
        if (v == u) next
        z_v <- z[v]
        
        # outgoing edge u -> v
        a_uv <- A[u, v]
        log_q <- log_q +
          a_uv * log(W[, z_v]) +
          (1 - a_uv) * log(1 - W[, z_v])
        
        # incoming edge v -> u
        a_vu <- A[v, u]
        log_q <- log_q +
          a_vu * log(W[z_v, ]) +
          (1 - a_vu) * log(1 - W[z_v, ])
      }
      log_q <- log_q - max(log_q)
      q <- exp(log_q); q <- q / sum(q)
      z[u] <- sample(1:K, 1, prob = q)
    }
    
    if (it > burn_in && ((it - burn_in) %% thin == 0)) {
      iter_save <- iter_save + 1
      Z_save[, iter_save] <- z
      W_save[, , iter_save] <- W
      p_save[, iter_save] <- p
    }
    
    if (verbose && it %% 500 == 0) {
      cat("Directed iteration", it, "of", n_iter, "\n")
    }
  }
  
  list(
    Z = Z_save,
    W = W_save,
    p = p_save,
    A = A,
    K = K,
    alpha = alpha,
    a = a,
    b = b,
    directed = TRUE
  )
}

##############################
##  SALSO-based summaries
##############################

## Undirected: SALSO partition + conjugate (p,W) given z_hat
sbm_posterior_summary_salso_undirected <- function(fit) {
  Z <- fit$Z   # n x S
  n <- nrow(Z)
  S <- ncol(Z)
  K <- fit$K
  A <- fit$A
  
  # SALSO partition: rows = draws, cols = items
  allocs <- t(Z)  # S x n
  z_hat <- as.integer(salso(allocs))
  
  # co-clustering matrix from original MCMC
  cocl <- matrix(0, n, n)
  for (s in 1:S) {
    z_s <- Z[, s]
    cocl <- cocl + outer(z_s, z_s, "==")
  }
  cocl <- cocl / S
  
  # conjugate posterior for p,W given z_hat
  n_k <- tabulate(z_hat, nbins = K)
  p_mean <- (fit$alpha + n_k) / sum(fit$alpha + n_k)
  
  stats <- sbm_block_stats_undirected(A, z_hat, K)
  m_mat <- stats$m
  M_mat <- stats$M
  
  W_mean <- (fit$a + m_mat) / (fit$a + fit$b + M_mat)
  W_mean[is.nan(W_mean)] <- fit$a / (fit$a + fit$b)  # if M=0, revert to prior mean
  
  list(
    z_hat  = z_hat,
    cocl   = cocl,
    p_mean = p_mean,
    W_mean = W_mean,
    m_mat  = m_mat,
    M_mat  = M_mat
  )
}

## Directed: SALSO partition + conjugate (p,W) given z_hat
sbm_posterior_summary_salso_directed <- function(fit) {
  Z <- fit$Z   # n x S
  n <- nrow(Z)
  S <- ncol(Z)
  K <- fit$K
  A <- fit$A
  
  # SALSO partition
  allocs <- t(Z)  # S x n
  z_hat <- as.integer(salso(allocs))
  
  # co-clustering matrix
  cocl <- matrix(0, n, n)
  for (s in 1:S) {
    z_s <- Z[, s]
    cocl <- cocl + outer(z_s, z_s, "==")
  }
  cocl <- cocl / S
  
  # conjugate posterior for p,W given z_hat
  n_k <- tabulate(z_hat, nbins = K)
  p_mean <- (fit$alpha + n_k) / sum(fit$alpha + n_k)
  
  stats <- sbm_block_stats_directed(A, z_hat, K)
  m_mat <- stats$m
  M_mat <- stats$M
  
  W_mean <- (fit$a + m_mat) / (fit$a + fit$b + M_mat)
  W_mean[is.nan(W_mean)] <- fit$a / (fit$a + fit$b)
  
  list(
    z_hat  = z_hat,
    cocl   = cocl,
    p_mean = p_mean,
    W_mean = W_mean,
    m_mat  = m_mat,
    M_mat  = M_mat
  )
}

##############################
##  Application 1: Zachary's karate club (undirected)
##############################

cat("=== Zachary's karate club ===\n")

g_karate <- make_graph("Zachary")
A_karate <- as.matrix(as_adjacency_matrix(g_karate))
diag(A_karate) <- 0

set.seed(123)
fit_karate <- sbm_gibbs_undirected(
  A = A_karate,
  K = 2,          # classic choice
  alpha = rep(1, 2),
  a = 1, b = 1,
  n_iter = 4000,
  burn_in = 2000,
  thin = 10,
  verbose = TRUE
)

sum_karate <- sbm_posterior_summary_salso_undirected(fit_karate)

cat("SALSO labels (karate):\n")
print(sum_karate$z_hat)
cat("Posterior mean W (karate):\n")
print(round(sum_karate$W_mean, 3))
cat("Posterior mean p (karate):\n")
print(round(sum_karate$p_mean, 3))

# quick plot
V(g_karate)$block <- factor(sum_karate$z_hat)
plot(g_karate, vertex.color = V(g_karate)$block,
     main = "Zachary's karate club - SBM blocks")

##############################
##  Application 2: UKfaculty (directed)
##############################

cat("\n=== UKfaculty (directed) ===\n")

data("UKfaculty")
A_uk <- as.matrix(as_adjacency_matrix(UKfaculty))
diag(A_uk) <- 0

set.seed(123)
fit_uk <- sbm_gibbs_directed(
  A = A_uk,
  K = 3,          # example choice; can adjust
  alpha = rep(1, 3),
  a = 1, b = 1,
  n_iter = 5000,
  burn_in = 2500,
  thin = 10,
  verbose = TRUE
)

sum_uk <- sbm_posterior_summary_salso_directed(fit_uk)

cat("SALSO labels (UKfaculty):\n")
print(sum_uk$z_hat)
cat("Posterior mean W (UKfaculty):\n")
print(round(sum_uk$W_mean, 3))
cat("Posterior mean p (UKfaculty):\n")
print(round(sum_uk$p_mean, 3))

# simple visualisation: collapse directions to undirected for plotting
g_uk <- graph_from_adjacency_matrix(A_uk)
V(g_uk)$block <- factor(sum_uk$z_hat)
plot(g_uk, vertex.color = V(g_uk)$block,
     main = "UKfaculty - SBM blocks (SALSO, plotted undirected)")
