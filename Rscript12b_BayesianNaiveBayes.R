## ------------------------------------------------------------
## Bayesian Naive Bayes (Dirichlet–Categorical) for Discrete X
## Application: Email triage (important / other / junk)
## ------------------------------------------------------------

set.seed(0)

# --------------------------
# 0) Utility: discretization
# --------------------------
disc_by_breaks <- function(x, breaks, right = FALSE) {
  # returns integers 1..K for K bins
  cut(x, breaks = breaks, right = right, include.lowest = TRUE, labels = FALSE)
}

# -------------------------------------------------
# 1) Simulate data (convincing email triage example)
# -------------------------------------------------
simulate_email_data <- function(n = 2000) {
  classes <- c("important", "other", "junk")
  J <- length(classes)
  
  # Class prior (true, unknown to the learner)
  pi_true <- c(0.20, 0.55, 0.25)  # proportions for important/other/junk
  
  # Draw class labels
  y <- factor(sample(classes, size = n, replace = TRUE, prob = pi_true), levels = classes)
  
  # Features:
  # F1: has_links (binary: 1=yes, 0=no), class-dependent Bernoulli
  p_links <- c(0.05, 0.20, 0.75)  # fewer links in important, many in junk
  has_links <- rbinom(n, 1, prob = p_links[as.integer(y)])
  
  # F2: email_length (numeric words) -> discretize later
  # Poisson means by class (short for junk, long for important)
  mu_len <- c(220, 120, 60)
  email_length <- rpois(n, lambda = mu_len[as.integer(y)])
  
  # F3: has_attachment (binary)
  p_attach <- c(0.40, 0.10, 0.02)
  has_attachment <- rbinom(n, 1, prob = p_attach[as.integer(y)])
  
  # F4: sender_domain_type (categorical with K=4)
  # 1=whitelist colleague, 2=corporate, 3=newsletter, 4=random
  phi_domain <- rbind(
    c(0.55, 0.30, 0.10, 0.05),  # important
    c(0.10, 0.45, 0.30, 0.15),  # other
    c(0.02, 0.08, 0.25, 0.65)   # junk
  )
  K_dom <- ncol(phi_domain)
  sender_domain <- vapply(
    as.integer(y),
    function(j) sample.int(K_dom, size = 1, prob = phi_domain[j, ]),
    integer(1)
  )
  
  # F5: hour_bucket (categorical with 4 parts)
  # 1=night(0-6), 2=morning(6-12), 3=afternoon(12-18), 4=evening(18-24)
  phi_hour <- rbind(
    c(0.05, 0.45, 0.35, 0.15),  # important arrive mostly morning/afternoon
    c(0.10, 0.35, 0.35, 0.20),
    c(0.25, 0.20, 0.30, 0.25)   # junk somewhat nocturnal
  )
  hour_bucket <- vapply(
    as.integer(y),
    function(j) sample.int(4, size = 1, prob = phi_hour[j, ]),
    integer(1)
  )
  
  # Discretize email_length into bins (global, not per-class)
  # breaks: [0, 50), [50, 100), [100, 200), [200, +inf)
  len_breaks <- c(-Inf, 50, 100, 200, Inf)
  length_bin <- disc_by_breaks(email_length, breaks = len_breaks)
  
  # Assemble data.frame of discrete features as integers 1..Nd per feature
  # - Convert binary {0,1} to {1,2} for a clean Dirichlet–Categorical coding
  X <- data.frame(
    links = has_links + 1L,           # 1=no, 2=yes
    length_bin = length_bin,          # 1..4
    attach = has_attachment + 1L,     # 1=no, 2=yes
    domain = sender_domain,           # 1..4
    hour = hour_bucket                # 1..4
  )
  
  list(X = X, y = y, classes = classes,
       meta = list(len_breaks = len_breaks))
}

sim <- simulate_email_data(n = 2500)
X <- sim$X; y <- sim$y; classes <- sim$classes

# Train/test split
set.seed(1)
idx <- sample.int(nrow(X), size = round(0.7 * nrow(X)))
train <- list(X = X[idx, , drop = FALSE], y = y[idx])
test  <- list(X = X[-idx, , drop = FALSE], y = y[-idx])

# --------------------------------------------------
# 2) Fit Bayesian Naive Bayes with Dirichlet priors
# --------------------------------------------------
# Model:
#  - Class prior:    pi ~ Dirichlet(alpha_pi)
#  - For each feature d and class j:
#        theta_{d,j} ~ Dirichlet(alpha_{d,j})
#  - Posterior means give closed-form predictive factors:
#        E[pi_j | D] = (alpha_pi_j + N_j) / sum_j (alpha_pi_j + N_j)
#        E[theta_{d,j,k} | D] = (alpha_{d,j,k} + count_{d,j,k}) / sum_k (...)

fit_nb_dirichlet <- function(X, y, alpha_pi = NULL, alpha_feat = NULL) {
  X <- as.data.frame(X)
  classes <- levels(y)
  J <- length(classes)
  D <- ncol(X)
  
  # Determine Nd per feature and ensure values are in 1..Nd
  Nd <- vapply(X, function(col) max(col), integer(1))
  
  # Hyperparameters
  if (is.null(alpha_pi)) alpha_pi <- rep(1, J)             # flat Dirichlet prior on classes
  if (is.null(alpha_feat)) {
    # alpha_feat[[d]][j, k] = prior for feature d, class j, category k
    alpha_feat <- lapply(seq_len(D), function(d) {
      matrix(1, nrow = J, ncol = Nd[d])                    # flat per feature
    })
  }
  
  # Sufficient statistics
  y_int <- as.integer(y)  # 1..J
  N_y <- tabulate(y_int, nbins = J)
  
  counts <- vector("list", D)
  for (d in seq_len(D)) {
    Kd <- Nd[d]
    tab_d <- array(0L, dim = c(J, Kd))
    for (j in seq_len(J)) {
      xdj <- X[y_int == j, d]
      if (length(xdj)) {
        tab_d[j, ] <- tabulate(xdj, nbins = Kd)
      }
    }
    counts[[d]] <- tab_d
  }
  
  # Posterior hyperparameters
  alpha_pi_post <- alpha_pi + N_y
  alpha_feat_post <- lapply(seq_len(D), function(d) alpha_feat[[d]] + counts[[d]])
  
  list(
    classes = classes, J = J, D = D, Nd = Nd,
    alpha_pi_post = alpha_pi_post,
    alpha_feat_post = alpha_feat_post
  )
}

model <- fit_nb_dirichlet(train$X, train$y)

# --------------------------------------------------------------
# 3) Posterior-predictive classification for new data X_new (NxD)
# --------------------------------------------------------------
predict_nb <- function(model, X_new, return_proba = FALSE) {
  X_new <- as.data.frame(X_new)
  N <- nrow(X_new)
  J <- model$J; D <- model$D
  Nd <- model$Nd
  
  # Posterior means for priors (class and conditional feature distributions)
  pi_bar <- model$alpha_pi_post / sum(model$alpha_pi_post)
  
  # For each feature d and class j, probability vector over categories k:
  theta_bar <- lapply(seq_len(D), function(d) {
    post <- model$alpha_feat_post[[d]]
    post / rowSums(post)
  })
  # theta_bar[[d]][j, k]
  
  # Compute unnormalized post-predictive p(y=j | x) ∝ E[pi_j] * Π_d E[theta_{d,j, x_d}]
  log_scores <- matrix(0, nrow = N, ncol = J)
  for (j in seq_len(J)) {
    # start with log class prior
    log_scores[, j] <- rep(log(pi_bar[j]), N)
    for (d in seq_len(D)) {
      xk <- X_new[[d]]
      # Safety: clamp to [1, Nd[d]] in case of out-of-range
      xk <- pmax(1L, pmin(as.integer(xk), Nd[d]))
      probs_dj <- theta_bar[[d]][j, ]
      log_scores[, j] <- log_scores[, j] + log(probs_dj[xk])
    }
  }
  
  # Normalize to probabilities
  # softmax row-wise
  max_row <- apply(log_scores, 1, max)
  exp_shift <- exp(log_scores - max_row)
  denom <- rowSums(exp_shift)
  proba <- exp_shift / denom
  
  if (return_proba) {
    colnames(proba) <- model$classes
    return(proba)
  } else {
    y_hat <- factor(model$classes[max.col(proba, ties.method = "first")], levels = model$classes)
    return(y_hat)
  }
}

# -------------------------
# 4) Evaluate on the test set
# -------------------------
y_hat <- predict_nb(model, test$X)
acc <- mean(y_hat == test$y)
cat(sprintf("Test accuracy: %.3f\n", acc))

# Confusion matrix
print(table(predicted = y_hat, truth = test$y))

# ---------------------------------------------------------
# 5) Classify brand-new points (hand-crafted illustrations)
# ---------------------------------------------------------
# Helper to build a single new x row:
# features: links(1=no,2=yes), length_bin(1..4), attach(1=no,2=yes), domain(1..4), hour(1..4)
new_points <- rbind(
  # likely important: no links, long, has attachment, known domain, morning
  c(1, 4, 2, 1, 2),
  # likely junk: links, short, no attachment, random domain, night
  c(2, 1, 1, 4, 1),
  # ambiguous/other: maybe newsletter, medium length, no attach, afternoon
  c(1, 3, 1, 3, 3)
)
colnames(new_points) <- names(train$X)

pred_labels <- predict_nb(model, new_points, return_proba = FALSE)
pred_proba  <- predict_nb(model, new_points, return_proba = TRUE)

cat("\nPredictions for hand-crafted new points:\n")
print(pred_labels)
print(round(pred_proba, 3))

# ---------------------------------------------------------
# 6) (Optional) Compact helper that wraps train + predict
# ---------------------------------------------------------
train_nb <- function(X, y, alpha_pi = NULL, alpha_feat = NULL) {
  fit_nb_dirichlet(X, y, alpha_pi, alpha_feat)
}

## ============================================================
## 7) Example: Run Bayesian Naive Bayes on the iris dataset
## ============================================================

# Load iris data
data(iris)

library(ggplot2)
ggplot(iris, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "Iris dataset: Petal Length vs. Petal Width by Species") +
  theme_minimal()

# Convert Species to a factor (it already is)
y_iris <- iris$Species

# Discretize continuous predictors into bins
# Here we use 4 equal-frequency bins per feature (quantiles)
X_iris <- as.data.frame(lapply(iris[, 1:4], function(col) {
  cut(col, breaks = quantile(col, probs = seq(0, 1, length.out = 5)),
      include.lowest = TRUE, labels = FALSE)
}))

# Split into train/test (70/30)
set.seed(123)
idx_iris <- sample.int(nrow(X_iris), size = round(0.7 * nrow(X_iris)))
train_iris <- list(X = X_iris[idx_iris, ], y = y_iris[idx_iris])
test_iris  <- list(X = X_iris[-idx_iris, ], y = y_iris[-idx_iris])

# Fit the model
model_iris <- fit_nb_dirichlet(train_iris$X, train_iris$y)

# Predict on test set
y_hat_iris <- predict_nb(model_iris, test_iris$X)

# Accuracy
acc_iris <- mean(y_hat_iris == test_iris$y)
cat(sprintf("\nIris dataset accuracy: %.3f\n", acc_iris))

# Confusion matrix
print(table(predicted = y_hat_iris, truth = test_iris$y))

# Predict probabilities for a few new points (randomly chosen from test set)
new_idx <- sample(seq_len(nrow(test_iris$X)), 3)
new_points_iris <- test_iris$X[new_idx, ]
pred_labels_iris <- predict_nb(model_iris, new_points_iris, return_proba = FALSE)
pred_proba_iris  <- predict_nb(model_iris, new_points_iris, return_proba = TRUE)

cat("\nPredictions for a few iris test points:\n")
print(pred_labels_iris)
print(round(pred_proba_iris, 3))
