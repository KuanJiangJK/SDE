rm(list = ls())

# Model parameters
c1 <- 1 
c2 <- 1
c3 <- 1
mu <- 6 # 
sigma <- 0.1 # diffusion
tau <- 5 # discrete delay in time steps
Dt <- 0.1
T <- 30
t <- seq(0, T, Dt)
Npaths <- 1e2
Nlen <- T / Dt

# Initialize matrix for Y
pmat_y <- matrix(0, ncol = Npaths, nrow = Nlen + 1)

# Simulate from tau+2 to Nlen (from the 1st to the tau+1'th of y,  y is assumed be 0)
for (j in (tau + 2):Nlen) {
  # delayed index is j - tau - 1
  delayed_index <- j - tau - 1
  delay_term <- c1 * sin(c2 * t[delayed_index] + c3) + pmat_y[delayed_index, ]
  noise_term <- sigma * rnorm(Npaths, 0, sqrt(Dt))
  pmat_y[j + 1, ] <- pmat_y[j, ] + mu * delay_term * Dt^2 + noise_term
}

# Plotting first 3 paths
par(mfrow = c(1, 1))
Nshow <- 3
colors <- c("blue", "red", "green")
matplot(t,
        pmat_y[, 1:Nshow],
        type = 'l', lty = 1, col = colors,
        xlab = 't', ylab = 'Y',
        las = 1, bty = 'n')

legend("bottomright",
       legend = c("y path 1", "y path 2", "y path 3"),
       col = colors,
       lty = 1,
       cex = 0.5,
       bty = "n")

# Define log-likelihood for a single path
LL_delay <- function(mu, y, Dt, tau, c1, c2, c3) {
  n <- length(y)
  idx_set <- (tau + 2):(n - 1) # valid indices (t in t+1)
  
  Ft <- c1 * sin(c2 * t[idx_set - tau - 1] + c3) + y[idx_set - tau - 1]
  dY <- y[idx_set + 1] - y[idx_set]
  mean_term <- mu * Ft * Dt^2
  
  sigma2_hat <- sum((dY - mean_term)^2) / (length(dY) * Dt)
  logLik <- -0.5 * length(dY) * log(2 * pi * sigma2_hat * Dt) - sum((dY - mean_term)^2) / (2 * sigma2_hat * Dt)
  return(logLik)
}

# Wrapper for vectorized optimization
LL_delay_wrapper <- function(mu, y, Dt, tau, c1, c2, c3) {
  if (length(mu) > 1) {
    sapply(mu, LL_delay, y = y, Dt = Dt, tau = tau, c1 = c1, c2 = c2, c3 = c3)
  } else {
    LL_delay(mu, y, Dt, tau, c1, c2, c3)
  }
}

#######################################
tau_candidates <- 1:30 # a vector of candidate tau to choose from
logLik_vec_total <- numeric(length(tau_candidates)) # ll store
estimate_mle <- array(dim=(c(length(tau_candidates),2,Npaths))) # mle estimate store
for (k in seq_along(tau_candidates)) {
  tau_k <- tau_candidates[k]
  loglik_tau_k <- numeric(Npaths) # 0 vector with lenght of Npaths
  
  for (i in 1:Npaths) {
    yi <- pmat_y[, i] # the i th path chosen from the pmat
    idx_set <- (tau_k + 2):(length(yi) - 1) # index of the effective yi
    if (length(idx_set) < 10) next # skips the current iteration of the inner for loop if there are fewer than 10 valid time points
    
    # Use likelihood-based optimization for mu
    opt <- optimize(LL_delay_wrapper, interval = c(0.01, 30), maximum = TRUE,
                    y = yi, Dt = Dt, tau = tau_k, c1 = 1, c2 = 1, c3 = 1)
    
    mu_hat <- opt$maximum
    
    # Compute Ft and deltaY with optimized mu
    Ft <- c1 * sin(c2 * t[idx_set - tau_k - 1] + c3) + yi[idx_set - tau_k - 1]
    dY <- yi[idx_set + 1] - yi[idx_set]
    
    # sigma2
    sigma2_hat <- sum((dY - mu_hat * Ft * Dt^2)^2) / (length(dY) * Dt)
    
    # calclate the likelihood
    ll <- -0.5 * length(dY) * log(2 * pi * sigma2_hat * Dt) -
      sum((dY - mu_hat * Ft * Dt^2)^2) / (2 * sigma2_hat * Dt)
    
    estimate_mle[k, 1, i] <- mu_hat
    estimate_mle[k, 2, i] <- sigma2_hat
    loglik_tau_k[i] <- ll
    # loglik_tau_k[i] <- ll / length(idx_set)
  }
  
  logLik_vec_total[k] <- mean(loglik_tau_k, na.rm = TRUE)
}


plot(tau_candidates, logLik_vec_total, type = "b", col = "blue", pch = 19,
     xlab = expression(tau), ylab = "Total Log-Likelihood",
     main = "Log-Likelihood Profile over Candidate Delays")
best_tau <- tau_candidates[which.max(logLik_vec_total)]
abline(v = best_tau, col = "red", lty = 2)
legend("bottomright", legend = paste("Best tau =", best_tau),
       col = "red", lty = 2, bty = "n")

par(mfrow = c(1, 2))
hist(estimate_mle[5,1,], breaks = 20, col = "lightblue", main = "mu", xlab = "")
hist(estimate_mle[5,2,], breaks = 20, col = "lightgreen", main = "sigma^2", xlab = "")
