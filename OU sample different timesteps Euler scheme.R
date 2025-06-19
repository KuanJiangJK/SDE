rm(list = ls())

## OU Generation function
OU_generate <- function(mu, nu, sigma, timestep = 0.01, Npaths = 1e5, timemin = 0, timemax = 4, plot = FALSE, Nshow = NULL){
  # genearte diffusion term.
  Nlen <- as.integer(round((timemax - timemin)/timestep)) # there are Nlen points to be generated per path
  # add as.integer and round function because of floating number issue (who knows what exactly happened? something might CS scientist know)
  diffusion <- matrix(sqrt(timestep) * sigma * rnorm(Npaths*Nlen, 0, 1), ncol = Npaths, nrow = Nlen) # We need to generate n paths (each one column), and we each path Nlen points (rows)
  pmat <- matrix(0, nrow = Nlen + 1, ncol = Npaths)  #s tarting values for all paths at 0, Nlen + 1 is to make it starts from the second row for all paths
  for (i in 1:Nlen){ #For all pathes, accumulate datapoints based on the present value.
    pmat[i+1, ] = pmat[i, ] + (nu * (mu - pmat[i,])) * timestep + diffusion[i,]
  }
  t <- seq(timemin, timemax, timestep)
  if (plot == TRUE){
    if (is.null(Nshow)){
      Nshow <- (Npaths*0.001)+1
    }
    matplot(t, pmat[, sample(seq(1, Npaths, 1), Nshow)], col = grey(.5, .5), type = 'l', lty = 1, xlab = 't', ylab = expression('X'['t']), las = 1, bty = 'n')
  }
  return(list(time = t, paths = pmat))
}


# Test
# OU <- OU_generate(mu = 0.5, nu = 0.5, sigma = 0.3, timestep = 0.01, Npaths = 1e3, timemin = 0, timemax = 300, plot = TRUE, 50)


########################################### Generate data with different timesteps
# Same time range, same sample size, different timesteps.
# generate data with different timestep (50 cols of data with different timestep (but same sample size))
seq_timestep <- seq(0.01, 0.5, 0.01)
n_sample <- 2500
seq_timemax <- n_sample * seq_timestep # determines the time range (0 to seq_timemax)
OU_sample <- matrix(NA, nrow = n_sample+1, ncol = length(seq_timestep))
for (i in 1: length(seq_timestep)){
  timestep <- seq_timestep[i]
  timemax <- seq_timemax[i]
  OU <- OU_generate(mu = 0.5, nu = 0.5, sigma = 1, timestep, Npaths = 1, timemin = 0, timemax)
  OU_sample[,i] <- OU$paths
}

## likelihood calculation
## Assume an OU drift function, as well as a constant diffusion term. OU drift = nu * (mu - x)
OU_drift <- function(x, theta) {
  nu <- theta[1]
  mu <- theta[2]
  return(nu * (mu - x))
} # Drift part is expressed as b(xi-1, theta), theta here is c(nu, mu)

Euler_LL_OU <- function(theta, X, delta) { # delta is the timestep
  # theta = c(nu, mu), I don't need to define it here, the function in optim will bring it to the optimization 
  X0 <- X[1:(length(X) - 1)]
  X1 <- X[2:length(X)]
  b_vals <- OU_drift(X0, theta) ## drift part values
  loglik <- sum((X1 - X0) * b_vals) - 0.5 * delta * sum(b_vals^2)
  return(-loglik)  # negative for minimization
}


# Initialize storage for estimates
nu_estimates <- numeric(length(seq_timestep))
mu_estimates <- numeric(length(seq_timestep))

# Estimate parameters for each timestep
for (i in seq_along(seq_timestep)) {
  X <- OU_sample[, i]
  delta <- seq_timestep[i]
  # optimization, starting value at 0.4 and 0.4
  result <- optim(par = c(0.4, 0.4), fn = Euler_LL_OU, X = X, delta = delta)
  # Store results
  nu_estimates[i] <- result$par[1]
  mu_estimates[i] <- result$par[2]
}

# Plotting
plot(seq_timestep, nu_estimates, type = 'o', pch = 16, col = 'blue',
     ylim = range(c(nu_estimates, mu_estimates, 0.5)),
     xlab = "Timestep", ylab = "Parameter Estimate",
     main = "Estimates of OU Parameters by Timestep")
lines(seq_timestep, mu_estimates, type = 'o', pch = 16, col = 'darkgreen')
abline(h = 0.5, col = 'red', lty = 2)
legend("bottomright", legend = c("Estimated ν", "Estimated μ", "True Value"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 1, 2), pch = 16)


######################################################################################
# same time range, different sample size, different timesteps.
# Parameters
true_nu <- 0.5
true_mu <- 0.5
true_sigma <- 1
timemax <- 2000
seq_timestep <- seq(0.02, 5, 0.02)

# Containers
nu_estimates <- numeric(length(seq_timestep))
mu_estimates <- numeric(length(seq_timestep))

# Loop over timestep values
for (i in seq_along(seq_timestep)) {
  timestep <- seq_timestep[i]
  OU <- OU_generate(mu = true_mu, nu = true_nu, sigma = true_sigma,
                    timestep = timestep, timemin = 0, timemax = timemax,
                    Npaths = 1, plot = FALSE)
  X <- OU$paths[, 1]
  delta <- timestep
  result <- optim(par = c(0.4, 0.4), fn = Euler_LL_OU, X = X, delta = delta)
  nu_estimates[i] <- result$par[1]
  mu_estimates[i] <- result$par[2]
}

# --- PLOT ---

plot(seq_timestep, nu_estimates, type = 'o', pch = 16, col = 'blue',
     ylim = range(c(nu_estimates, mu_estimates, 0.5)),
     xlab = "Timestep (Δ)", ylab = "Parameter Estimate",
     main = "Estimates of OU Parameters by Timestep (Fixed Total Time)")
lines(seq_timestep, mu_estimates, type = 'o', pch = 16, col = 'darkgreen')
abline(h = 0.5, col = 'red', lty = 2)
legend("bottomright", legend = c("Estimated ν", "Estimated μ", "True Value"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 1, 2), pch = 16)

plot(seq_timestep, abs(nu_estimates-true_nu), type = 'o', pch = 16, col = 'blue',
     ylim = range(c(abs(nu_estimates- true_nu), abs(mu_estimates-true_mu), 0)),
     xlab = "Timestep (Δ)", ylab = "Parameter Estimate",
     main = "Estimates of OU Parameters by Timestep (Fixed Total Time)")
lines(seq_timestep, abs(mu_estimates-true_mu), type = 'o', pch = 16, col = 'darkgreen')
abline(h = 00, col = 'red', lty = 2)
legend("bottomright", legend = c("Estimated ν", "Estimated μ", "True Value"),
       col = c("blue", "darkgreen", "red"), lty = c(1, 1, 2), pch = 16)

