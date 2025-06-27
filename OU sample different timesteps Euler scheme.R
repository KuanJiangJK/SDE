rm(list = ls())

# True parameters
nu <- 1
mu <- 0.5
sigma <- 0.5
timemin <- 0
timemax <- 20

# Simulation settings
seq_timestep <- seq(0.001, 0.5, 0.001)
timestep <- 0.001 # Time step being used to generate the first set of the OU process
t <- seq(timemin, timemax, timestep)
x <- numeric(length(t))
x[1] <- 0 # initial value

# Simulate one path using the exact OU transition
for (i in 1:(length(t) - 1)) {
  x[i + 1] <- rnorm(1,
                    mean = x[i] * exp(-nu * timestep) + mu * (1 - exp(-nu * timestep)),
                    sd = sqrt((sigma^2 * 0.5 / nu) * (1 - exp(-2 * nu * timestep)))
  )
}

# Plot path
plot(t, x, col = grey(0.5, 0.5), type = "l")
abline(h = mu, col = "red")

# Drift and Euler functions
OU_drift <- function(x, theta) theta[1] * (theta[2] - x)

dcEuler <- function(x, dt, x0, theta, drift) {
  dd <- drift(x0, theta)
  (x - x0) * dd - 0.5 * dt * dd^2
}

# Estimation results
nu_euler <- numeric(length(seq_timestep)) # to store results from the Euler approximation
mu_euler <- numeric(length(seq_timestep)) # to store results from the Euler approximation
nu_exact <- numeric(length(seq_timestep)) # to store results from the exact likelihood
mu_exact <- numeric(length(seq_timestep)) # to store results from the exact likelihood

# Loop over timesteps
for (i in seq_along(seq_timestep)) {
  dt <- seq_timestep[i] # the current timestep
  step <- round(dt / timestep) # steps needed to extract one subsample
  idx <- seq(1, length(x), by = step)
  x_sub <- x[idx] # extract the subsample based on the indexation
  n <- length(x_sub) # sample size
  
  # ---- Euler pseudo-likelihood ----
  negloglik_euler <- function(par) {
    theta <- par
    -sum(dcEuler(x_sub[2:n], dt, x_sub[1:(n - 1)], theta, OU_drift))
  }
  
  result_euler <- tryCatch( 
    optim(par = c(1.5, 1), fn = negloglik_euler,
          method = "L-BFGS-B", lower = c(0, 0)),
    error = function(e) list(par = c(NA, NA))
  )
  nu_euler[i] <- result_euler$par[1]
  mu_euler[i] <- result_euler$par[2]
  
  # ---- Exact transition likelihood ----
  negloglik_exact <- function(par) {
    theta_nu <- par[1]
    theta_mu <- par[2]
    mean_vec <- x_sub[1:(n - 1)] * exp(-theta_nu * dt) +
      theta_mu * (1 - exp(-theta_nu * dt))
    var <- (sigma^2 / (2 * theta_nu)) * (1 - exp(-2 * theta_nu * dt))
    
    -sum(dnorm(x_sub[2:n], mean = mean_vec, sd = sqrt(var), log = TRUE)) ## from the true transition distribution, after log-transformation, and summed up.
  }
  
  result_exact <- tryCatch(
    optim(par = c(1.5, 1), fn = negloglik_exact,
          method = "L-BFGS-B", lower = c(1e-6, 0)),
    error = function(e) list(par = c(NA, NA))
  )
  nu_exact[i] <- result_exact$par[1]
  mu_exact[i] <- result_exact$par[2]
}

# ---- Plot comparison ----
plot(seq_timestep, nu_euler, type = "l", col = "blue",
     ylim = range(0,2),
     xlab = "Timestep", ylab = "Estimated Parameters",
     main = "Euler vs Exact Likelihood Estimates")
lines(seq_timestep, mu_euler, col = "red")
lines(seq_timestep, nu_exact, col = "blue", lty = 2)
lines(seq_timestep, mu_exact, col = "red", lty = 2)
abline(h = c(nu, mu), col = c("blue", "red"), lty = 3)
legend("bottomright",
       legend = c("nu (Euler)", "mu (Euler)", "nu (Exact)", "mu (Exact)", "True nu", "True mu"),
       col = c("blue", "red", "blue", "red", "blue", "red"), # colors
       lty = c(1, 1, 2, 2, 3, 3), # line types 
       bty = "n", # o or n.type of legend box. n is better, transparent.
       cex = 0.5) # size of legend

