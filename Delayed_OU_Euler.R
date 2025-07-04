rm(list = ls())

# Model parameters for OU
mu <- 1.5
nu <- .6
sigma <- 1.5
k <- 5 # time delayed 

# Simulation parameters
Dt <- 0.001
T <- 40 # Time up to which paths should be simulated
t <- seq(0, T, Dt)
Npaths <- 1e2 # number of sample paths to simulate
Nlen <- T/Dt

set.seed(123)
# generate delayed OU process based on Euler discretization
pmat <- matrix(0, ncol = Npaths, nrow = Nlen+1)
for (i in 1: k){ ## for the first k observation, they are affected by X0 = 0
   pmat[i+1, ] <- pmat[i, ]+ mu*nu*Dt + rnorm(Npaths, 0, sigma*sqrt(Dt)) ## first k only affected by X0.
}
for(j in (k+1):Nlen){
  pmat[j+1, ] <- pmat[j, ] + mu*(nu - pmat[j-k, ])*Dt + rnorm(Npaths, 0, sigma*sqrt(Dt))  # Start from k+1, affected by x1
}

Nshow <- 3
matplot(t, pmat[, 1:Nshow] + matrix(c(rep(0, Nlen/2), rep(1, Nlen/2+1)), byrow = FALSE, nrow = Nlen+1, ncol = Nshow), col = grey(.5, .5), type = 'l', lty = 1, xlab = 't', ylab = expression('X'['t']), las = 1, bty = 'n')

## LL function, with fixed k
LL_euler <- function(mu, k, x, Dt) {
  n <- length(x) - 1 # eventually x will be applied to the funciton by column so it would be a vector not a matrix.
  xmin   <- x[1:n]
  xplus <- x[2:(n + 1)]
  inc <- xplus - xmin
  inckt <-inc[(k+1): n] # inc from k+1 to n, length = n-k
  
  nu_hat <- (sum(inc) + mu*Dt*sum(xplus[1:(n-k)]))/(n*mu*Dt)

  rss1 <- sum(((inc - mu*nu_hat*Dt)^2)[1:k])
  rss2 <- sum((inckt - mu*nu_hat*Dt + mu*Dt*xplus[1:(n-k)])^2)
  sigma2_hat <- (rss1+rss2)/n
  
  logLik1 <- -0.5 * k * log(2 * pi * sigma2_hat * Dt) - rss1^2 / (2 * sigma2_hat * Dt)
  logLik2 <- -0.5 * (n - k) * log(2 * pi * sigma2_hat * Dt) - rss2^2 / (2 * sigma2_hat * Dt)
  logLik <- logLik1 + logLik2
  return(logLik)
}

k_candidates <- 1:5
MLEmu_mat <- matrix(NA, nrow = length(k_candidates), ncol = ncol(pmat))
for (i in 1: length(k_candidates)){
  k <- k_candidates[i]
  MLEmu_mat[i,] <- apply(pmat, 2, function(x){
    optimize(LL_euler, c(0,15), k=k, x = x, Dt = Dt)$minimum
  })
}

