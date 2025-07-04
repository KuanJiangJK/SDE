# Udo's code (annotated). add Euler likelihood approximation part (fix mu, profiling over nu & sigma2).
# Model parameters for OU
mu <- 1.5
nu <- .6
sigma <- 1.5

# Simulation parameters
Dt <- 0.01
T <- 40 # Time up to which paths should be simulated
t <- seq(0, T, Dt)
Npaths <- 1e2 # number of sample paths to simulate
Nlen <- T/Dt

# Path simulation with true transition density
set.seed(123)
pmat <- matrix(0, ncol = Npaths, nrow = Nlen+1)
for(i in 1:Nlen){
  pmat[i+1, ] <- pmat[i, ]*exp(-mu*Dt) + nu*(1-exp(-mu*Dt)) + rnorm(Npaths, 0, sigma/(sqrt(2*mu))*sqrt(1-exp(-2*mu*Dt)))  # generate paths matrix of OU-process
}

# Visualise the first Nshow sample paths
Nshow <- 3
matplot(t, pmat[, 1:Nshow] + matrix(c(rep(0, Nlen/2), rep(1, Nlen/2+1)), byrow = FALSE, nrow = Nlen+1, ncol = Nshow), col = grey(.5, .5), type = 'l', lty = 1, xlab = 't', ylab = expression('X'['t']), las = 1, bty = 'n')


LL <- function(mu, x, Dt){ # mu is reversion rate, x is a vector of data, dt is time increment.
  n <- length(x)-1
  xplus <- x[2:(n+1)] # x without the first point
  xmin <- x[1:n] # x without the last point
  if(length(mu) > 1){ # if mu is a vector, then compute the log-likelihood for each value in the vector using sapply, this enables vectorized likelihood evaluations
      logLik <- sapply(mu, LL, x, Dt) # sapply means simplified form of lapply(). useful for list object. 
      # this sapply applies the function LL recursively to each value in mu, passing along the same x and Dt. logLik[i] <- LL(mu[i], x, Dt) for every i. The return is a vector of log-likelihoods, for each mu's respectively.
      # sapply() applies the function LL() to each element of the vector mu, keeping x and Dt fixed!!!!
  }else{ # single value case.
    beta <- exp(-Dt*mu) # component of the transition density (let mu be fixed). thus the mean part of the transition denstiy becomes beta*x + nu*(1-beta), variance becomes sigma^2/2*mu*(1-beta^2)
    nu <- 1/n*sum(xplus-xmin*beta)/(1-beta) # nu expressed in terms of known mu (transformed to be beta) (profiling likelihood)
    sigma2 <- 1/n * 2*mu/(1-beta^2) * sum((xplus - xmin * beta + nu*(1-beta))^2) # sigma expressed in terms of known mu (transformed to be beta), and sigma2 (profiling likelihood)
    logLik <- -sum(dnorm(xplus-xmin*beta-nu*(1-beta), mean = 0, sd = sqrt(.5*sigma2/mu*(1-beta^2)), log = TRUE)) # so here the likelihood becomes a function of mu, to be optimized
  }
  return(logLik)
}

LL_euler <- function(mu, x, Dt) {
  n <- length(x) - 1
  xmin   <- x[1:n]
  xplus <- x[2:(n + 1)]
  
  if (length(mu) > 1) {
    return(sapply(mu, LL_euler, x = x, Dt = Dt))
  } else {
    nu_hat <- mean(xmin) + mean((xplus - xmin) / Dt) / mu
    x_exp <- xmin - mu * (xmin - nu_hat) * Dt  
    residuals <- xplus - x_exp
    sigma2_hat <- sum(residuals^2) / (n * Dt)
    logLik <- -0.5 * n * log(2 * pi * sigma2_hat * Dt) - sum(residuals^2) / (2 * sigma2_hat * Dt)
    return(logLik)
  }
}


# Subsample pmat and obtain ML estimates
KSubSteps <- 50 # Maximum of multiples of Delta t for which data are subsampled
MLEmu <- MLEnu <- MLEsig <- matrix(0, nrow = ncol(pmat), ncol = KSubSteps) #Estimated paras for each path at each subsampling level
MLEmu_euler  <- MLEnu_euler <- MLEsig_euler <- matrix(0, nrow = ncol(pmat), ncol = KSubSteps)


for(k in 1:KSubSteps){
  x <- pmat[seq(1,nrow(pmat), k), ]
  n <- nrow(x)
  Dt_k <- k*Dt
  # true MLE
  MLEmu[,k] <- apply(x, 2, 
                     function(x, Dt){
                       optimize(LL, c(0, 15), x, Dt)$minimum}, # c(0, 15) is the input for x (here mu)
                     # $minimum: the x value that minimizes (or maximizes); $objective: the negative log-likelihood value at that x value
                     Dt_k) # apply(X, MARGIN, FUN, ...), Dt here serves as the argument for the second inline function
  betaMat <- matrix(exp(-Dt_k*MLEmu[,k]), byrow = TRUE, ncol = ncol(x), nrow = nrow(x)-1)
  MLEnu[,k] <- 1/(n-1)*colSums(x[2:n, ] - x[1:(n-1), ]*betaMat)/(1-exp(-Dt_k*MLEmu[,k]))
  MLEsig[,k] <- 1/(n-1) * 2*MLEmu[,k]/(1-exp(-2*Dt_k*MLEmu[,k])) * colSums((x[2:n, ] - x[1:(n-1), ]*betaMat + MLEnu[,k]*(1-betaMat))^2)
  
  # Euler MLE
  MLEmu_euler[, k] <- apply(x, 2, 
                            function(x, Dt) {
                              optimize(LL_euler, interval = c(0.01, 15), x, Dt)$minimum}, Dt_k)
  MLEnu_euler[, k] <- colMeans(x[1:(n-1),]) + colMeans(x[2:n,]-x[1:(n-1),])/ MLEmu_euler[,k]
  rep_MLEmu <- matrix(MLEmu_euler[,k], byrow = TRUE, nrow = (nrow(x)-1), ncol = length(MLEnu[,k])) # this line replicates a vector into a matrix with same rows. So that I can do minus between matrix. (said to be an efficient way, from Stackoverflow)
  rep_MLEnu <- matrix(MLEnu_euler[,k], byrow = TRUE, nrow = (nrow(x)-1), ncol = length(MLEmu[,k])) # this line replicates a vector into a matrix with same rows. So that I can do minus between matrix
  x_exp <- x[1:(n-1),] - rep_MLEmu * (x[1:(n-1),] - rep_MLEnu)*Dt_k
  residuals <- x[2:n,] - x_exp
  MLEsig_euler[, k] <- colSums(residuals^2) / ((n - 1) * Dt_k)
}

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
boxplot(MLEmu, main = "Exact MLE mu", las = 1, ylim = range(c(MLEmu, MLEmu_euler)))
abline(h = mu, col = 2)
boxplot(MLEmu_euler, main = "Euler MLE mu", las = 1, ylim = range(c(MLEmu, MLEmu_euler)))
abline(h = mu, col = 2)

boxplot(MLEnu, main = "Exact MLE nu", las = 1, ylim = range(c(MLEnu, MLEnu_euler)))
abline(h = nu, col = 2)
boxplot(MLEnu_euler, main = "Euler MLE nu", las = 1, ylim = range(c(MLEnu, MLEnu_euler)))
abline(h = nu, col = 2)

boxplot(MLEsig, main = "Exact MLE sigma^2", las = 1, ylim = range(c(MLEsig, MLEsig_euler)))
abline(h = sigma^2, col = 2)
boxplot(MLEsig_euler, main = "Euler MLE sigma^2", las = 1, ylim = range(c(MLEsig, MLEsig_euler)))
abline(h = sigma^2, col = 2)

# Plot true MLEs as function of (subsampled) time step size
par(mfrow=c(3,1))
par(mar = c(4, 4, 2, 1))
matplot(t(MLEmu), type = 'l', col = grey(.5, .5), las = 1, lty = 1)
matplot(t(MLEnu), type = 'l', col = grey(.5, .5), las = 1, lty = 1)
matplot(t(MLEsig), type = 'l', col = grey(.5, .5), las = 1, lty = 1)

