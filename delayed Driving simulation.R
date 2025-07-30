rm(list = ls())
# Model parameters
c1 <- 1 
c2 <- 1
c3 <- 1
mu <- -2 # 
sigma <- 1 # diffusion
tau <- 8 # discrete delay in time steps
Dt <- 1/60
T <- 120
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
  pmat_y[j + 1, ] <- pmat_y[j, ] + mu * delay_term * Dt + noise_term
}

# Plotting first 3 paths
par(mfrow = c(1, 1))
Nshow <- 3
colors <- c("blue", "red", "green")
matplot(t,
        pmat_y[, 1:Nshow],
        type = 'l', lty = 1, col = colors,
        xlab = 't', ylab = 'Y',
        las = 1, bty = 'n',
        ylim = c(-3,3))

legend("bottomright",
       legend = c("y path 1", "y path 2", "y path 3"),
       col = colors,
       lty = 1,
       cex = 0.5,
       bty = "n")

#################### MLE estimator by Udo ########
MLE <- function(tau, x){
  f <- c1 * sin(c2 * (t[1:(length(t)-tau-1)]) + c3)
  xmin <- x[tau:(Nlen-1)]
  xplus <- x[(tau+1):Nlen]
  xdel <- x[1:(Nlen - tau)]
  muhat <- sum((xplus-xmin)*(xdel + f)) / (sum((xdel + f)^2)*Dt)
  sigmasqhat <- 1/((Nlen-tau)*Dt)*sum((xplus-xmin-muhat*(xdel+f)*Dt)^2)
  return(c(muhat, sigmasqhat)) #modify: xdel-f to xdel+f
}

LL <- function(x, theta, tau){
  mu <- theta[1]
  sigmasq <- theta[2]
  
  idx <- (tau+2):(length(x)-1)
  if(length(idx) > 1){
    LL <- sum(sapply(idx, foo, mu, sigmasq, x, tau))
  }else{
    LL <- foo(idx, foo, mu, sigmasq, x, tau)
  }
  return(LL)
}

foo <- function(idx, mu, sigmasq, x, tau){
  ll = dnorm(x[idx+1], mean = x[idx] + mu*(x[idx-tau-1] + sin(c2 * t[idx-tau-1] + c3))*Dt, sd = sqrt(sigmasq*Dt), log = TRUE)
  return(ll) # modification made here idx-tau -> idx-tau-1
}

# For a single path
# LLs <- numeric(20)
# x <- pmat_y[, 2]
# for(i in 1:20){
#   theta <- MLE(i, x)
#   LLs[i] <- LL(x, theta, i)
# }
# plot(1:20, LLs, xlab = 'tau', ylab = 'LL', type = 'l')
theta <- list()
LLs <- list()
path_est <- function(x, i){
  theta <- MLE(i, x)
  LLs <- LL(x, theta, i)
  solution <- data.frame("mu" = theta[1],"sigma2" = theta[2], "likelihood"= LLs, "tau" = i)
}
# for one candidate i # test
# i = 1 # test
# total_est <- apply(pmat_y, MARGIN = 2, function(x) path_est(x, i)) # estimate for only 1 candidate i
# total_est <- do.call(rbind, apply(pmat_y, MARGIN = 2, function(x) path_est(x, i))) # estimate for only 1 candidate i
total_est_tau <- lapply(1:20, function(i) do.call(rbind, apply(pmat_y, MARGIN = 2, function(x) path_est(x, i))))
total_est_tau <- do.call(rbind, total_est_tau)
total_est_tau$path <- rep(1:Npaths, 20)
path_tau <- lapply(split(total_est_tau, total_est_tau$path), function(x) x$tau[which.max(x$likelihood)])
table(do.call(rbind, path_tau))
est_tau8 <- subset(total_est_tau, tau == 8)
median(est_tau8$mu) # -1.945386
median(est_tau8$sigma2) # 1.005818
quantile(est_tau8$mu, c(0.025, 0.975))
#     2.5%     97.5% 
#   -2.151817 -1.727323 
quantile(est_tau8$sigma2, c(0.025, 0.975))
#     2.5%     97.5% 
#   0.9725701 1.0373283

library(ggplot2)
ggplot(total_est_tau, aes(x = factor(tau), y = likelihood))+
  geom_boxplot() +
  labs(x = "candidate tau", y = "Log-likelihood")
ggplot(total_est_tau, aes(x = factor(tau), y = mu)) +
  geom_boxplot() +
  labs(x = "Candidate tau", y = "Estimated mu")
ggplot(total_est_tau, aes(x = factor(tau), y = sigma2)) +
  geom_boxplot() +
  labs(x = "Candidate tau", y = "Estimated sigma2")



