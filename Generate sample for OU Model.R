set.seed(123)
# Model parameters for constant drift model
# dXt = v dt + sigma dwt

v <- .5 # Constant drift
sigma <- 1 # diffusion

# Model parameters for OU
mu <- .5
nu <- 1
sigma <- 1

# Simulation parameters
Dt <- .01 # time step, delta t
T <- 4 # Time up to which paths should be simulated
t <- seq(0, T, Dt) # time vector, from 0 to 4, step = 0.01
Npaths <- 1e5 # number of sample paths to simulate
Nlen <- T/Dt # number of time steps

# Step 1 Generate a matrix with Gaussian RVs that represent the diffusion part of increments: sigma * (W_{t+Delta t} - W_{t})

# pmat <- matrix(rnorm(Npaths*Nlen, 0, sigma), ncol = Npaths, nrow = Nlen) # better to add the scaling? according to Bronian motion?
pmat <- matrix(sqrt(Dt) * rnorm(Npaths*Nlen, 0, sigma), ncol = Npaths, nrow = Nlen)

# Step 2 Add constant offsets that represent the drift part of increments: v dt
pmat <- pmat + v*Dt

# Step 3 Add initial value
pmat <- rbind(rep(0, Npaths), pmat)
# Step 4 Compute paths by adding increments
pmat <- apply(pmat, 2, cumsum) # 2 means by columns
# cumsum: cumulative sum. Xt = X0 + SIGMA(from 1 to n)(v*dt + sigma*N(0,1)), full trajectories

# Note: To implement the OU process you should modify Step 2 and integrate it with Step 4, simply using cumsum won't work anymore

# Visualise the first Nshow sample paths
Nshow <- 1e2
matplot(t, pmat[, 1:Nshow], col = grey(.5, .5), type = 'l', lty = 1, xlab = 't', ylab = expression('X'['t']), las = 1, bty = 'n')


#####################
# Simulate for OU
# dXt = nu*(mu - Xt) dt + sigma dwt
# Model parameters for OU
mu <- .5 # long-term mean
nu <- 1 # mean reversion speed
sigma <- 1 # diffusion

# Simulation parameters
Dt <- .01 # time step, delta t
T <- 4 # Time up to which paths should be simulated
t <- seq(0, T, Dt) # time vector, from 0 to 4, step = 0.01
Npaths <- 1e5 # number of sample paths to simulate
Nlen <- T/Dt # number of time steps

# Step 1 Generate a matrix with Gaussian RVs that represent the diffusion part of increments: sigma * (W_{t+Delta t} - W_{t})
# diffusion <- matrix(rnorm(Npaths*Nlen, 0, sigma), ncol = Npaths, nrow = Nlen) 
diffusion <- matrix(sqrt(Dt) * rnorm(Npaths*Nlen, 0, sigma), ncol = Npaths, nrow = Nlen)

# Step 2 Add constant offsets that represent the drift part of increments: v dt
pmat <- matrix(0, nrow = Nlen + 1, ncol = Npaths)  #starting values, all 0s, Nlen + 1 means one extra time point (initial 0).
for (i in 1:Nlen){ #For all pathes, accumulate datapoints based on the present value.
  pmat[i+1, ] = pmat[i, ] + (nu * (mu - pmat[i,])) * Dt + diffusion[i,]
  print(i)
}
# Padding the difffusion with one extra zero (so I can add them up)

# Visualise the first Nshow sample paths
Nshow <- 3
matplot(t, pmat[, 1:Nshow], col = grey(.5, .5), type = 'l', lty = 1, xlab = 't', ylab = expression('X'['t']), las = 1, bty = 'n')






