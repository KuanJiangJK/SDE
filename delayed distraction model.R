# simulate driving distraction delayed model
rm(list = ls())

# Model parameters
# dXt = c1 * sin(c2*t + c3) + Yt
# dYt = mu * (Xt-tau - Xt-tau-1)dt + sigma * dW
c1 <- 1 # parameters for xt
c2 <- 1 # parameters for xt
c3 <- 1 # parameters for xt
mu <- 1 # parameters for yt
sigma <- 0.5 # parameters for yt
tau <- 5 # time steps delayed
Dt <- 0.1
T <- 10 # Time up to which paths should be simulated
t <- seq(0, T, Dt)
Npaths <- 1e2 # number of sample paths to simulate
Nlen <- T/Dt
##
pmat_x <- matrix(0, ncol = Npaths, nrow = Nlen+1)
pmat_y <- matrix(0, ncol = Npaths, nrow = Nlen+1)

# generate from the 1st data to the tau+1'th  data
# tau+1'th data point determined by tau'th data point, and the 0th difference of figures
# tau+2'th data point determined by the tau+1'th data point, and the 1st differnce of figures (x1-x0)
for (i in 1: (tau+1)){  
  pmat_y[i+1, ] <- pmat_y[i, ]+ sigma * rnorm(Npaths, 0, sqrt(Dt)) # the observation on ti, the i'th observation
  pmat_x[i+1, ] <- pmat_x[i, ]+ (c1*sin(c2*t[i]+c3) + pmat_y[i,]) * Dt
}
for(j in (tau+2):Nlen){
  pmat_y[j+1, ] <- pmat_y[j, ] + mu*(pmat_x[(j-tau), ] - pmat_x[(j-tau-1), ])*Dt + sigma*rnorm(Npaths, 0, sqrt(Dt))  # the observation on tj, the j'th observation
  pmat_x[j+1, ] <- pmat_x[j, ]+ (c1*sin(c2*t[j]+c3) + pmat_y[j,]) * Dt
}

par(mfrow = c(1, 1))
Nshow <- 3
colors <- c("blue", "red", "green")
matplot(t,
        pmat_x[, 1:Nshow]
        # + matrix(c(rep(0, Nlen/2), rep(1, Nlen/2+1)),
        #                            byrow = FALSE, nrow = Nlen+1, ncol = Nshow)
        ,
        type = 'l', lty = 1, col = colors,
        xlab = 't', ylab = 'X',
        las = 1, bty = 'n', ylim = range(c(pmat_x[,1:Nshow], pmat_y[,1:Nshow])))

# Add second group y
matlines(t,
         pmat_y[, 1:Nshow] 
         # + matrix(c(rep(0, Nlen/2), rep(1, Nlen/2+1)),
         #                            byrow = FALSE, nrow = Nlen+1, ncol = Nshow)
         ,
         lty = 2, col = colors)

legend("bottomright",
       legend = c("x path 1", "x path 2", "x path 3",
                  "y path 1", "y path 1", "y path 1"),
       col = rep(colors, 2),
       lty = c(rep(1, 3), rep(2, 3)),
       cex = 0.5,
       bty = "n")



