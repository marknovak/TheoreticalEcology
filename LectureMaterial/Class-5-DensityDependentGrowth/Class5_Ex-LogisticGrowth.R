qm <- 2
# Simulating Logistic Growth
########################
# Part A. The discrete logistic
########################
# By simulating the following discrete logistic growth model and plotting the dynamics,
# convince yourself that there are two equilibrium population sizes for logistic growth:
# One in which the population always achieves its carrying capacity regardless of its
# starting population size, and one that is dependent upon the initial population size.

# Define a function for the deterministic discrete logistic whose input parameters
# are N(0), total time T, r_d, and K
DlogisK <- function(N0, T, rd, K) {
  N <- numeric()
  N[1] <- N0
  for (t in 1:T) {
    N[t + 1] <- N[t] + rd * (1 - N[t] / K) * N[t]
  }
  return(N)
}

# Try out various combinations of N0, T, r_d, and K (but make sure to keep r_d less than 1.5)

out <- DlogisK(N0 = 0.01,
               T = 30,
               rd = 0.5,
               K = 100)
quartz()
plot(
  out,
  type = 'o',
  xlab = 'Time',
  ylab = 'N',
  ylim = c(0, 150)
)

# overlay dynamics for different parameter values. For example, N0=1...
out2 <- DlogisK(N0 = 1,
                T = 30,
                rd = 0.5,
                K = 100)
points(out2, type = 'o', col = 'grey')

# What happens when N0 > K?


# What does a plot of Nt+1 - Nt versus Nt look like?
# Start with N[1]=0.01 again.

##################################################################
##################################################################
##################################################################
# Pause here in class
##################################################################
##################################################################
##################################################################

#########################
#########################
# The continuous logistic
#########################
#########################
# Let's clear all the results of our previous simulations using the command:
rm(list = ls()) # clears workspace
# (I put this line of code at the top of all my scripts to avoid mistakes.)

# If it's not already installed on your computer, install the R-package called 'deSolve'.
# This only needs to be done once (unless you want to update the installed package).

# install.packages('deSolve')

# To use the installed package you first have to 'load' it into the library.
library(deSolve)

# Then define the function for the continuous logistic
# Notice that this function has three inputs:  
# the starting population size'N0',
# the total time 'T', 
# and a third input variable 'p' representing the two parameters
# of our equation: r and a.

ClogisK <- function(t, y, p) {
  N <- y[1]
  with(as.list(p), {
    dNdt <- r * N * (1 - N / K)
    return(list(dNdt))
  })
}

# We therefore have to include both r and a in the params variable:
params <- c(r = 0.1, K = 100)

# Then specify the starting abundance and the time over which we want to simulate dynamics:
N0 <- c(N = 1)
t <- 1:100

# Now run your simulation using the ode function of deSolve
out <- ode(
  y = N0,
  times = t,
  func = ClogisK,
  parms = params
)

# Convert this output to a data.frame to ease plotting and data extraction
out <- data.frame(out)
quartz() # Windows users may need to comment-out this line
plot(out$time, out$N)

# Plot the change in N between time t and t+1 as a function N(t)
# To do this we must exclude the last observed N
N <- out$N
dNdt <- N[-1] - N[-length(N)]
quartz()
plot(dNdt ~ N[-length(N)])

# Put these last two plots side-by-side
quartz(height = 4 * qm, width = 7 * qm)
par(mfrow = c(1, 2))
plot(out$time, out$N, type = 'l', lwd = 4)
plot(dNdt ~ N[-length(N)],
     type = 'l',
     lwd = 4,
     xlab = 'N')

################
# Try out various combinations of N0, T, r, and alpha, saving the results into new variables,
# out2, N2, and dNdt2, for example, to avoid overwritting your previous simulations.

# Different per capita growth rate
out2 <-
  data.frame(ode(
    y = c(N = 1),
    times = t,
    func = ClogisK,
    parms = c(r = 0.2, K = 100)
  ))
N2 <- out2$N
dNdt2 <- N2[-1] - N2[-length(N2)]

# Different carrying capacity
out3 <-
  data.frame(ode(
    y = c(N = 1),
    times = t,
    func = ClogisK,
    parms = c(r = 0.2, K = 110)
  ))
N3 <- out3$N
dNdt3 <- N3[-1] - N3[-length(N3)]

# Plot the results of out, out2, and out3 all together on the same graph.
# Note:  You may need to adjust the minimum and maximum values of xlim and ylim
# to see all your data.
quartz(height = 4 * qm, width = 7 * qm)
par(mfrow = c(1, 2))
plot(
  N ~ time,
  data = out,
  type = 'l',
  lwd = 4,
  ylim = c(0, 150)
)
points(
  N ~ time,
  data = out2,
  type = 'l',
  lwd = 4,
  col = 'red'
)
points(
  N ~ time,
  data = out3,
  type = 'l',
  lwd = 4,
  col = 'blue'
)

plot(
  dNdt ~ N[-length(N)],
  type = 'l',
  lwd = 4,
  xlab = 'N',
  ylim = c(-2, 6),
  xlim = c(0, 150)
)
points(dNdt2 ~ N2[-length(N2)],
       type = 'l',
       lwd = 4,
       col = 'red')
points(dNdt3 ~ N3[-length(N3)],
       type = 'l',
       lwd = 4,
       col = 'blue')
abline(0, 0, col = 'grey')

##########
# How do the values of r, K, and N0 affect the maximum value of dNdt?

#########################################################################################
#########################################################################################
#########################################################################################
