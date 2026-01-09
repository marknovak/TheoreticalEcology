rm(list = ls()) # clears workspace
quartz.options(pointsize = 12 * qm)
###########################################################
###########################################################
# Geometric growth - Arithmetic (natural) vs. Log-scale
###########################################################
###########################################################

# N vs. t
quartz(height = 8, width = 20)
lab <- c(1, 2, 4, 8, 16, 32)

lambda <-  2
N <- numeric()
N[1] <- 1
T <- 1:5
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}
par(mfrow = c(1, 2))
par(
  mar = c(5, 5, 1.5, 1),
  lwd = 2,
  cex.lab = 1.2,
  mgp = c(2, 0.5, 0),
  tcl = -0.3
)
#---- Natural scale
plot(
  N ~ x,
  xlab = 'Time Steps',
  ylab = 'Population size - N(t)',
  type = 'p',
  pch = 19,
  axes = F
)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 1,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 1,
         lwd = 1,
         col = 'grey')
points(x, N)
box(lwd = 2)
axis(1, at = x, label = x - 1)
axis(2, at = lab, labels = lab)

#---- Log-scale
plot(
  N ~ x,
  log = 'y',
  xlab = 'Time Steps',
  ylab = 'log N(t)',
  pch = 19,
  axes = F
)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 1,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 1,
         lwd = 1,
         col = 'grey')
box(lwd = 2)
axis(1, at = x, label = x - 1)
axis(2, at = lab, labels = lab)

# N vs. t;  Different N(0)
quartz(height = 8, width = 20)
lab <- c(1, 2, 4, 8, 16, 32, 64, 128, 256)
par(mfrow = c(1, 2))
par(
  mar = c(5, 5, 1.5, 1),
  lwd = 2,
  cex.lab = 1.2,
  mgp = c(2, 0.5, 0),
  tcl = -0.3
)
lambda <-  2
N <- numeric()
N[1] <- 1
T <- 1:7
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}

#---- Natural scale
plot(
  N ~ x,
  xlab = 'Time Steps',
  ylab = 'Population size - N(t)',
  pch = 19,
  axes = F,
  ylim = c(1, max(lab))
)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 1,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 1,
         lwd = 1,
         col = 'grey')
box(lwd = 2)
axis(1, at = x, label = x - 1)
axis(2, at = lab, labels = lab)
lambda <-  2
N <- numeric()
N[1] <- 2
T <- 1:7
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}
points(N ~ x,
       pch = 18)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 2,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 2,
         lwd = 1,
         col = 'grey')

#---- log scale
lambda <- 2
N <- numeric()
N[1] <- 1
T <- 1:7
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}
plot(
  N ~ x,
  log = 'y',
  xlab = 'Time Steps',
  ylab = 'log N(t)',
  pch = 19,
  axes = F,
  ylim = c(1, max(lab))
)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 2,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 2,
         lwd = 1,
         col = 'grey')
box(lwd = 2)
axis(1, at = x, label = x - 1)
axis(2, at = lab, labels = lab)
lambda <-  2
N <-  numeric()
N[1] <- 2
T <- 1:7
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}
points(N ~ x, pch = 18)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 2,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 2,
         lwd = 1,
         col = 'grey')

# lambda vs N
quartz(height = 8, width = 20)
par(mfrow = c(1, 2))
par(
  mar = c(5, 5, 1.5, 1),
  lwd = 2,
  cex.lab = 1.2,
  mgp = c(2, 0.5, 0),
  tcl = -0.3
)
plot(
  N[-1] / N[-length(N)] ~ N[-length(N)],
  type = 'b',
  pch = 19,
  ylab = 'Per capita growth rate',
  xlab = 'N(t)',
  ylim = c(1, 3),
  axes = FALSE
)
axis(2)
axis(1, at = N)
box(lwd = 1)
plot(
  N[-1] / N[-length(N)] ~ N[-length(N)],
  log = 'x',
  type = 'b',
  pch = 19,
  ylab = 'Per capita growth rate',
  xlab = 'log N(t)',
  ylim = c(1, 3),
  axes = FALSE
)
axis(2)
axis(1, at = N)
box(lwd = 1)

# # delta N vs N
# quartz(height=8,width=20)
# par(mfrow=c(1,2));par(mar=c(5,5,1.5,1),lwd=2,cex.lab=1.2,mgp=c(2,0.5,0),tcl=-0.3)
# plot(N[-1]-N[-length(N)]~N[-length(N)],type='b',pch=19,ylab='Change in popn size',xlab='N(t)',ylim=c(1,150),axes=FALSE)
# axis(2); axis(1,at=N); box(lwd=1)
# plot(N[-1]-N[-length(N)]~N[-length(N)],log='x',type='b',pch=19,ylab='Change in popn size',xlab='log N(t)',ylim=c(1,150),axes=FALSE)
# axis(2); axis(1,at=N); box(lwd=1)

####################################################################################
####################################################################################
# Numerical approximation of Euler's constant
####################################################################################
####################################################################################
n <- seq(0, 5, 1)
N0 <- 1
rd <- 1

N1 <- N0 * (1 + rd / n) ^ n
e <- exp(rd)
print(e)

est.e <- max(N1 / N0)
print(est.e)

quartz(height = 8, width = 20)
par(mar = c(5, 5, 1, 1),
    lwd = 2,
    cex.lab = 1)
par(mar = c(5, 5, 1, 1),
    lwd = 2,
    cex = 1)
plot(
  n,
  N1 / N0,
  ylab = expression(N[1] / N[0] == N[1] / 1),
  type = 'b',
  pch = 19,
  ylim = c(1, 3)
)
abline(h = e,
       col = 'red',
       lwd = 3,
       lty = 2)

text(max(n) / 2, 2, paste('For n = ', max(n)))
text(max(n) / 2, 1.7,
     bquote((1 + frac("r"["d"], "n")) ^ "n" == .(round(N1[length(N1)] / N0, 6))))
text(max(n) / 2, 1.2,
     bquote('Compare: when r = 1, e' ^ 'r' == .(round(e, 6))))

##########
n <- seq(0, 20, 1)
N1 <- N0 * (1 + rd / n) ^ n

est.e <- max(N1 / N0)

plot(
  n,
  N1 / N0,
  ylab = expression(N[1] / N[0] == N[1] / 1),
  type = 'b',
  pch = 19,
  ylim = c(1, 3)
)
abline(h = e,
       col = 'red',
       lwd = 3,
       lty = 2)

text(max(n) / 2, 2, paste('For n = ', max(n)))
text(max(n) / 2, 1.7,
     bquote((1 + frac("r"["d"], "n")) ^ "n" == .(round(N1[length(N1)] / N0, 6))))
text(max(n) / 2, 1.2,
     bquote('Compare: when r = 1, e' ^ 'r' == .(round(e, 6))))

#########
n <- seq(0, 200, 1)
N1 <- N0 * (1 + rd / n) ^ n

est.e <- max(N1 / N0)

plot(
  n,
  N1 / N0,
  ylab = expression(N[1] / N[0] == N[1] / 1),
  type = 'b',
  pch = 19,
  ylim = c(1, 3)
)
abline(h = e,
       col = 'red',
       lwd = 3,
       lty = 2)
text(max(n) / 2, 2, paste('For n = ', max(n)))
text(max(n) / 2, 1.7,
     bquote((1 + frac("r"["d"], "n")) ^ "n" == .(round(N1[length(N1)] / N0, 6))))
text(max(n) / 2, 1.2,
     bquote('Compare: when r = 1, e' ^ 'r' == .(round(e, 6))))

########################################################################
########################################################################
########################################################################

lab <- c(1, 2, 4, 8, 16, 32)
lambda <-  2
N <- numeric()
N[1] <- 1
T <- 1:5
x <- c(T, max(T) + 1)
for (t in T) {
  N[t + 1] <- lambda * N[t]
}
par(
  mar = c(5, 5, 1.5, 1),
  lwd = 2,
  cex.lab = 1.2,
  mgp = c(2, 0.5, 0),
  tcl = -0.3
)
#---- Natural scale
plot(
  N ~ x,
  xlab = 'Time Steps',
  ylab = 'Population size - N(t)',
  type = 'p',
  pch = 19,
  axes = F
)
segments(x[-length(x)], N[-length(N)], 
         x[-1], N[-length(N)], 
         lty = 1,
         lwd = 1,
         col = 'grey')
segments(x[-1], N[-length(x)], 
         x[-1], N[-1], 
         lty = 1,
         lwd = 1,
         col = 'grey')
points(x, N)
box(lwd = 2)
axis(1, at = x, label = x - 1)
axis(2, at = lab, labels = lab)

r <- log(lambda)
curve(N[1] * exp(r * (x-1)),
      add = TRUE)

message(paste0('lambda = ', lambda, 
               ', r_d = ', lambda-1, 
               ', r = ', round(r,4)))

        