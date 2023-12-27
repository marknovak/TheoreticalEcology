library(deSolve)


CAllee <- function(t, y, p) {
  N <- y[1]
  with(as.list(p), {
    dNdt <- r * N * (1 - N / K) * (N / A - 1)
    return(list(dNdt))
  })
}


t <- 1:1500

# Now run your simulation using the ode function of deSolve
out1 <-
  data.frame(ode(
    y = c(N = 35.1),
    times = t,
    func = CAllee,
    parms = c(r = 0.01, K = 100, A = 35)
  ))
out2 <-
  data.frame(ode(
    y = c(N = 34.9),
    times = t,
    func = CAllee,
    parms = c(r = 0.01, K = 100, A = 35)
  ))
out3 <-
  data.frame(ode(
    y = c(N = 120),
    times = t,
    func = CAllee,
    parms = c(r = 0.01, K = 100, A = 35)
  ))

N1 <- out1$N
dNdt1 <- N1[-1] - N1[-length(N1)]
N2 <- out2$N
dNdt2 <- N2[-1] - N2[-length(N2)]
N3 <- out3$N
dNdt3 <- N3[-1] - N3[-length(N3)]

# Put plots side-by-side
par(mfrow = c(2, 2))

plot(out1$time,
     N1,
     type = 'l',
     lwd = 4,
     ylim = c(0, 120))
plot(
  dNdt1 ~ N1[-length(N1)],
  type = 'l',
  lwd = 4,
  xlab = 'N',
  xlim = c(0, 110)
)

plot(out2$time,
     N2,
     type = 'l',
     lwd = 4,
     ylim = c(0, 120))
plot(
  dNdt2 ~ N2[-length(N2)],
  type = 'l',
  lwd = 4,
  xlab = 'N',
  xlim = c(0, 110)
)


par(
  mfrow = c(1, 2),
  cex.lab = 2,
  mgp = c(1.5, 0.6, 0),
  las = 1
)
plot(
  out1$time,
  N1,
  type = 'n',
  lwd = 3,
  ylim = c(0, 120),
  ylab = 'N',
  xlab = 'Time',
  axes = FALSE
)
points(out2$time, N2, type = 'l', lwd = 3)
points(out3$time, N3, type = 'l', lwd = 3)
abline(h = 35, lty = 2, lwd = 1)
axis(1, at = 0, labels = 0)
axis(2, at = c(0, 35, 100), labels = c(0, 'A', 'K'))
box(lwd = 2)

par(las = 1)
plot(
  dNdt1 ~ N1[-length(N1)],
  type = 'l',
  lwd = 3,
  xlab = 'N',
  xlim = c(0, 110),
  ylim = c(-0.1, 0.3),
  ylab = "dN/dt",
  axes = FALSE
)
points(dNdt2 ~ N2[-length(N2)], type = 'l', lwd = 3)
points(dNdt3 ~ N3[-length(N3)], type = 'l', lwd = 3)
abline(h = 0, lwd = 1, lty = 2)
axis(1, at = c(0, 35, 100), labels = c(0, 'A', 'K'))
axis(2, at = 0, labels = 0)
box(lwd = 2)



plot(
  dNdt1 ~ N1[-length(N1)],
  type = 'l',
  lwd = 3,
  xlab = 'N',
  xlim = c(0, 110),
  ylim = c(-0.1, 0.3),
  ylab = "dN/dt",
  axes = FALSE
)
points(dNdt2 ~ N2[-length(N2)], type = 'l', lwd = 3)
points(dNdt3 ~ N3[-length(N3)], type = 'l', lwd = 3)
abline(h = 0, lwd = 1, lty = 2)
axis(1, at = c(0, 35, 100), labels = c(0, 'A', 'K'))
axis(2, at = 0, labels = 0)
box(lwd = 2)

points(
  35,
  0,
  cex = 2,
  lwd = 3,
  pch = 21,
  bg = 'white'
)
points(
  100,
  0,
  cex = 2,
  lwd = 3,
  pch = 21,
  bg = 'black'
)
points(
  0,
  0,
  cex = 2,
  lwd = 3,
  pch = 21,
  bg = 'black'
)