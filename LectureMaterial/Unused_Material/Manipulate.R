library(deSolve)
library(manipulate)

# Wrap simulation and plots into function to pass to manipulate
myPlot <- function(r.val, K.val, a.val, 
                   h.val, e.val, m.val){
  
# Define model
LVmod <- function(times, init, parms) {
  with(as.list(c(init, parms)), {
    dR <- r * R * (1 - R/K) - a * R * C / (1 + a * h * R)
    dC <- e *  a * R * C / (1 + a * h * R) - m * C
    return(list(c(dR, dC)))
  })
}

parms <- c(
          r = r.val,
          K = K.val,
          a = a.val, 
          h = h.val,
          e = e.val,
          m = m.val
          )

init <- c(R = 0.1,
          C = 0.1)

times <- seq(0, 800, by = 1)

out <- as.data.frame(ode(func = LVmod,
                         y = init,
                         parms = parms,
                         times = times))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate instantaneous fraction feeding
fracFeed <- function(R, parms){
  with(as.list(parms), {
       fF = h * a * R / (1 + a * h * R)
    return(fF)
  })
}

out$fF <- fracFeed(out$R, parms)

# Expected fraction feeding at steady-state
fFss <- with(as.list(parms), {
  return(m * h / e)
  })

# Analytical boundaries between cycles, fixed point, and extinction
fFbounds <- with(as.list(parms), {
  fFcyc = 1 - m / (a * e * K) - 1 / (a * h * K)
  fFext = 1 - m / (a * e * K) 
  return(c(fFcyc = fFcyc, fFext = fFext))
  })

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2, 1),
    las = 1,
    mar = c(3.5, 3, 0.5, 0.5),
    mgp = c(1.5, 0.2, 0),
    tcl = -0.1,
    cex.lab = 1.25)

matplot(out[, c('R', 'C')],
        type = 'l',
        xlab = 'Time',
        ylab = 'Abundance')

legend('topright',
       c('Resource', 'Consumer'),
       lty = c(1, 2),
       col = c(1, 2),
       # bty = 'n',
       cex = 0.8,
       horiz = FALSE)
       
plot(out$time,
     out$fF,
     ylim = c(0, 1),
     type = 'l',
     xlab = 'Time',
     ylab = 'Fraction feeding')
abline(h = fFss,
       col = 'blue',
       lty = 2)
abline(h = fFbounds['fFcyc'],
       col = 'red',
       lty = 4)
abline(h = fFbounds['fFext'],
       col = 'orange',
       lty = 4)
legend('topright',
       c('Instantaneous', 'Steady state', 
       'Cycle boundary','Extinction boundary'),
       lty = c(1, 2, 4, 4),
       col = c('black', 'blue', 'red', 'orange'),
       cex = 0.8)

} # end of plotting function


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

manipulate(
  myPlot(r.val, K.val, a.val, h.val, e.val, m.val),
  r.val = slider(0, 0.5, initial = 0.2, step = 0.05),
  K.val = slider(0, 100, initial = 10, step = 1),
  a.val = slider(0, 0.5, initial = 0.15, step = 0.05),
  h.val = slider(0, 2, initial = 1, step = 0.1),
  e.val = slider(0, 1, initial = 0.2, step = 0.05),
  m.val = slider(0, 0.2, initial = 0.05, step = 0.01)
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
