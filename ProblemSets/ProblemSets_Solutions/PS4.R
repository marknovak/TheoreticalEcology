#################################################################
######### Theoretical Ecology - Problem set  ####################
#################################################################
# In[1]:= du1dtau = u1*(1 - u1 - a12*u2)
# du2dtau = \[Rho] * u2*(1 - a21*u1 - u2)
#
# Out[1]= u1 (1 - u1 - a12 u2)
#
# Out[2]= (1 - a21 u1 - u2) u2 \[Rho]
#
# In[3]:= isou2 = Solve[du1dtau == 0, u2]
#
# Out[3]= {{u2 -> (1 - u1)/a12}}
#
# In[4]:= isou1 = Solve[du2dtau == 0, u1]
#
# Out[4]= {{u1 -> (1 - u2)/a21}}
#
# In[13]:= u1 /. isou1
#
# Out[13]= {(1 - u2)/a21}
#
# In[16]:= Solve[(u1 /. isou1) == u1, u2]
#
# Out[16]= {{u2 -> 1 - a21 u1}}

##################################################################################
##################################################################################
rm(list = ls())

source('VectorField.R')
library(deSolve)

isou2 <- function(x) {
  (1 - x) / a12
}
isou1 <- function(x) {
  1 - a21 * x
}


rho = 1 # species have equal density-independent growth rates

par(
  mfrow = c(2, 2),
  mar = c(1,1,1,1),
  pty = 's',
  xaxs = 'i',
  yaxs = 'i',
  tcl = 0.1,
  mgp = c(1, 0.1, 0),
  cex.axis = 0.5
)
# species 1 dominant
a21 = 1.5
a12 = 0.5

plotVectorField(function(u1, u2) {
  c(u1 * (1 - u1 - a12 * u2), 
    rho * u2 * (1 - u2 - a21 * u1))
}, c(0, 2), c(0, 2))
curve(isou2, add = T)
curve(isou1, add = T, lty = 2)
legend(
  'topright',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title(main = 'Dominance by u1', xlab = 'u1', ylab = 'u2')

# species 2 dominant
a21 = 0.5
a12 = 1.5

plotVectorField(function(u1, u2) {
  c(u1 * (1 - u1 - a12 * u2), 
    rho * u2 * (1 - u2 - a21 * u1))
}, c(0, 2), c(0, 2))
curve(isou2, add = T)
curve(isou1, add = T, lty = 2)
legend(
  'topright',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title(main = 'Dominance by u2', xlab = 'u1', ylab = 'u2')

# coexistence
a21 = 0.5
a12 = 0.7

plotVectorField(function(u1, u2) {
  c(u1 * (1 - u1 - a12 * u2), 
    rho * u2 * (1 - u2 - a21 * u1))
}, c(0, 2), c(0, 2))
curve(isou2, add = T)
curve(isou1, add = T, lty = 2)
legend(
  'topright',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title(main = 'Coexistence', xlab = 'u1', ylab = 'u2')

# priority effect
a21 = 1.5
a12 = 1.7

plotVectorField(function(u1, u2) {
  c(u1 * (1 - u1 - a12 * u2), 
    rho * u2 * (1 - u2 - a21 * u1))
}, c(0, 2), c(0, 2))
curve(isou2, add = T)
curve(isou1, add = T, lty = 2)
legend(
  'topright',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title(main = 'Priority effect', xlab = 'u1', ylab = 'u2')

##################################################################################
model <- function(t, x, params) {
  u1 <- x[1]
  
  u2 <- x[2]
  
  with(as.list(parameters), {
    du1dt <- u1 * (1 - u1 - a12 * u2)
    du2dt <- rho * u2 * (1 - u2 - a21 * u1)
    res <- c(du1dt, du2dt)
    list(res)
  })
}

T <- 20
t <- seq(0, T, by = 1)
xstart <- c(u1 = 0.1, u2 = 0.1)

par(mfrow = c(3, 2))

# species 1 dominant
parameters <- c(a21 = 1.5, a12 = 0.5, rho = 1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1 ~ time,
     data = out,
     type = 'l',
     ylab = 'u',
     ylim = c(0, 1.5))
lines(out$time, out$u2, lty = 2)
legend(
  'topleft',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title('Dominance by u1')

# species 2 dominant
parameters <- c(a21 = 0.5, a12 = 1.5, rho = 1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1 ~ time,
     data = out,
     type = 'l',
     ylab = 'u',
     ylim = c(0, 1.5))
lines(out$time, out$u2, lty = 2)
legend(
  'topleft',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title('Dominance by u2')

# coexistence
parameters <- c(a21 = 0.5, a12 = 0.7, rho = 1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1 ~ time,
     data = out,
     type = 'l',
     ylab = 'u',
     ylim = c(0, 1.5))
lines(out$time, out$u2, lty = 2)
legend(
  'topleft',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title('Coexistence')

# priority effect - starting abundance 1
xstart <- c(u1 = 0.1, u2 = 0.1)
parameters <- c(a21 = 1.5, a12 = 1.7, rho = 1)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1 ~ time,
     data = out,
     type = 'l',
     ylab = 'u',
     ylim = c(0, 1.5))
lines(out$time, out$u2, lty = 2)
legend(
  'topleft',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title('Priority effect')

# priority effect - starting abundance 2
xstart <- c(u1 = 0.5, u2 = 0.1)
parameters <- c(a21 = 1.5, a12 = 1.7, rho = 2)
out <- as.data.frame(ode(xstart, t, model, parameters))
plot(u1 ~ time,
     data = out,
     type = 'l',
     ylab = 'u',
     ylim = c(0, 1.5))
lines(out$time, out$u2, lty = 2)
legend(
  'topleft',
  legend = c('u1', 'u2'),
  lty = c(1, 2),
  bg = 'white'
)
title('Priority effect')

##################################################################################
##################################################################################
##################################################################################
